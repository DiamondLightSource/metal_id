import logging
import shutil
from dataclasses import dataclass
from pathlib import Path

import gemmi
import numpy as np


@dataclass
class PeakData:
    density: float
    rmsd: float
    xyz: tuple[float, float, float]


def read_pha(pha_file: Path) -> np.ndarray:
    """Read an ANODE .pha file (PHS format: h k l F FOM phi sigF) and return an array
    with columns h, k, l, F*FOM, phi (degrees), as expected by Mtz.set_data.

    The indices are read from fixed width columns (3I4) as large negative indices can
    run into one another, the remaining columns are whitespace separated.
    """
    rows = []
    with open(pha_file) as fh:
        for line in fh:
            if not line.strip():
                continue
            hkl = (int(line[0:4]), int(line[4:8]), int(line[8:12]))
            f, fom, phi = map(float, line[12:].split()[:3])
            rows.append((*hkl, f * fom, phi))
    return np.array(rows, dtype=np.float32)


def pha_to_map(
    pha_file: Path,
    structure: gemmi.Structure,
    exact_size: list[int] | None = None,
    sample_rate: float = 3.0,
) -> gemmi.FloatGrid:
    """Calculate a map from a .pha file using the unit cell and space group of the
    structure. If exact_size is given the map is calculated on a grid of that size,
    otherwise the grid spacing is set from the resolution and sample_rate.
    """
    mtz = gemmi.Mtz(with_base=True)
    mtz.spacegroup = structure.find_spacegroup()
    mtz.set_cell_for_all(structure.cell)
    mtz.add_dataset("anode")
    mtz.add_column("F", "F")
    mtz.add_column("PHI", "P")
    mtz.set_data(read_pha(pha_file))
    if exact_size:
        return mtz.transform_f_phi_to_map("F", "PHI", exact_size=exact_size)
    return mtz.transform_f_phi_to_map("F", "PHI", sample_rate=sample_rate)


def move_near_model(
    pos: gemmi.Position, structure: gemmi.Structure, ns: gemmi.NeighborSearch
) -> gemmi.Position:
    """Return the symmetry equivalent of pos that is closest to an atom in the model"""
    mark = ns.find_nearest_atom(pos)
    if mark is None:
        return pos
    atom = mark.to_cra(structure[0]).atom
    image = structure.cell.find_nearest_image(atom.pos, pos, gemmi.Asu.Any)
    return structure.cell.find_nearest_pbc_position(atom.pos, pos, image.sym_idx)


def make_double_diff_map_and_get_peaks(
    pha_above: Path,
    pha_below: Path,
    pdb_file: Path,
    map_out: Path,
    rmsd_threshold: float,
    max_peaks: int,
) -> list[PeakData]:
    """Generate a double difference map by subtracting the anomalous map below the metal
    absorption edge from the map above the edge. Any peaks above the rmsd_threshold are
    then found and the coordinates and peak heights in units of rmsd and map density are
    returned, largest first, with coordinates moved to the symmetry copy closest to the model.
    """
    structure = gemmi.read_structure(str(pdb_file))

    map_diff = pha_to_map(pha_above, structure)
    map_below = pha_to_map(
        pha_below, structure, exact_size=[map_diff.nu, map_diff.nv, map_diff.nw]
    )
    map_diff.array[:] -= map_below.array

    ccp4 = gemmi.Ccp4Map()
    ccp4.grid = map_diff
    ccp4.update_ccp4_header()
    ccp4.write_ccp4_map(str(map_out))
    logging.info(f"Double difference map written to {map_out}")

    logging.info("\n## Finding peaks in double difference map ##")
    rms = float(np.std(map_diff.array))
    logging.info(f"Map rmsd: {rms:.4f}")
    blobs = gemmi.find_blobs_by_flood_fill(
        map_diff,
        cutoff=rmsd_threshold * rms,
        min_volume=0.0,
        min_score=0.0,
        min_peak=rmsd_threshold * rms,
    )
    blobs.sort(key=lambda blob: blob.peak_value, reverse=True)

    ns = gemmi.NeighborSearch(structure[0], structure.cell, 5).populate()
    peaks: list[PeakData] = []
    for blob in blobs:
        pos = move_near_model(blob.centroid, structure, ns)
        # Skip symmetry copies of peaks that have already been found
        if any(
            structure.cell.find_nearest_image(
                gemmi.Position(*peak.xyz), pos, gemmi.Asu.Any
            ).dist()
            < 1.0
            for peak in peaks
        ):
            continue
        if len(peaks) == max_peaks:
            logging.warning(
                f"Found more peaks than the set maximum of {max_peaks} - storing only the largest {max_peaks}"
            )
            break
        peaks.append(
            PeakData(
                density=round(blob.peak_value, 2),
                rmsd=round(blob.peak_value / rms, 2),
                xyz=(round(pos.x, 2), round(pos.y, 2), round(pos.z, 2)),
            )
        )

    return peaks


def calc_double_diff_maps(
    pdb_file: Path,
    pha_above: Path,
    pha_below: Path,
    output_dir: Path,
    peak_threshold: float,
    max_peaks: int,
) -> dict[str, Path]:
    # Check file inputs
    for file_type, file_path in [
        ("AnoDe map above", pha_above),
        ("AnoDe map below", pha_below),
    ]:
        if not file_path.is_file():
            logging.error(f"Could not find {file_type}, expected at: {file_path}")
            return {}

    logging.info("Copying input files to working directory")
    for file, filename in [
        (pdb_file, "final.pdb"),
        (pha_above, "above.pha"),
        (pha_below, "below.pha"),
    ]:
        shutil.copyfile(file, output_dir / filename)

    logging.info("Making double difference map")
    logging.info(f"Using {pdb_file} as reference coordinates for map")
    map_out = output_dir / "diff.map"

    peak_data = make_double_diff_map_and_get_peaks(
        pha_above,
        pha_below,
        pdb_file,
        map_out,
        peak_threshold,
        max_peaks,
    )

    files_out = {"pdb": output_dir / "final.pdb", "map": map_out}

    if not peak_data:
        logging.info(
            f"\nNo peaks found above the threshold of {peak_threshold} rmsd.\n"
        )

    else:
        # Print the extracted information
        logging.info(
            f"\nThe largest peaks (up to a maximum of {max_peaks} peaks) found above the threshold of {peak_threshold} rmsd:"
        )

        peak_file = output_dir / "found_peaks.dat"

        with open(peak_file, "w") as fh:
            for peak_num, peak in enumerate(peak_data, start=1):
                line = f"Peak {peak_num}: Electron Density = {peak.density} e/Å^3, RMSD = {peak.rmsd}, XYZ = {peak.xyz}"
                logging.info(line)
                fh.write(line + "\n")

        files_out["peaks"] = peak_file

    return files_out
