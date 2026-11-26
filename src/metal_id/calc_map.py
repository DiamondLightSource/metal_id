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
    with columns h, k, l, F*FOM, phi (degrees).

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


def coefficients_to_map(
    coefficients: np.ndarray,
    structure: gemmi.Structure,
    exact_size: list[int] | None = None,
    sample_rate: float = 3.0,
) -> gemmi.FloatGrid:
    """Calculate a map from an array with columns h, k, l, F, phi (degrees) using the
    unit cell and space group of the structure. If exact_size is given the map is
    calculated on a grid of that size, otherwise the grid spacing is set from the
    resolution and sample_rate.
    """
    mtz = gemmi.Mtz(with_base=True)
    mtz.spacegroup = structure.find_spacegroup()
    mtz.set_cell_for_all(structure.cell)
    mtz.add_dataset("map")
    mtz.add_column("F", "F")
    mtz.add_column("PHI", "P")
    mtz.set_data(coefficients.astype(np.float32))
    if exact_size:
        return mtz.transform_f_phi_to_map("F", "PHI", exact_size=exact_size)
    return mtz.transform_f_phi_to_map("F", "PHI", sample_rate=sample_rate)


def grid_size(grid: gemmi.FloatGrid) -> list[int]:
    return [grid.nu, grid.nv, grid.nw]


def anode_maps(
    pha_above: Path, pha_below: Path, structure: gemmi.Structure
) -> tuple[gemmi.FloatGrid, gemmi.FloatGrid]:
    """Calculate the above and below anomalous maps, on the same grid, from ANODE .pha
    files. ANODE normalises the data, so the maps are in units of σ.
    """
    map_above = coefficients_to_map(read_pha(pha_above), structure)
    map_below = coefficients_to_map(
        read_pha(pha_below), structure, exact_size=grid_size(map_above)
    )
    return map_above, map_below


def read_map_coefficients(
    mtz_file: Path, f_label: str = "FAN", phi_label: str = "PHAN"
) -> tuple[np.ndarray, np.ndarray]:
    """Read map coefficients from an mtz file, returning an array with columns h, k, l,
    F, phi and an array of the resolution of each reflection. Reflections without
    coefficients are left out.
    """
    mtz = gemmi.read_mtz_file(str(mtz_file))
    f = mtz.column_with_label(f_label).array
    phi = mtz.column_with_label(phi_label).array
    ok = np.isfinite(f) & np.isfinite(phi)
    coefficients = np.column_stack([mtz.make_miller_array(), f, phi])[ok]
    return coefficients, mtz.make_d_array()[ok]


def hkl_keys(coefficients: np.ndarray) -> np.ndarray:
    """A unique integer for each h, k, l (indices must be within ±1000)"""
    hkl = coefficients[:, :3].astype(np.int64) + 1024
    return (hkl[:, 0] * 2048 + hkl[:, 1]) * 2048 + hkl[:, 2]


def refmac_anomalous_maps(
    mtz_above: Path, mtz_below: Path, structure: gemmi.Structure
) -> tuple[gemmi.FloatGrid, gemmi.FloatGrid]:
    """Calculate the above and below anomalous maps, on the same grid, from the refmac
    FAN/PHAN coefficients. Only reflections present in both datasets, to the resolution
    of the lower resolution dataset, are used, so that the maps can be subtracted. The
    data were scaled against the model by refmac, so the maps are on an approximately
    absolute scale (e/Å^3).
    """
    coeffs_above, d_above = read_map_coefficients(mtz_above)
    coeffs_below, d_below = read_map_coefficients(mtz_below)
    d_min = max(d_above.min(), d_below.min())
    _, in_above, in_below = np.intersect1d(
        hkl_keys(coeffs_above), hkl_keys(coeffs_below), return_indices=True
    )
    keep = d_above[in_above] >= d_min
    in_above, in_below = in_above[keep], in_below[keep]
    logging.info(
        f"Anomalous map coefficients: above {len(coeffs_above)} reflections (to {d_above.min():.2f} Å), "
        f"below {len(coeffs_below)} (to {d_below.min():.2f} Å). "
        f"Using the {len(in_above)} reflections common to both, to {d_min:.2f} Å"
    )

    map_above = coefficients_to_map(coeffs_above[in_above], structure)
    map_below = coefficients_to_map(
        coeffs_below[in_below], structure, exact_size=grid_size(map_above)
    )
    return map_above, map_below


def write_map(grid: gemmi.FloatGrid, map_file: Path):
    ccp4 = gemmi.Ccp4Map()
    ccp4.grid = grid
    ccp4.update_ccp4_header()
    ccp4.write_ccp4_map(str(map_file))


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


def find_peaks(
    grid: gemmi.FloatGrid,
    structure: gemmi.Structure,
    rmsd_threshold: float,
    max_peaks: int,
) -> list[PeakData]:
    """Find peaks above rmsd_threshold in the map and return their coordinates and peak
    heights in units of rmsd and map density, largest first, with coordinates moved to
    the symmetry copy closest to the model.
    """
    rms = float(np.std(grid.array))
    logging.info(f"Map rmsd: {rms:.4f}")
    blobs = gemmi.find_blobs_by_flood_fill(
        grid,
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
    map_above: gemmi.FloatGrid,
    map_below: gemmi.FloatGrid,
    output_dir: Path,
    peak_threshold: float,
    max_peaks: int,
) -> dict[str, Path]:
    """Subtract the anomalous map below the absorption edge from the map above it, write
    out the maps and find the peaks in the double difference map.
    """
    structure = gemmi.read_structure(str(pdb_file))
    if pdb_file != output_dir / "final.pdb":
        shutil.copyfile(pdb_file, output_dir / "final.pdb")

    write_map(map_above, output_dir / "anom_above.map")
    write_map(map_below, output_dir / "anom_below.map")

    logging.info("Making double difference map")
    logging.info(f"Using {pdb_file} as reference coordinates for map")
    map_diff = map_above.clone()
    map_diff.array[:] -= map_below.array
    map_out = output_dir / "diff.map"
    write_map(map_diff, map_out)
    logging.info(f"Double difference map written to {map_out}")

    logging.info("\n## Finding peaks in double difference map ##")
    peak_data = find_peaks(map_diff, structure, peak_threshold, max_peaks)

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
