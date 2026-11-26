import logging
import shutil
import sys
import urllib.request
from pathlib import Path

import gemmi

from metal_id.metal_id_helpers import run_ccp4

ESSENTIAL_LABELS = ["IMEAN", "SIGIMEAN", "I(+)", "SIGI(+)", "I(-)", "SIGI(-)"]


def check_mtz(mtz_file: Path) -> float:
    """
    Check that an mtz file has the intensity columns metal_id needs and return the
    wavelength of the data. Exits if a column or the wavelength is missing.
    """
    mtz = gemmi.read_mtz_file(str(mtz_file))
    column_labels = [column.label for column in mtz.columns]
    for label in ESSENTIAL_LABELS:
        if label not in column_labels:
            logging.error(
                f"Input MTZ file {mtz_file} missing essential column label '{label}' - cannot run metal_id"
            )
            sys.exit(1)

    # Use the dataset that holds the anomalous data, as the base dataset (HKL_base)
    # may have a placeholder wavelength
    wavelength = mtz.column_with_label("I(+)").dataset.wavelength
    if wavelength <= 0:
        logging.error(f"No wavelength found in {mtz_file} - cannot run metal_id")
        sys.exit(1)
    return wavelength


def mtz_symmetry(mtz_file: Path) -> tuple[gemmi.UnitCell, gemmi.SpaceGroup | None]:
    mtz = gemmi.read_mtz_file(str(mtz_file), with_data=False)
    return mtz.cell, mtz.spacegroup


def model_symmetry(
    model_file: Path,
) -> tuple[gemmi.UnitCell, gemmi.SpaceGroup | None]:
    structure = gemmi.read_structure(str(model_file))
    return structure.cell, structure.find_spacegroup()


def symmetry_difference(
    sym_1: tuple[gemmi.UnitCell, gemmi.SpaceGroup | None],
    sym_2: tuple[gemmi.UnitCell, gemmi.SpaceGroup | None],
    rel_length_tolerance: float,
    angle_tolerance: float,
) -> str | None:
    """
    Compare two (unit cell, space group) pairs and return a description of why they
    are not compatible, or None if they are.

    Space groups are compatible if they have the same point group and lattice centring,
    as pointless can then reindex the data into the same space group. Unit cells are
    compared after Niggli reduction, so that cells that differ only by the choice of
    axes (alternative indexing) are treated as the same.
    """
    (cell_1, sg_1), (cell_2, sg_2) = sym_1, sym_2
    if (
        sg_1 is None
        or sg_2 is None
        or not cell_1.is_crystal()
        or not cell_2.is_crystal()
    ):
        return "missing space group or unit cell"
    if (sg_1.point_group_hm(), sg_1.centring_type()) != (
        sg_2.point_group_hm(),
        sg_2.centring_type(),
    ):
        return f"different point group or lattice ({sg_1.hm} and {sg_2.hm})"

    reduced = []
    for cell, sg in [(cell_1, sg_1), (cell_2, sg_2)]:
        gv = gemmi.GruberVector(cell, sg.centring_type())
        gv.niggli_reduce()
        reduced.append(gv.cell_parameters())
    lengths_1, lengths_2 = reduced[0][:3], reduced[1][:3]
    angles_1, angles_2 = reduced[0][3:], reduced[1][3:]
    if any(
        abs(l1 - l2) > rel_length_tolerance * l1 for l1, l2 in zip(lengths_1, lengths_2)
    ) or any(abs(a1 - a2) > angle_tolerance for a1, a2 in zip(angles_1, angles_2)):
        return f"unit cells differ ({format_cell(cell_1)} and {format_cell(cell_2)})"
    return None


def format_cell(cell: gemmi.UnitCell) -> str:
    return " ".join(f"{p:.2f}" for p in cell.parameters)


def fetch_pdb(code: str, output_dir: Path) -> Path | None:
    """
    Download a model from the PDB and save it in PDB format. Returns None if the
    download fails.
    """
    pdb_file = output_dir / f"{code.upper()}.pdb"
    for url, fmt in [
        (f"https://files.rcsb.org/download/{code.upper()}.pdb", "pdb"),
        (f"https://files.rcsb.org/download/{code.upper()}.cif", "cif"),
    ]:
        try:
            with urllib.request.urlopen(url, timeout=30) as response:
                text = response.read().decode()
        except OSError as error:
            logging.info(f"Could not download {url}: {error}")
            continue
        if fmt == "pdb":
            pdb_file.write_text(text)
        else:
            structure = gemmi.read_structure_string(text)
            structure.write_pdb(str(pdb_file))
        logging.info(f"Downloaded {code} from {url} to {pdb_file}")
        return pdb_file
    return None


def find_matching_model(
    models: list[Path],
    data_symmetries: dict[str, tuple[gemmi.UnitCell, gemmi.SpaceGroup | None]],
    rel_length_tolerance: float,
    angle_tolerance: float,
) -> Path | None:
    """Return the first model whose symmetry matches all of the datasets"""
    for model in models:
        model_sym = model_symmetry(model)
        problems = {
            name: symmetry_difference(
                model_sym, data_sym, rel_length_tolerance, angle_tolerance
            )
            for name, data_sym in data_symmetries.items()
        }
        if not any(problems.values()):
            logging.info(f"Model {model} matches the symmetry of both datasets")
            return model
        for name, problem in problems.items():
            if problem:
                logging.info(f"Model {model} does not match the {name} data: {problem}")
    return None


def reindex_to_model(mtz_file: Path, model: Path, work_dir: Path) -> Path:
    """
    Run pointless to put the data in the same space group and indexing as the model.
    Returns the path of the reindexed mtz file.
    """
    hklout = work_dir / "reindexed.mtz"
    logging.info(f"Running pointless to reindex {mtz_file} to match {model}")
    pointless_log = run_ccp4(
        "pointless",
        ["HKLIN", mtz_file, "HKLOUT", hklout, "XYZIN", model],
        work_dir / "pointless.log",
        keywords=["TOLERANCE 5", "END"],
    )
    if "Incompatible symmetries" in pointless_log or not hklout.is_file():
        logging.error(
            f"pointless could not match {mtz_file} to {model}, see {work_dir / 'pointless.log'}"
        )
        sys.exit(1)
    return hklout


def copy_model(model: Path, output_dir: Path) -> Path:
    """Copy a model to the output directory in PDB format, as ANODE needs a PDB file"""
    pdb_file = output_dir / f"{model.stem}.pdb"
    if pdb_file.exists() and pdb_file.samefile(model):
        return pdb_file
    if model.suffix.lower() == ".pdb":
        shutil.copyfile(model, pdb_file)
    else:
        gemmi.read_structure(str(model)).write_pdb(str(pdb_file))
    return pdb_file
