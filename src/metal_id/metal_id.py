import argparse
import shutil
import sys
import logging.config
from pathlib import Path

import gemmi
from metal_id_helpers import (
    ensure_unique_directory,
    PDBFileOrCode,
    mtz_exists,
    run_dimple,
    generate_coot_viewer_script,
)
from prepare_data import (
    check_mtz,
    copy_model,
    fetch_pdb,
    find_matching_model,
    mtz_symmetry,
    reindex_to_model,
    symmetry_difference,
)
from anomalous import run_anode, run_ctruncate, run_refmac_anomalous
from calc_map import anode_maps, calc_double_diff_maps, refmac_anomalous_maps
from viewer import gen_html_metal_id
from version import get_version

# Command line interface
parser = argparse.ArgumentParser(
    prog="metal_id",
    description="Locate a given element from data collected above and below the element's absorption edge.",
    usage="%(prog)s mtz_above mtz_below pdb [additional pdbs] [options]",
    epilog="metal_id is still in development, if you have any issues with it, please contact Phil Blowey at philip.blowey@diamond.ac.uk",
)

parser.add_argument(
    "--version",
    action="version",
    version=f"%(prog)s {get_version()}",
)
parser.add_argument(
    "mtz_above",
    type=mtz_exists,
    help="Path to mtz file containing data collected above the absorption edge",
)
parser.add_argument(
    "mtz_below",
    type=mtz_exists,
    help="Path to mtz file containing data collected below the absorption edge",
)
parser.add_argument(
    "pdb",
    type=PDBFileOrCode,
    nargs="+",
    help="Path to pdb file(s) and/or 4 character PDB codes",
)
parser.add_argument(
    "-o",
    "--output",
    type=Path,
    help="Path to output directory. If the path already exists, a numerical suffix will be added. Default is metal_id",
    default=Path("metal_id"),
)
parser.add_argument(
    "--peak-threshold",
    type=float,
    help="Set the peak height threshold that the peaks must exceed in order to be detected",
    default=5.0,
)
parser.add_argument(
    "--max-peaks",
    type=int,
    help="Set the maximum number of peaks to detect",
    default=10,
)
parser.add_argument(
    "--method",
    choices=["refmac", "anode"],
    help="Program used to calculate the anomalous maps. refmac scales the data against the model, so the maps are on an approximately absolute scale; ANODE maps are normalised. Default is refmac",
    default="refmac",
)
parser.add_argument(
    "--cell-tolerance",
    type=float,
    help="Maximum difference in unit cell lengths (in percent) for the datasets and model to be treated as matching. If the model does not match the data, dimple is run to get a model that does. Default is 1.0",
    default=1.0,
)
parser.add_argument(
    "--max-r-factor",
    type=float,
    help="Stop if the R factor from refmac against either dataset is above this value. Default is 0.35",
    default=0.35,
)

args = parser.parse_args()

mtz_above = args.mtz_above.resolve()
mtz_below = args.mtz_below.resolve()
pdb = args.pdb
output_dir = args.output.resolve()
peak_threshold = args.peak_threshold
max_peaks = args.max_peaks
method = args.method
rel_length_tolerance = args.cell_tolerance / 100
angle_tolerance = 1.0

if output_dir.exists():
    output_dir = ensure_unique_directory(output_dir)

output_dir.mkdir(parents=True)

METAL_ID_LOGGING = {
    "version": 1,
    "disable_existing_loggers": False,
    "formatters": {
        "default": {"format": "%(message)s"},
    },
    "handlers": {
        "console": {"class": "logging.StreamHandler", "formatter": "default"},
        "file": {
            "class": "logging.FileHandler",
            "formatter": "default",
            "filename": f"{output_dir / 'metal_id.log'}",
            "mode": "w",
            "encoding": "utf-8",
        },
    },
    "root": {"level": "INFO", "handlers": ["console", "file"]},
}

logging.config.dictConfig(METAL_ID_LOGGING)

logging.info(f"metal_id version: {get_version()}")
logging.info(f"Running command: metal_id {' '.join(sys.argv[1::])}")

logging.info(f"Outputting to: {output_dir}")

for arg, arg_name in [
    (mtz_above, "mtz_above"),
    (mtz_below, "mtz_below"),
    ([str(pdb_file_or_code.value) for pdb_file_or_code in pdb], "pdb"),
    (output_dir, "output_dir"),
    (method, "method"),
]:
    logging.info(f"{arg_name} = {arg}")

# Check the data

logging.info("\n### Checking above and below data ###\n")
mtz_files = {"above": mtz_above, "below": mtz_below}
wavelengths = {name: check_mtz(mtz_file) for name, mtz_file in mtz_files.items()}
data_symmetries = {name: mtz_symmetry(mtz_file) for name, mtz_file in mtz_files.items()}
for name in mtz_files:
    cell, spacegroup = data_symmetries[name]
    logging.info(
        f"{name}: wavelength {wavelengths[name]:.4f} Å, space group "
        f"{spacegroup.hm if spacegroup else 'unknown'}, cell {' '.join(f'{p:.2f}' for p in cell.parameters)}"
    )

problem = symmetry_difference(
    data_symmetries["above"],
    data_symmetries["below"],
    rel_length_tolerance,
    angle_tolerance,
)
if problem:
    logging.error(
        f"ERROR: MTZ files have incompatible symmetry ({problem}) - cannot run metal_id"
    )
    sys.exit(1)
logging.info("Above and below data have compatible symmetry")

# Find a model that matches the data, running dimple to get one if needed

logging.info("\n### Checking model(s) against the data ###\n")
models = []
# Models for dimple: PDB codes that couldn't be downloaded are left for dimple to fetch
dimple_models = []
for pdb_file_or_code in pdb:
    if pdb_file_or_code.is_file:
        model = pdb_file_or_code.value
    else:
        model = fetch_pdb(pdb_file_or_code.value, output_dir)
    if model is None:
        logging.warning(f"Could not download {pdb_file_or_code.value}")
        dimple_models.append(str(pdb_file_or_code.value))
    else:
        models.append(model)
        dimple_models.append(str(model))

model = find_matching_model(
    models, data_symmetries, rel_length_tolerance, angle_tolerance
)

if model is None:
    logging.info(
        "\nNo model matches the symmetry of the data, running dimple on the 'above' data to get one"
    )
    dimple_dir = output_dir / "dimple"
    dimple_output = run_dimple(mtz_above, dimple_models, dimple_dir)
    logging.info(f"Captured output from dimple: \n {dimple_output.stdout}")
    if dimple_output.returncode != 0 or not (dimple_dir / "final.pdb").is_file():
        logging.error(f"dimple failed, see {dimple_dir / 'dimple.log'}")
        sys.exit(1)
    model = find_matching_model(
        [dimple_dir / "final.pdb"],
        data_symmetries,
        rel_length_tolerance,
        angle_tolerance,
    )
    if model is None:
        logging.error(
            "The model from dimple does not match both datasets - cannot run metal_id"
        )
        sys.exit(1)

model = copy_model(model, output_dir)
logging.info(f"Using model {model} for both datasets")

# Reindex the data against the model and scale it against the model with refmac

refmac_mtzs = {}
pha_files = {}
for name, mtz_file in mtz_files.items():
    logging.info(f"\n### Processing the '{name}' data ###\n")
    work_dir = output_dir / name
    work_dir.mkdir()
    reindexed = reindex_to_model(mtz_file, model, work_dir)
    truncated = run_ctruncate(reindexed, work_dir)
    refmac_mtzs[name], r_factor = run_refmac_anomalous(
        truncated, model, wavelengths[name], work_dir
    )
    logging.info(f"R factor of the model against the '{name}' data: {r_factor:.4f}")
    if r_factor > args.max_r_factor:
        logging.error(
            f"R factor against the '{name}' data is above {args.max_r_factor} - the model does not fit the data well enough. "
            "Refine the model against the data first (e.g. with dimple) and rerun metal_id with the refined model"
        )
        sys.exit(1)
    shutil.copy(refmac_mtzs[name], output_dir / f"refmac_{name}.mtz")
    if method == "anode":
        pha_files[name] = run_anode(reindexed, model, work_dir)

logging.info("\n### Calculating map of element location ###\n")
structure = gemmi.read_structure(str(model))
if method == "refmac":
    map_above, map_below = refmac_anomalous_maps(
        refmac_mtzs["above"], refmac_mtzs["below"], structure
    )
else:
    map_above, map_below = anode_maps(pha_files["above"], pha_files["below"], structure)

files_out = calc_double_diff_maps(
    model, map_above, map_below, output_dir, peak_threshold, max_peaks
)
if files_out:
    logging.info("\n### Generating coot viewer script ###\n")
    coot_viewer_script = generate_coot_viewer_script(
        files_out["pdb"], files_out["map"], peak_threshold
    )
    logging.info(f"Coot script for viewing results '{coot_viewer_script}' generated\n")
    logging.info("To use run command:\n")
    logging.info(f"coot -s {coot_viewer_script}\n")

if "peaks" in files_out:
    logging.info("\n### Generating viewer file ###\n")
    gen_html_metal_id(
        results_directory=output_dir,
        pdb_file=files_out["pdb"],
        mtz_file=output_dir / "refmac_below.mtz",
        diff_map_file=files_out["map"],
        peak_file=files_out["peaks"],
        isovalue=peak_threshold,
    )
    logging.info(
        f"Molstar visualization generated to {output_dir / 'metal_id.html'} \n"
    )

logging.info("\n### End of script ###")
