import argparse
import shutil
import sys
import logging.config
from pathlib import Path
from metal_id_helpers import (
    ensure_unique_directory,
    PDBFileOrCode,
    mtz_exists,
    run_dimple,
)
from scaling import scale_data
from calc_map import calc_double_diff_maps

# Command line interface
parser = argparse.ArgumentParser(
    prog="metal_id",
    description="Locate a given element from data collected above and below the element's absorption edge.",
    usage="%(prog)s mtz_above mtz_below pdb [additional pdbs] [options]",
    epilog="metal_id is still in development, if you have any issues with it, please contact Phil Blowey at philip.blowey@diamond.ac.uk",
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
    nargs="*",
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

args = parser.parse_args()

mtz_above = args.mtz_above.resolve()
mtz_below = args.mtz_below.resolve()
pdb = args.pdb
output_dir = args.output.resolve()
peak_threshold = args.peak_threshold
max_peaks = args.max_peaks
fcolumn_label = args.fcolumn


for arg, arg_name in [
    (mtz_above, "mtz_above"),
    (mtz_below, "mtz_below"),
    (pdb, "pdb"),
    (output_dir, "output_dir"),
]:
    logging.info(f"{arg_name} = {arg}")

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

logging.info(f"Running command: {' '.join(sys.argv)}")

logging.info(f"Outputting to: {output_dir}")

# Scale the data

logging.info("\n### Scaling above data relative to the below data ###\n")
mtz_above, mtz_below = scale_data(mtz_above, mtz_below, output_dir, fcolumn_label)
logging.info(f"Scaled above data written to file {mtz_above}")

pdb_files_and_codes = []
for pdb_file_or_code in pdb:
    if pdb_file_or_code.is_file:
        shutil.copy(pdb_file_or_code.value, output_dir)
        pdb_file_or_code.value = output_dir / pdb_file_or_code.value.name
    pdb_files_and_codes.append(str(pdb_file_or_code.value))

dimple_dir_above = output_dir / "dimple_above"

logging.info("\n### Running dimple on the 'above' data ###\n")
dimple_output = run_dimple(mtz_above, pdb_files_and_codes, dimple_dir_above)
logging.info(f"Captured output from dimple: \n {dimple_output.stdout}")

dimple_dir_below = output_dir / "dimple_below"
pdb_above = dimple_dir_above / "final.pdb"
pha_above = dimple_dir_above / "anode.pha"

logging.info("\n### Running dimple on the 'below' data ###\n")
dimple_output = run_dimple(mtz_below, pdb_above, dimple_dir_below)
logging.info(f"\nCaptured output from dimple: \n {dimple_output.stdout}")

pdb_below = dimple_dir_below / "final.pdb"
pha_below = dimple_dir_below / "anode.pha"

logging.info("### Calculating map of element location ###\n")
calc_double_diff_maps(
    pdb_above, pdb_below, pha_above, pha_below, output_dir, peak_threshold, max_peaks
)

logging.info("\n### End of script ###")
