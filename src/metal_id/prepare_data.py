import logging
import subprocess
import sys

from iotbx import mtz


def calc_amplitudes(mtz_obj, mtz_file, output_dir):
    """
    Use truncate to calculate amplitudes from IMEAN data.

    Takes mtz object and mtz file name and runs truncate to calculate
    amplitudes if they are missing. Returns the file name of the output file
    with suffix: "_amplit.mtz", or the input file name if no amplitudes were needed.
    """
    if "F" not in mtz_obj.column_labels():
        logging.info(f"Amplitude data not in {mtz_file}, running TRUNCATE to calculate")
        amplit_file = output_dir / f"{mtz_file.stem}_amplit.mtz"
        truncate_script = [
            f"truncate hklin {mtz_file} hklout {amplit_file} <<END-TRUNCATE",
            "labin IMEAN=IMEAN SIGIMEAN=SIGIMEAN",
            "labout F=F SIGF=SIGF",
            "NOHARVEST",
            "END",
            "END-TRUNCATE",
        ]
        ccp4_command(truncate_script, "truncate.log", output_dir)
        return amplit_file
    else:
        return mtz_file


def ccp4_command(script, output, output_dir):
    command = "\n".join(script)
    logging.info(f"Running command:\n{command}")
    result = subprocess.run(
        command,
        shell=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
    )
    with open(output_dir / output, "w") as log_file:
        log_file.write(result.stdout)
    return result.stdout


def prepare_data(mtz_above, mtz_below, output_dir):
    """
    Check that the above and below data are compatible and ready for dimple.

    Runs pointless to put the above data into the same space group and indexing as the
    below data, checks that both files have the required intensity columns, and
    calculates amplitudes with truncate if they are missing. Returns the paths of the
    prepared above and below mtz files.

    No scaling is done here: dimple runs ANODE on the intensities, and SHELXC/ANODE
    normalise the anomalous data, so any scale applied beforehand has no effect on the
    anomalous maps.
    """
    prepare_dir = output_dir / "prepare_data"
    prepare_dir.mkdir()

    # Ensure that the mtz files have compatible symmetry and put them into the same space group using pointless
    hklout = prepare_dir / f"{mtz_above.stem}_reindexed.mtz"
    pointless_command = [
        f"pointless hklin {mtz_above} hklout {hklout} hklref {mtz_below}"
    ]

    logging.info(
        "Running pointless to ensure above and below data are in the same space group"
    )
    pointless_log = ccp4_command(pointless_command, "pointless.log", prepare_dir)

    if "Incompatible symmetries" in pointless_log:
        logging.error(
            "ERROR: MTZ files have incompatible symmetry - cannot run metal_id"
        )
        sys.exit(1)
    # Update mtz_above to the reindexed file path
    mtz_above = hklout

    # Read in mtz files as iotbx mtz objects
    obj_below = mtz.object(str(mtz_below))
    obj_above = mtz.object(str(mtz_above))

    # Check that files meet minimum requirements
    essential_labels = ["IMEAN", "SIGIMEAN", "I(+)", "SIGI(+)", "I(-)", "SIGI(-)"]
    for mtz_object in [obj_above, obj_below]:
        column_labels = mtz_object.crystals()[1].datasets()[0].column_labels()
        for label in essential_labels:
            if label not in column_labels:
                logging.error(
                    f"Input MTZ file missing essential column label '{label}' - cannot run metal_id"
                )
                sys.exit(1)

    # Calculate structure factors if needed using truncate
    mtz_below = calc_amplitudes(obj_below, mtz_below, prepare_dir)
    mtz_above = calc_amplitudes(obj_above, mtz_above, prepare_dir)

    return mtz_above, mtz_below
