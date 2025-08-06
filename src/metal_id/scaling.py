import logging
import subprocess
import sys

from iotbx import mtz


def calc_amplitudes(mtz_obj, mtz_file, output_dir):
    """
    Use truncate to calculate amplitudes from IMEAN data.

    Takes mtz object and mtz file name and runs truncate to calculate
    amplitudes. Returns a new mtz_object and file name for the output file
    with suffix: "_amplit.mtz"
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

        amplit_obj = mtz.object(str(amplit_file))
        return amplit_obj, amplit_file
    else:
        return mtz_obj, mtz_file


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


def scale_data(mtz_above, mtz_below, output_dir):
    scaling_dir = output_dir / "scaling"
    scaling_dir.mkdir()

    # Ensure that the mtz files have compatible symmetry and put them into the same space group using pointless
    hklout = scaling_dir / f"{mtz_above.stem}_reindexed.mtz"
    pointless_command = [
        f"pointless hklin {mtz_above} hklout {hklout} hklref {mtz_below}"
    ]

    logging.info(
        "Running pointless to ensure above and below data are in the same space group"
    )
    pointless_log = ccp4_command(pointless_command, "pointless.log", scaling_dir)

    if "Incompatible symmetries" in pointless_log:
        logging.error(
            "ERROR: MTZ files have incompatible symmetry - cannot run metal_id"
        )
        sys.exit()
    # Update mtz_der to the reindexed file path
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
                sys.exit()

    # Calculate structure factors if needed using truncate
    obj_below, mtz_below = calc_amplitudes(obj_below, mtz_below, scaling_dir)
    obj_above, mtz_above = calc_amplitudes(obj_above, mtz_above, scaling_dir)

    # Add the F and SIGF data from one file to the other with cad
    mtz_combi = scaling_dir / f"{mtz_above.stem}_combined.mtz"
    # Get list of column headers excluding hkl
    all_labels_der = obj_above.crystals()[1].datasets()[0].column_labels()
    # List of labels to look for in the file
    selected_labels = essential_labels + ["F", "SIGF", "FreeR_flag"]
    # Filtered list of labels that exist in the file
    selected_labels_der = [
        label for label in all_labels_der if label in selected_labels
    ]
    # Convert list to cad input format
    labin_der = [
        f"E{_i}={_label}" for _i, _label in enumerate(selected_labels_der, start=1)
    ]
    cad_script = [
        f"cad hklin1 {mtz_above} hklin2 {mtz_below} hklout {mtz_combi} <<END-CAD",
        "TITLE Add data for scaling",
        "LABIN FILE 2 E1=F E2=SIGF",
        f"LABIN FILE 1 {' '.join(labin_der)}",
        "LABOUT FILE 2 E1=Fscale E2 = SIGFscale",
        "DNAME FILE_NUMBER 2 ALL=refData",
        "SYSAB_KEEP",
        "END",
        "END-CAD",
    ]

    logging.info("Adding reference columns to data for scaling with CAD")
    ccp4_command(cad_script, "cad.log", scaling_dir)

    # Scale the above data using the data added by cad
    mtz_combined_scaled = scaling_dir / f"{mtz_above.stem}_combined_scaled.mtz"
    scaleit_script = [
        f"scaleit hklin {mtz_combi} hklout {mtz_combined_scaled} <<END-SCALEIT",
        "TITLE Scale data using added ref data",
        "LABIN FP=Fscale SIGFP=SIGFscale FPH1=F SIGFPH1=SIGF",
        "AUTO",
        "WEIGHT",
        "REFINE SCALE",
        "END",
        "END-SCALEIT",
    ]

    logging.info("Scaling data with SCALEIT")
    ccp4_command(scaleit_script, "scaleit.log", scaling_dir)

    # Remove the reference columns used for scaling from the mtz file using mtzutils
    logging.info(f"Removing scaling columns from {mtz_combined_scaled} with mtzutils")
    mtz_scaled = scaling_dir / f"{mtz_above.stem}_scaled.mtz"
    mtzutil_script = [
        f"mtzutils hklin {mtz_combined_scaled} hklout {mtz_scaled} <<END-MTZUTILS",
        "EXCLUDE Fscale SIGFscale",
        "END",
        "END-MTZUTILS",
    ]

    logging.info("Running mtzutils to remove reference columns from scaled data")
    ccp4_command(mtzutil_script, "mtzutils.log", scaling_dir)

    return mtz_scaled, mtz_below
