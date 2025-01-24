import logging
import shutil
import subprocess

from iotbx import mtz


def find_cols_from_type(obj, type, file="mtz_file"):
    """Get the label and corresponding sigma label for a specifed data type from an mtz file.

    Function takes an iotbx mtz.object and mtz column data-type identifier as input and returns
    the corresponding column heading and sigma column heading. Function will only return the
    first column heading of that type.
    """
    col_types = obj.column_types()
    col_labs = obj.column_labels()
    if type in col_types:
        indices = [_index for _index, _value in enumerate(col_types) if _value == type]
        if len(indices) > 1:
            logging.warning(
                f"Multiple {mtz.column_type_legend[type]} data columns found in {file}, using the first one"
            )
        col_lab = col_labs[indices[0]]
    else:
        logging.error(
            f"Could not find {mtz.column_type_legend[type]} data column in {file}"
        )
        return None

    if (sig_col_lab := "SIG" + col_lab) not in col_labs:
        logging.error(f"Could not find {sig_col_lab} data in {file}")
        return None
    return col_lab, sig_col_lab


def calc_amplitudes(mtz_obj, mtz_file, output_dir):
    """
    Use truncate to calculate amplitudes from IMEAN data.

    Takes mtz object and mtz file name and runs truncate to calculate
    amplitudes. Returns a new mtz_object and file name for the output file
    with suffix: "_amplit.mtz"
    """
    if "F" not in mtz_obj.column_types():
        _col_lab, _sig_col_lab = find_cols_from_type(mtz_obj, "J", mtz_file)
        logging.info(f"Amplitude data not in {mtz_file}, running TRUNCATE to calculate")
        amplit_file = output_dir / f"{mtz_file.stem}_amplit.mtz"
        truncate_script = [
            f"truncate hklin {mtz_file} hklout {amplit_file} <<END-TRUNCATE",
            f"labin IMEAN={_col_lab} SIGIMEAN={_sig_col_lab}",
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

    for input_file in [mtz_above, mtz_below]:
        shutil.copy(input_file, scaling_dir)

    mtz_above = scaling_dir / mtz_above.name
    mtz_below = scaling_dir / mtz_below.name

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
        logging.error("mtz files have incompatible symmetry")
        return False
    # Update mtz_der to the reindexed file path
    mtz_above = hklout

    # Read in mtz files
    obj_below = mtz.object(str(mtz_below))
    obj_above = mtz.object(str(mtz_above))

    # Calculate structure factors if needed using truncate
    obj_below, mtz_below = calc_amplitudes(obj_below, mtz_below, scaling_dir)
    obj_above, mtz_above = calc_amplitudes(obj_above, mtz_above, scaling_dir)

    col_labs = {}
    col_params = [
        ("F_nat", "SIGF_nat", obj_below, "F", mtz_below),
        ("F_der", "SIGF_der", obj_above, "F", mtz_above),
        ("DANO_der", "SIGDANO_der", obj_above, "D", mtz_above),
    ]

    for _val, _sigval, obj, type, file in col_params:
        try:
            col_labs[_val], col_labs[_sigval] = find_cols_from_type(obj, type, file)
        except TypeError:
            logging.error(f"{_val} and/or {_sigval} missing from {file}")
            return False

    # Add the F and SIGF data from one file to the other with cad
    mtz_combi = scaling_dir / f"{mtz_above.stem}_combined.mtz"
    # Get list of column headers excluding hkl
    col_labs_der = obj_above.crystals()[1].datasets()[0].column_labels()
    # Convert list to cad input format
    labin_der = [f"E{_i+1}={_label}" for _i, _label in enumerate(col_labs_der)]
    cad_script = [
        f"cad hklin1 {mtz_above} hklin2 {mtz_below} hklout {mtz_combi} <<END-CAD",
        "TITLE Add data for scaling",
        f"LABIN FILE 2 E1={col_labs['F_nat']} E2={col_labs['SIGF_nat']}",
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
        f"LABIN FP=Fscale SIGFP=SIGFscale FPH1={col_labs['F_der']} SIGFPH1={col_labs['SIGF_der']} DPH1={col_labs['DANO_der']} SIGDPH1={col_labs['SIGDANO_der']}",
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
