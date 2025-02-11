from pathlib import Path
import argparse
import subprocess
import logging


def ensure_unique_directory(path):
    path = Path(path)
    counter = 1
    new_path = path

    while new_path.exists():
        new_path = path.with_name(f"{path.stem}_{counter}{path.suffix}")
        counter += 1

    return new_path


def file_or_code(value):
    if Path(value).is_file():
        return Path(value)
    elif len(value) == 4 and value.isalnum():
        return value
    else:
        raise argparse.ArgumentTypeError(
            f"Invalid pdb input: '{value}', must be valid four-character pdb code or an existing file path."
        )


def mtz_exists(input_file):
    if not Path(input_file).is_file():
        raise argparse.ArgumentTypeError(f"File not found: '{input_file}'.")
    return Path(input_file)


# Class for handling PDB as a file path or a four character code.
class PDBFileOrCode:
    def __init__(self, file_or_code):
        if Path(file_or_code).is_file():
            self.value = Path(file_or_code)
            self.is_file = True
        elif len(file_or_code) == 4 and file_or_code.isalnum():
            self.value = str(file_or_code)
            self.is_file = False
        else:
            raise ValueError(
                f"Invalid input '{file_or_code}'. pdb must be given as a valid file or 4-character pdb code"
            )


def run_dimple(mtz, pdb, dimple_dir):
    if isinstance(pdb, list):
        pdb = " ".join(pdb)
    dimple_command = f"dimple {mtz} {pdb} {dimple_dir} --anode -fpng"

    logging.info(f"Running dimple with command:\n\n{dimple_command}\n")
    dimple_output = subprocess.run(
        dimple_command,
        shell=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
    )

    return dimple_output
