import logging
import math
import os
import re
import shutil
import subprocess

from iotbx import pdb


def view_as_quat(p1, p2):
    """
    Calculate a quaternion representing the rotation necessary to orient a viewer's
    perspective from an initial view direction towards a desired view direction,
    given by the positions p1 and p2, respectively.

    Parameters:
    - p1: tuple or list representing the initial rotation centre (x, y, z) of the viewer.
    - p2: tuple or list representing the desired position (x, y, z) towards which
        the viewer should orient.

    Returns:
    - Quaternion: A tuple representing the quaternion (w, x, y, z) that represents
    the rotation necessary to align the initial view direction with the desired
    view direction. If either p1 or p2 is None, returns the default identity
    quaternion (0., 0., 0., 1.), indicating no rotation.
    """
    if p1 is None or p2 is None:
        return (0.0, 0.0, 0.0, 1.0)
    # Find and normalise direction vector from p1 to p2
    d = (p2[0] - p1[0], p2[1] - p1[1], p2[2] - p1[2])
    length = math.sqrt(d[0] ** 2 + d[1] ** 2 + d[2] ** 2)
    d = (d[0] / length, d[1] / length, d[2] / length)
    # Cross product of d and (0, 0, -1) to view down the direction vector.
    prod = (d[1], -d[0], 0)
    # Generate and normalise quaternion from cross product
    quat = (prod[0], prod[1], prod[2], 1 - d[2])
    qlen = math.sqrt(sum(a * a for a in quat))
    return (quat[0] / qlen, quat[1] / qlen, quat[2] / qlen, quat[3] / qlen)


def are_pdbs_similar(file_1, file_2):
    """
    Determine if two pdb files have the same crystal symmetry, the same number
    and type of atoms and sufficiently similar unit cell and atomic coordinates
    within the defined tolerances
    """

    def read_pdb(file):
        """
        Read a pdb file to get crystal symmetry, atom names and atom coordinates
        """
        pdb_obj = pdb.input(file)
        sym = pdb_obj.crystal_symmetry()
        atoms = pdb_obj.atoms()
        atom_names = atoms.extract_name()
        list_atoms = list(atom_names)
        atom_coords = atoms.extract_xyz()
        list_coords = list(atom_coords)
        return sym, list_atoms, list_coords

    # Read pdb files
    sym_1, atoms_1, coords_1 = read_pdb(str(file_1))
    sym_2, atoms_2, coords_2 = read_pdb(str(file_2))

    # Use default if none set tolerances
    tolerances = {
        "rel_cell_length": 0.01,
        "abs_cell_angle": 1.0,
        "abs_coord_diff": 5.0,  # Units Å
    }

    # Compare symmetry
    is_similar_sym = sym_1.is_similar_symmetry(
        sym_2,
        relative_length_tolerance=tolerances["rel_cell_length"],
        absolute_angle_tolerance=tolerances["abs_cell_angle"],
    )
    if not is_similar_sym:
        logging.error("PDB file symmetries are too different")
        return False

    # Compare atom type/number
    if atoms_1 != atoms_2:
        logging.error("Different number or type of atoms in pdb files")

    # Compare atom coordinates
    combined_coords = zip(coords_1, coords_2)
    for xyz_1, xyz_2 in combined_coords:
        # Calculate the distance between xyz_1 and xyz_2
        diff = abs(
            (
                (xyz_1[0] - xyz_2[0]) ** 2
                + (xyz_1[1] - xyz_2[1]) ** 2
                + (xyz_1[2] - xyz_2[2]) ** 2
            )
            ** 0.5
        )
        if diff > tolerances["abs_coord_diff"]:
            logging.error(
                f"PDB atom coordinates have difference > tolerance ({tolerances['abs_coord_diff']} Å"
            )
            return False
    return True


def make_double_diff_map_and_get_peaks(
    map_above,
    map_below,
    working_directory,
    pdb_file,
    map_out,
    rmsd_threshold,
    max_peaks,
):
    """Creates and calls a script in coot to generate a double difference map from the anomalous maps from above and below
    the metal absorption edge. Any peaks above the rmsd_threshold are then found and the coordinates and peak heights in
    units of rmsd and electron density are returned

    """
    coot_script = [
        "#!/usr/bin/env coot",
        "# python script for coot - generated by metal_ID",
        "set_nomenclature_errors_on_read('ignore')",
        f"read_pdb('{pdb_file}')",
        f"map_above = read_phs_and_make_map_using_cell_symm_from_previous_mol('{map_above}')",
        f"map_below = read_phs_and_make_map_using_cell_symm_from_previous_mol('{map_below}')",
        "map_diff = difference_map(map_above, map_below, 1)",
        f"difference_map_peaks(3, 0, {rmsd_threshold}, 0.0, 1, 0, 0)",
        f"export_map(map_diff, '{map_out}')",
        "coot_real_exit(0)",
    ]
    coot_script_path = working_directory / "coot_diff_map.py"
    with open(coot_script_path, "w") as script_file:
        for line in coot_script:
            script_file.write(line + "\n")
    logging.info(f"Running coot script {coot_script_path} to create diff.map")
    coot_command = f"coot --no-guano --no-graphics -s {coot_script_path}"
    result = subprocess.run(
        coot_command,
        shell=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
    )

    with open(working_directory / "metal_id.log", "w") as log_file:
        log_file.write(result.stdout)

    logging.info("\n## Finding peaks in double difference map ##")
    # Regex pattern to match lines containing peaks from coot output in format: "0 dv: 77.94 n-rmsd: 42.52 xyz = (     24.08,     12.31,     28.48)"
    pattern = (
        r"\s*\d+\s+dv:\s*([\d.]+)\s+n-rmsd:\s*([\d.]+)\s+xyz\s*=\s*\(\s*([\d., -]+)\)"
    )
    # Extract peaks from coot output
    matches = re.finditer(pattern, result.stdout)
    electron_densities = []
    rmsds = []
    peak_coords = []
    for match in matches:
        if len(peak_coords) == max_peaks:
            logging.warning(
                f"Found more peaks than the set maximum of {max_peaks} - storing only the largest {max_peaks}"
            )
            break
        density = float(match.group(1))
        rmsd = float(match.group(2))
        xyz = tuple(map(float, match.group(3).split(",")))
        electron_densities.append(density)
        rmsds.append(rmsd)
        peak_coords.append(xyz)

    return peak_coords, electron_densities, rmsds


def find_protein_centre(pdb_file):
    """Runs find-blobs on a pdb file then uses regex to get the protein centre of mass coordinates from the output
    Despite the name, find-blobs is not being used to find blobs here. Returns (x,y,z) coordinates of the protein centre.
    """
    find_blobs_command = f"find-blobs -c {pdb_file}"
    result = subprocess.run(
        find_blobs_command,
        shell=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
    )
    # Regex pattern for extracting coords from find-blobs output in format "Protein mass center: xyz = (     12.37,     23.89,     32.69)"
    pattern = r"Protein mass center: xyz = \(\s*([-+]?\d*\.\d+|\d+\.\d*)\s*,\s*([-+]?\d*\.\d+|\d+\.\d*)\s*,\s*([-+]?\d*\.\d+|\d+\.\d*)\s*\)"
    match = re.search(pattern, result.stdout)
    assert match, "Protein mass center not found"
    centre = tuple(map(float, match.groups()))
    return centre


def render_diff_map_peaks(
    output_directory,
    pdb_file,
    diff_map,
    peak_threshold,
    peak_coords,
):
    """Plots protein molecule coordinates (pdb file) and difference map in coot, applied a threshold for displaying the map
    then renders images centred on the peak_coords, viewing towards the protein centre.
    """
    logging.info("Finding protein centre")
    protein_centre = find_protein_centre(pdb_file)

    logging.info(f"Protein mass centre at: {protein_centre}")

    render_script = [
        "#!/usr/bin/env coot",
        "# python script for coot - generated by metal_ID",
        "set_nomenclature_errors_on_read('ignore')",
        f"read_pdb('{pdb_file}')",
        f"read_ccp4_map('{diff_map}', 1)",
        f"set_contour_level_in_sigma(1, {peak_threshold})",
    ]
    render_paths = []
    for _i, peak in enumerate(peak_coords, start=1):
        quat = view_as_quat(peak, protein_centre)
        # Use relative path as explicit paths can exceed render command length limit
        render_path = str(output_directory / f"peak_{_i}.r3d")
        mini_script = [
            f"set_rotation_centre{peak}",
            "set_zoom(30.0)",
            f"set_view_quaternion{quat}",
            "graphics_draw()",
            f"raster3d('{str(render_path)}')",
        ]
        render_script.extend(mini_script)
        render_paths.append(render_path)
    render_script.append("coot_real_exit(0)")

    render_script_path = output_directory / "coot_render.py"
    with open(render_script_path, "w") as script_file:
        for line in render_script:
            script_file.write(line + "\n")
    logging.info(f"Running coot rendering script {render_script_path}")
    render_command = f"coot --no-guano --no-graphics -s {render_script_path}"
    subprocess.run(
        render_command,
        shell=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
    )

    # Convert r3d files to pngs
    logging.info("Converting r3d files to pngs")
    for render_path in render_paths:
        render_png_path = f"{os.path.splitext(render_path)[0]}.png"
        logging.info(f"Converting {render_path} to {render_png_path}")
        r3d_command = f"cat {render_path} | render -png {render_png_path}"
        subprocess.run(
            r3d_command,
            shell=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
        )


def calc_double_diff_maps(
    pdb_above, pdb_below, pha_above, pha_below, output_dir, peak_threshold, max_peaks
):
    pdb_files = [pdb_above, pdb_below]

    # Check file inputs
    for file_type, file_path in [
        ("AnoDe map above", pha_above),
        ("AnoDe map below", pha_below),
    ]:
        if not file_path.is_file():
            logging.error(f"Could not find {file_type}, expected at: {file_path}")
            return False

    # Check pdb files
    logging.info(
        f"Checking pdb files for similarity. Files: {pdb_files[0]}, {pdb_files[1]}"
    )
    pdbs_are_similar = are_pdbs_similar(
        pdb_files[0],
        pdb_files[1],
    )

    if not pdbs_are_similar:
        logging.error("PDB files are not similar enough, not running metal_id")
        return False
    logging.info("PDB files are similar enough, continuing with metal_id")
    pdb_file = pdb_files[0]

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

    peak_coords, electron_densities, rmsds = make_double_diff_map_and_get_peaks(
        pha_above,
        pha_below,
        output_dir,
        pdb_file,
        map_out,
        peak_threshold,
        max_peaks,
    )

    peak_data = list(zip(electron_densities, rmsds, peak_coords))

    # Print the extracted information
    logging.info(
        f"\nThe largest peaks (up to a maximum of {max_peaks} peaks) found above the threshold of {peak_threshold} rmsd:"
    )

    peak_file = output_dir / "found_peaks.dat"

    with open(peak_file, "w") as fh:
        for peak_num, (density, rmsd, xyz) in enumerate(peak_data, start=1):
            line = f"Peak {peak_num}: Electron Density = {density} e/Å^3, RMSD = {rmsd}, XYZ = {xyz}"
            logging.info(line)
            fh.write(line + "\n")

    logging.info("\n## Rendering images of peaks ##\n")
    render_diff_map_peaks(output_dir, pdb_file, map_out, peak_threshold, peak_coords)
