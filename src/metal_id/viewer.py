from __future__ import annotations

from pathlib import Path

import gemmi
import molviewspec as mvs
import numpy as np
import os
import pandas as pd
import shutil
import re


def parse_peaks(peak_file):
    atoms = []

    # Look for matches in format:
    # Peak 1: Electron Density = 77.93 e/Å^3, RMSD = 42.52, XYZ = (24.08, 12.31, 28.48)
    pattern = re.compile(
        r"Peak\s+(?P<peak>\d+):\s*"
        r"Electron Density\s*=\s*(?P<density>[-\d.]+)\s*e/Å\^3,\s*"
        r"RMSD\s*=\s*(?P<rmsd>[-\d.]+),\s*"
        r"XYZ\s*=\s*\((?P<x>[-\d.]+),\s*(?P<y>[-\d.]+),\s*(?P<z>[-\d.]+)\)"
    )

    with open(peak_file, "r") as f:
        for line in f:
            match = pattern.search(line)
            if match:
                xf, yf, zf = map(
                    float, (match.group("x"), match.group("y"), match.group("z"))
                )
                height = float(match.group("rmsd"))
                atoms.append([xf, yf, zf, height])

    atoms_df = pd.DataFrame(atoms, columns=["x", "y", "z", "Height"])

    return atoms_df


def save_cropped_maps(
    pdb_file, map_file, atoms_df, peaks, radius, prefix, results_directory
):
    st = gemmi.read_structure(pdb_file.as_posix())
    cell = st.cell

    tmpdir = Path(results_directory) / "tmp_molviewspec"
    tmpdir.mkdir(parents=False, exist_ok=True)

    for i in range(peaks):
        m = gemmi.read_ccp4_map(map_file.as_posix(), setup=True)
        grid = m.grid

        mask = grid.clone()
        mask.fill(0.0)

        map_out = f"{prefix}{i + 1}.map"  #
        center = gemmi.Position(
            atoms_df["x"][i], atoms_df["y"][i], atoms_df["z"][i]
        )  # peak loc

        mask.set_points_around(center, radius, 1.0, use_pbc=True)  # spherical mask in Å

        dl = gemmi.Position(radius, radius, radius)  # box d/2
        box = gemmi.FractionalBox()
        box.extend(cell.fractionalize(center - dl))
        box.extend(cell.fractionalize(center + dl))

        grid.array[:] *= mask.array
        m.set_extent(box)
        outfile = str(f"{tmpdir}/{map_out}")
        m.write_ccp4_map(outfile)


def find_camera_pos(structure: gemmi.Structure):
    atoms = [
        atom
        for model in structure
        for chain in model
        for residue in chain
        for atom in residue
    ]
    a, b, c, d = gemmi.find_best_plane(atoms)
    # normal = gemmi.Vec3(a, b, c)

    targetp = gemmi.Position(0, 0, 0)
    # The center can be taken as the centroid of the input atoms projected onto the plane
    for atom in atoms:
        targetp += atom.pos
    targetp /= len(atoms)
    targetp

    camera_pos = targetp - (120 * gemmi.Position(a, b, c))

    return targetp.tolist(), camera_pos.tolist()


def mtz_to_map(mtz_file, map_file, label="FWT", ph_label="PHWT"):
    # Load MTZ
    mtz = gemmi.read_mtz_file(mtz_file.as_posix())
    grid = mtz.transform_f_phi_to_map(label, ph_label, sample_rate=3)

    # Write CCP4 map
    ccp4 = gemmi.Ccp4Map()
    ccp4.grid = grid
    ccp4.update_ccp4_header()
    ccp4.write_ccp4_map(map_file.as_posix())


def gen_html_metal_id(results_directory, isovalue=5):
    pdb_file = results_directory / "final.pdb"
    mtz_file = results_directory / "dimple_below.mtz"
    diff_map = results_directory / "diff.map"
    peak_file = results_directory / "found_peaks.dat"

    atoms_df = parse_peaks(peak_file)
    peaks = len(atoms_df)  # set peaks to length of atoms_df
    heights = atoms_df["Height"]

    # convert dimple mtz to map format
    map_file = results_directory / "final.map"
    mtz_to_map(mtz_file, map_file)

    save_cropped_maps(
        pdb_file,
        diff_map,
        atoms_df,
        peaks=peaks,
        radius=3,
        prefix="box",
        results_directory=results_directory,
    )

    save_cropped_maps(
        pdb_file,
        map_file,
        atoms_df,
        peaks=peaks,
        radius=6,
        prefix="fbox",
        results_directory=results_directory,
    )  # 2fofc

    st = gemmi.read_structure(pdb_file.as_posix())
    cell = st.cell

    snapshot_list = []

    for j in range(peaks):
        globals()["map_file" + str(j + 1)] = (
            f"{results_directory}/tmp_molviewspec/box{j + 1}.map"
        )
        globals()["fmap_file" + str(j + 1)] = (
            f"{results_directory}/tmp_molviewspec/fbox{j + 1}.map"
        )

        with open(globals()["map_file" + str(j + 1)], mode="rb") as f:
            globals()["map_data" + str(j + 1)] = f.read()

        with open(globals()["fmap_file" + str(j + 1)], mode="rb") as f:
            globals()["fmap_data" + str(j + 1)] = f.read()

    targetp, camerap = find_camera_pos(st)
    ns = (
        gemmi.NeighborSearch(st[0], cell, 5).populate(include_h=False).populate()
    )  # neighbour search

    for i in range(peaks + 1):
        if i == 0:  # main page is different to all the others
            builder = mvs.create_builder()

            structure = (
                builder.download(url=pdb_file).parse(format="pdb").model_structure()
            )  # symmetry_mates_structure()
            structure.component(selector="polymer").representation(
                type="surface", size_factor=0.9
            ).opacity(opacity=0.2).color(color="#AABDF1")
            structure.component(selector="polymer").representation().opacity(
                opacity=0.25
            ).color(custom={"molstar_color_theme_name": "chain_id"})
            structure.component(selector="ligand").representation(
                type="ball_and_stick"
            ).color(custom={"molstar_color_theme_name": "element-symbol"})
            structure.component(selector="ligand").representation(
                type="surface"
            ).opacity(opacity=0.1).color(
                custom={"molstar_color_theme_name": "element-symbol"}
            )

            for j in range(peaks):
                peakcoords = np.array(
                    [atoms_df["x"][j], atoms_df["y"][j], atoms_df["z"][j]]
                ).tolist()
                labelcoords = np.array(
                    [atoms_df["x"][j], atoms_df["y"][j], atoms_df["z"][j]]
                ).tolist()
                builder.primitives(opacity=0.1).sphere(
                    center=peakcoords,
                    radius=1,
                    color="#da21fa",
                    tooltip=f"peak {j + 1}",
                ).label(position=labelcoords, text=f"{j + 1}", label_size=5)

            for k in range(peaks):
                ccp4 = builder.download(url=globals()["map_file" + str(k + 1)]).parse(
                    format="map"
                )
                ccp4.volume().representation(
                    type="isosurface",
                    relative_isovalue=isovalue,
                    show_wireframe=True,
                    show_faces=False,
                ).color(color="#da21fa").opacity(opacity=0.25)

            builder.camera(position=camerap, target=targetp, up=[0, 0, 1])

            globals()["snapshot" + str(i + 1)] = builder.get_snapshot(
                title="Main View",
                description=f"## Metal_ID Results: \n ### Summary \n - Anomalous double difference map shown at {isovalue}σ, magenta for the top {peaks} sites listed in 'found_peaks.dat'",
                transition_duration_ms=700,
                linger_duration_ms=5000,
                key="Main",
            )

            snapshot_list.append(globals()["snapshot" + str(i + 1)])

        else:
            builder = mvs.create_builder()

            structure = (
                builder.download(url=pdb_file).parse(format="pdb").model_structure()
            )  # symmetry_mates_structure()
            structure.component(selector="polymer").representation(
                type="surface", size_factor=0.9
            ).opacity(opacity=0.2).color(color="#AABDF1")
            structure.component(selector="polymer").representation().opacity(
                opacity=0.25
            ).color(custom={"molstar_color_theme_name": "chain_id"})
            structure.component(selector="ligand").representation(
                type="ball_and_stick"
            ).color(custom={"molstar_color_theme_name": "element-symbol"})
            structure.component(selector="ligand").representation(
                type="surface"
            ).opacity(opacity=0.1).color(
                custom={"molstar_color_theme_name": "element-symbol"}
            )

            for j in range(peaks):
                peakcoords = np.array(
                    [atoms_df["x"][j], atoms_df["y"][j], atoms_df["z"][j]]
                ).tolist()
                labelcoords = np.array(
                    [atoms_df["x"][j], atoms_df["y"][j], atoms_df["z"][j]]
                ).tolist()
                builder.primitives(opacity=0.1).sphere(
                    center=peakcoords,
                    radius=1,
                    color="#da21fa",
                    tooltip=f"peak {j + 1}",
                ).label(position=labelcoords, text=f"{j + 1}", label_size=2)  # spheres

            nearest_atom_mark = ns.find_nearest_atom(
                gemmi.Position(
                    atoms_df["x"][i - 1], atoms_df["y"][i - 1], atoms_df["z"][i - 1]
                )
            )
            residue = mvs.ComponentExpression(
                atom_id=st[0][nearest_atom_mark.chain_idx][
                    nearest_atom_mark.residue_idx
                ][nearest_atom_mark.atom_idx].serial
            )  # nearest atom id
            structure.component(
                selector=residue,
                custom={
                    "molstar_show_non_covalent_interactions": True,
                    "molstar_non_covalent_interactions_radius_ang": 5.0,
                },
            )

            for k in range(peaks):
                ccp4 = builder.download(url=globals()["map_file" + str(k + 1)]).parse(
                    format="map"
                )
                isosurface_peak = (
                    ccp4.volume()
                    .representation(
                        type="isosurface",
                        relative_isovalue=isovalue,
                        show_wireframe=True,
                        show_faces=False,
                    )
                    .color(color="#da21fa")
                    .opacity(opacity=0.25)
                )
                if k == i - 1:
                    isosurface_peak.focus()  # focus on the current peak site

                ccp42 = builder.download(url=globals()["fmap_file" + str(k + 1)]).parse(
                    format="map"
                )  # 2fo-fc
                ccp42.volume().representation(
                    type="isosurface",
                    relative_isovalue=1.5,
                    show_wireframe=True,
                    show_faces=False,
                ).color(color="#2f78d7").opacity(opacity=0.25)

            globals()["snapshot" + str(i + 1)] = builder.get_snapshot(
                title=f"Site {i}",
                description=f"## Anode Results: \n ### Site {i} \n - Displaying unique site {i}, height {heights[i - 1]} σ \n - Anomolous double difference map {isovalue}σ, magenta \n - 2FO-FC at 1.5σ, blue \n \n [Back to Main Summary Page](#Main)",
                transition_duration_ms=700,
                linger_duration_ms=5000,
                key="site1",
            )

            snapshot_list.append(globals()["snapshot" + str(i + 1)])

    with open(pdb_file) as f:
        pdb_data = f.read()

    data_dict = {}
    for i in range(peaks):
        key = globals()["map_file" + str(i + 1)]
        value = globals()["map_data" + str(i + 1)]
        data_dict[key] = value

        key = globals()["fmap_file" + str(i + 1)]
        value = globals()["fmap_data" + str(i + 1)]
        data_dict[key] = value

    data_dict[pdb_file.as_posix()] = pdb_data

    states = mvs.States(
        snapshots=snapshot_list, metadata=mvs.GlobalMetadata(description="anode")
    )
    html = mvs.molstar_widgets.molstar_html(states, data=data_dict, ui="stories")
    with open(f"{results_directory}/metal_id.html", "w") as f:
        f.write(html)

    # clean up
    tmpdir = results_directory / "tmp_molviewspec"
    shutil.rmtree(str(tmpdir))
    os.remove(map_file)


gen_html_metal_id(
    results_directory=Path(
        "/scratch/dwe15129_scratch_space/metal_ID/2023_test_data/processed/metal_id"
    ),
    isovalue=5,
)
