from __future__ import annotations

from collections import namedtuple
from pathlib import Path

import gemmi
import molviewspec as mvs
import os
import shutil
import re

Peak = namedtuple("Peak", ["x", "y", "z", "height"])


def parse_peaks(peak_file: Path) -> list[Peak]:
    peaks = []

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
                peaks.append(Peak(xf, yf, zf, height))

    return peaks


def save_cropped_maps(
    pdb_file: Path,
    map_file: Path,
    peak: Peak,
    radius: float,
    filename: str,
    tmpdir: Path,
):
    structure = gemmi.read_structure(pdb_file.as_posix())
    cell = structure.cell

    map = gemmi.read_ccp4_map(map_file.as_posix(), setup=True)
    grid = map.grid

    mask = grid.clone()
    mask.fill(0.0)

    map_out = f"{filename}.map"  #
    center = gemmi.Position(peak.x, peak.y, peak.z)

    mask.set_points_around(center, radius, 1.0, use_pbc=True)  # spherical mask in Å

    dl = gemmi.Position(radius, radius, radius)  # box d/2
    box = gemmi.FractionalBox()
    box.extend(cell.fractionalize(center - dl))
    box.extend(cell.fractionalize(center + dl))

    grid.array[:] *= mask.array
    map.set_extent(box)
    outfile = str(f"{tmpdir}/{map_out}")
    map.write_ccp4_map(outfile)
    with open(outfile, "rb") as f:
        map_data = f.read()
    return map_data


def find_camera_pos(structure: gemmi.Structure):
    atoms = [
        atom
        for model in structure
        for chain in model
        for residue in chain
        for atom in residue
    ]
    a, b, c, d = gemmi.find_best_plane(atoms)

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


def generate_spheres(builder: mvs.Builder, peaks: list[Peak]):
    for peak_num, peak in enumerate(peaks, start=1):
        peakcoords = [peak.x, peak.y, peak.z]
        builder.primitives(opacity=0.1).sphere(
            center=peakcoords,
            radius=1,
            color="#da21fa",
            tooltip=f"peak {peak_num}",
        ).label(position=peakcoords, text=f"{peak_num}", label_size=2)


def generate_isosurfaces_and_focus_on_current_peak(
    builder: mvs.Builder,
    peaks: list[Peak],
    focus_peak_num: int,
    isovalue=5,
):
    for peak_num, _ in enumerate(peaks, start=1):
        # Add anomalous double difference map for the current peak
        ccp4 = builder.download(url=f"map{peak_num}").parse(format="map")
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
        if peak_num == focus_peak_num:
            isosurface_peak.focus()  # focus on the current peak site

        # Add 2FO-FC map for the current peak
        ccp4_2fo_fc = builder.download(url=f"fmap{peak_num}").parse(format="map")
        ccp4_2fo_fc.volume().representation(
            type="isosurface",
            relative_isovalue=1.5,
            show_wireframe=True,
            show_faces=False,
        ).color(color="#2f78d7").opacity(opacity=0.25)


def gen_html_metal_id(
    results_directory: Path,
    pdb_file: Path,
    mtz_file: Path,
    diff_map_file: Path,
    peak_file: Path,
    isovalue: float = 5.0,
):
    peaks = parse_peaks(peak_file)

    # convert dimple mtz to map format
    map_file = results_directory / "final.map"
    mtz_to_map(mtz_file, map_file)

    tmpdir = Path(results_directory) / "tmp_molviewspec"
    tmpdir.mkdir(parents=False, exist_ok=True)

    data_dict = {}
    for peak_num, peak in enumerate(peaks, start=1):
        # Anomalous double difference map
        map_data = save_cropped_maps(
            pdb_file,
            diff_map_file,
            peak,
            radius=3,
            filename=f"box{peak_num}",
            tmpdir=tmpdir,
        )
        # 2FO-FC map
        fmap_data = save_cropped_maps(
            pdb_file,
            map_file,
            peak,
            radius=6,
            filename=f"fbox{peak_num}",
            tmpdir=tmpdir,
        )
        data_dict[f"map{peak_num}"] = map_data
        data_dict[f"fmap{peak_num}"] = fmap_data

    st = gemmi.read_structure(pdb_file.as_posix())
    cell = st.cell

    # Store map files and data in dictionaries indexed by peak number
    snapshot_list = []

    targetp, camerap = find_camera_pos(st)
    ns = gemmi.NeighborSearch(st[0], cell, 5).populate(include_h=False).populate()

    """ Create main page """
    builder = mvs.create_builder()

    structure = (
        builder.download(url="pdb").parse(format="pdb").model_structure()
    )  # symmetry_mates_structure()
    structure.component(selector="polymer").representation(
        type="surface", size_factor=0.9
    ).opacity(opacity=0.2).color(color="#AABDF1")
    structure.component(selector="polymer").representation().opacity(
        opacity=0.25
    ).color(custom={"molstar_color_theme_name": "chain_id"})
    structure.component(selector="ligand").representation(type="ball_and_stick").color(
        custom={"molstar_color_theme_name": "element-symbol"}
    )
    structure.component(selector="ligand").representation(type="surface").opacity(
        opacity=0.1
    ).color(custom={"molstar_color_theme_name": "element-symbol"})

    for peak_num, peak in enumerate(peaks, start=1):
        peakcoords = [peak.x, peak.y, peak.z]
        builder.primitives(opacity=0.1).sphere(
            center=peakcoords,
            radius=1,
            color="#da21fa",
            tooltip=f"peak {peak_num}",
        ).label(position=peakcoords, text=f"{peak_num}", label_size=5)

        ccp4 = builder.download(url=f"map{peak_num}").parse(format="map")
        ccp4.volume().representation(
            type="isosurface",
            relative_isovalue=isovalue,
            show_wireframe=True,
            show_faces=False,
        ).color(color="#da21fa").opacity(opacity=0.25)

    builder.camera(position=camerap, target=targetp, up=[0, 0, 1])

    snapshot_main = builder.get_snapshot(
        title="Main View",
        description=f"## Metal_ID Results: \n ### Summary \n - Anomalous double difference map shown at {isovalue}σ, magenta for the top {len(peaks)} sites listed in 'found_peaks.dat'",
        transition_duration_ms=700,
        linger_duration_ms=5000,
        key="Main",
    )
    snapshot_list.append(snapshot_main)

    """ Create individual peak pages """
    for peak_num, peak in enumerate(peaks, start=1):
        builder = mvs.create_builder()

        structure = (
            builder.download(url="pdb").parse(format="pdb").model_structure()
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
        structure.component(selector="ligand").representation(type="surface").opacity(
            opacity=0.1
        ).color(custom={"molstar_color_theme_name": "element-symbol"})

        generate_spheres(builder, peaks)
        nearest_atom_mark = ns.find_nearest_atom(gemmi.Position(peak.x, peak.y, peak.z))
        residue = mvs.ComponentExpression(
            atom_id=st[0][nearest_atom_mark.chain_idx][nearest_atom_mark.residue_idx][
                nearest_atom_mark.atom_idx
            ].serial
        )  # nearest atom id
        structure.component(
            selector=residue,
            custom={
                "molstar_show_non_covalent_interactions": True,
                "molstar_non_covalent_interactions_radius_ang": 5.0,
            },
        )

        generate_isosurfaces_and_focus_on_current_peak(
            builder, peaks, peak_num, isovalue=isovalue
        )

        snapshot = builder.get_snapshot(
            title=f"Site {peak_num}",
            description=f"## Metal_ID Results: \n ### Site {peak_num} \n - Displaying unique site {peak_num}, height {peak.height} σ \n - Anomalous double difference map {isovalue}σ, magenta \n - 2FO-FC at 1.5σ, blue \n \n [Back to Main Summary Page](#Main)",
            transition_duration_ms=700,
            linger_duration_ms=5000,
            key="site1",
        )

        snapshot_list.append(snapshot)

    with open(pdb_file) as f:
        pdb_data = f.read()

    data_dict["pdb"] = pdb_data

    states = mvs.States(
        snapshots=snapshot_list, metadata=mvs.GlobalMetadata(description="metal_id")
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
    pdb_file=Path(
        "/scratch/dwe15129_scratch_space/metal_ID/2023_test_data/processed/metal_id/final.pdb"
    ),
    mtz_file=Path(
        "/scratch/dwe15129_scratch_space/metal_ID/2023_test_data/processed/metal_id/dimple_below.mtz"
    ),
    diff_map_file=Path(
        "/scratch/dwe15129_scratch_space/metal_ID/2023_test_data/processed/metal_id/diff.map"
    ),
    peak_file=Path(
        "/scratch/dwe15129_scratch_space/metal_ID/2023_test_data/processed/metal_id/found_peaks.dat"
    ),
)
