"""Synthetic structures, mtz files and stand-ins for the CCP4 programs, for testing
metal_id without CCP4.

The test structure is a cloud of water O atoms with one Ca, one Zn and one S. The stand-in
programs calculate their anomalous map coefficients from the model, weighting each
anomalous scatterer by its tabulated f'' at the wavelength they are given, so a
double difference map from them has a peak at the Ca site, whose f'' changes across
the Ca K edge between the two test wavelengths.
"""

import shutil
from pathlib import Path

import gemmi
import numpy as np

CELL = (40.0, 45.0, 50.0, 90.0, 90.0, 90.0)
SPACEGROUP = "P 21 21 21"
# Either side of the Ca K edge (3.07 Å)
WAVELENGTH_ABOVE = 3.024
WAVELENGTH_BELOW = 3.17908
# Placeholder wavelength that some programs write for the HKL_base dataset
BASE_WAVELENGTH = 9.999

CA_SITE = (15.0, 20.0, 25.0)
ZN_SITE = (26.0, 14.0, 30.0)
S_SITE = (24.0, 27.0, 18.0)
ANOMALOUS_ELEMENTS = ("Ca", "Zn", "S")


def make_atom(name: str, element: str, xyz, b_iso: float = 20.0) -> gemmi.Atom:
    atom = gemmi.Atom()
    atom.name = name
    atom.element = gemmi.Element(element)
    atom.pos = gemmi.Position(*xyz)
    atom.occ = 1.0
    atom.b_iso = b_iso
    return atom


def make_residue(name: str, seq_num: int, atoms: list[gemmi.Atom], het: bool):
    residue = gemmi.Residue()
    residue.name = name
    residue.seqid = gemmi.SeqId(seq_num, " ")
    residue.het_flag = "H" if het else "A"
    for atom in atoms:
        residue.add_atom(atom)
    return residue


def make_structure(cell=CELL, spacegroup=SPACEGROUP, n_waters: int = 100):
    """A structure with n_waters water O atoms at random positions near the middle of
    the cell, plus a Met SD, a Ca ion and a Zn ion at fixed sites. Waters, rather than
    residues with bonds between their atoms, so that refinement (e.g. by dimple) has no
    bond restraints to pull the random positions together."""
    rng = np.random.default_rng(0)
    sites = [np.array(site) for site in (CA_SITE, ZN_SITE, S_SITE)]
    placed = []
    chain = gemmi.Chain("A")
    for _ in range(n_waters):
        # Keep atoms apart, and away from the anomalous scatterers
        while True:
            xyz = rng.uniform(8, 32, size=3)
            if all(np.linalg.norm(xyz - site) > 3.5 for site in sites) and all(
                np.linalg.norm(xyz - other) > 2.8 for other in placed
            ):
                break
        placed.append(xyz)
    chain.add_residue(make_residue("MET", 1, [make_atom("SD", "S", S_SITE)], het=False))
    chain.add_residue(make_residue("CA", 2, [make_atom("CA", "Ca", CA_SITE)], het=True))
    chain.add_residue(make_residue("ZN", 3, [make_atom("ZN", "Zn", ZN_SITE)], het=True))
    for seq_num, xyz in enumerate(placed, start=4):
        chain.add_residue(
            make_residue("HOH", seq_num, [make_atom("O", "O", xyz)], het=True)
        )

    structure = gemmi.Structure()
    structure.cell = gemmi.UnitCell(*cell)
    structure.spacegroup_hm = spacegroup
    model = gemmi.Model("1")
    model.add_chain(chain)
    structure.add_model(model)
    structure.setup_entities()
    return structure


def write_model(path: Path, cell=CELL, spacegroup=SPACEGROUP) -> Path:
    make_structure(cell, spacegroup).write_pdb(str(path))
    return path


def structure_factors(
    structure: gemmi.Structure, d_min: float, wavelength: float | None = None
) -> tuple[np.ndarray, np.ndarray]:
    """
    Calculate structure factors for the unique reflections to d_min. If a wavelength is
    given, only the anomalous scatterers are included, each with its occupancy scaled
    by f''(wavelength) / Z, so the result approximates the anomalous scattering.
    Returns the Miller indices and complex structure factors.
    """
    structure = structure.clone()
    if wavelength is not None:
        energy = 12398.42 / wavelength
        for residue in structure[0]["A"]:
            for i in reversed(range(len(residue))):
                atom = residue[i]
                if atom.element.name not in ANOMALOUS_ELEMENTS:
                    del residue[i]
                    continue
                z = atom.element.atomic_number
                atom.occ = gemmi.cromer_liberman(z=z, energy=energy)[1] / z
        structure.remove_empty_chains()
    calc = gemmi.DensityCalculatorX()
    calc.d_min = d_min
    calc.set_grid_cell_and_spacegroup(structure)
    calc.put_model_density_on_grid(structure[0])
    asu_data = gemmi.transform_map_to_f_phi(calc.grid).prepare_asu_data(dmin=d_min)
    return np.array(asu_data.miller_array), np.array(asu_data.value_array)


def write_mtz(
    path: Path,
    hkl: np.ndarray,
    columns: dict[str, tuple[str, np.ndarray]],
    wavelength: float,
    cell=CELL,
    spacegroup=SPACEGROUP,
    base_wavelength: float = BASE_WAVELENGTH,
) -> Path:
    """Write an mtz file with a base dataset and one dataset holding the columns"""
    mtz = gemmi.Mtz(with_base=True)
    mtz.spacegroup = gemmi.SpaceGroup(spacegroup)
    mtz.set_cell_for_all(gemmi.UnitCell(*cell))
    mtz.datasets[0].wavelength = base_wavelength
    dataset = mtz.add_dataset("SAD")
    dataset.wavelength = wavelength
    for label, (column_type, _) in columns.items():
        mtz.add_column(label, column_type)
    data = np.column_stack([hkl, *(values for _, values in columns.values())])
    mtz.set_data(data.astype(np.float32))
    mtz.write_to_file(str(path))
    return path


def write_data_mtz(
    path: Path,
    wavelength: float,
    d_min: float = 2.0,
    cell=CELL,
    spacegroup=SPACEGROUP,
    base_wavelength: float = BASE_WAVELENGTH,
    labels=("IMEAN", "SIGIMEAN", "I(+)", "SIGI(+)", "I(-)", "SIGI(-)"),
) -> Path:
    """Write an mtz file with merged intensity columns like those metal_id takes"""
    hkl, f = structure_factors(make_structure(cell, spacegroup), d_min)
    intensity = np.abs(f) ** 2
    column_types = {"IMEAN": "J", "I(+)": "K", "I(-)": "K"}
    columns = {
        label: (
            column_types.get(label, "Q" if label == "SIGIMEAN" else "M"),
            intensity if not label.startswith("SIG") else np.sqrt(intensity),
        )
        for label in labels
    }
    return write_mtz(path, hkl, columns, wavelength, cell, spacegroup, base_wavelength)


def anomalous_structure_factors(
    structure: gemmi.Structure, d_min: float, wavelength: float
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Calculate F(+) and F(-) for the unique reflections to d_min, with f' and f'' for
    every atom from its tabulated values at the wavelength. Returns the Miller indices
    and complex F(+) and F(-).
    """
    cell = structure.cell
    spacegroup = structure.find_spacegroup()
    hkl = gemmi.make_miller_array(cell, spacegroup, d_min)
    stol2 = cell.calculate_1_d2_array(hkl) / 4
    energy = 12398.42 / wavelength

    ops = list(spacegroup.operations())
    f_plus = np.zeros(len(hkl), dtype=complex)
    f_minus = np.zeros(len(hkl), dtype=complex)
    for cra in structure[0].all():
        atom = cra.atom
        z = atom.element.atomic_number
        fp, fpp = gemmi.cromer_liberman(z=z, energy=energy)
        f0 = np.array([atom.element.it92.calculate_sf(x) for x in stol2])
        f = (f0 + fp + 1j * fpp) * atom.occ * np.exp(-atom.b_iso * stol2)
        frac = np.array(cell.fractionalize(atom.pos).tolist())
        for op in ops:
            xyz = (
                np.array(op.rot) @ frac / gemmi.Op.DEN
                + np.array(op.tran) / gemmi.Op.DEN
            )
            phase = 2 * np.pi * (hkl @ xyz)
            f_plus += f * np.exp(1j * phase)
            f_minus += f * np.exp(-1j * phase)
    return hkl, f_plus, f_minus


def write_anomalous_data_mtz(
    path: Path,
    wavelength: float,
    d_min: float = 2.0,
    structure: gemmi.Structure | None = None,
    seed: int = 0,
) -> Path:
    """
    Write an mtz file with realistic merged anomalous intensities calculated from the
    structure (by default the test structure) at the wavelength, with 2% noise, for
    testing with the real CCP4 programs.
    """
    structure = structure or make_structure()
    hkl, f_plus, f_minus = anomalous_structure_factors(structure, d_min, wavelength)
    rng = np.random.default_rng(seed)
    columns = {}
    intensities = {}
    for sign, f in [("+", f_plus), ("-", f_minus)]:
        intensity = np.abs(f) ** 2
        sigma = 0.02 * intensity + 0.01 * intensity.mean()
        intensities[sign] = (intensity + rng.normal(0, 1, len(f)) * sigma, sigma)
    # Centric reflections have no anomalous difference
    centric = gemmi.GroupOps(
        list(structure.find_spacegroup().operations())
    ).centric_flag_array(hkl)
    i_mean = (intensities["+"][0] + intensities["-"][0]) / 2
    sig_mean = np.sqrt(intensities["+"][1] ** 2 + intensities["-"][1] ** 2) / 2
    columns["IMEAN"] = ("J", i_mean)
    columns["SIGIMEAN"] = ("Q", sig_mean)
    for sign in ("+", "-"):
        intensity, sigma = intensities[sign]
        columns[f"I({sign})"] = ("K", np.where(centric, i_mean, intensity))
        columns[f"SIGI({sign})"] = ("M", np.where(centric, sig_mean, sigma))
    return write_mtz(
        path,
        hkl,
        columns,
        wavelength,
        tuple(structure.cell.parameters),
        structure.spacegroup_hm,
    )


def mtz_wavelength(mtz_file: Path) -> float:
    mtz = gemmi.read_mtz_file(str(mtz_file))
    return mtz.column_with_label("I(+)").dataset.wavelength


def write_refmac_mtz(path: Path, model: Path, wavelength: float, d_min: float) -> Path:
    """Write an mtz file like refmac's, with FWT/PHWT and FAN/PHAN from the model"""
    structure = gemmi.read_structure(str(model))
    hkl, f_all = structure_factors(structure, d_min)
    hkl_anom, f_anom = structure_factors(structure, d_min, wavelength)
    assert np.array_equal(hkl, hkl_anom)
    columns = {
        "FWT": ("F", np.abs(f_all)),
        "PHWT": ("P", np.angle(f_all, deg=True)),
        "FAN": ("F", np.abs(f_anom)),
        "PHAN": ("P", np.angle(f_anom, deg=True)),
    }
    return write_mtz(
        path,
        hkl,
        columns,
        wavelength,
        tuple(structure.cell.parameters),
        structure.spacegroup_hm,
    )


def write_pha(path: Path, model: Path, wavelength: float, d_min: float) -> Path:
    """Write an ANODE-style .pha file (h k l F FOM phi sigF) from the model"""
    hkl, f_anom = structure_factors(gemmi.read_structure(str(model)), d_min, wavelength)
    with open(path, "w") as fh:
        for (h, k, l), value in zip(hkl, f_anom):  # noqa: E741
            phi = np.angle(value, deg=True)
            fh.write(f"{h:4d}{k:4d}{l:4d} {abs(value):10.3f} 1.000 {phi:8.2f} 1.00\n")
    return path


class FakeCCP4:
    """
    Stand-in for run_ccp4 that imitates the CCP4 programs metal_id runs, and records
    the calls. Outputs are calculated from the model (see the module docstring).
    """

    def __init__(self, r_factor: float = 0.2):
        self.r_factor = r_factor
        self.calls: list[tuple[str, list[str], list[str]]] = []

    def __call__(self, program, args, log_file, keywords=None, cwd=None):
        args = [str(arg) for arg in args]
        self.calls.append((program, args, keywords or []))
        cwd = Path(cwd) if cwd else Path.cwd()
        log = getattr(self, program.replace("5", ""))(args, keywords or [], cwd)
        Path(log_file).write_text(log)
        return log

    def calls_to(self, program: str):
        return [call for call in self.calls if call[0] == program]

    @staticmethod
    def _arg(args: list[str], name: str) -> str:
        return args[args.index(name) + 1]

    def pointless(self, args, keywords, cwd):
        shutil.copyfile(self._arg(args, "HKLIN"), self._arg(args, "HKLOUT"))
        return "fake pointless\n"

    def ctruncate(self, args, keywords, cwd):
        shutil.copyfile(self._arg(args, "-hklin"), self._arg(args, "-hklout"))
        return "fake ctruncate\n"

    def refmac(self, args, keywords, cwd):
        wavelength = float(
            next(k for k in keywords if k.startswith("ANOMALOUS WAVELENGTH")).split()[
                -1
            ]
        )
        hklin = gemmi.read_mtz_file(self._arg(args, "HKLIN"))
        write_refmac_mtz(
            Path(self._arg(args, "HKLOUT")),
            Path(self._arg(args, "XYZIN")),
            wavelength,
            hklin.resolution_high(),
        )
        shutil.copyfile(self._arg(args, "XYZIN"), self._arg(args, "XYZOUT"))
        return f"fake refmac\nOverall R factor                     =     {self.r_factor:.4f}\n"

    def mtz2sca(self, args, keywords, cwd):
        # Remember which mtz file the .sca file came from, for the fake anode
        (cwd / args[1]).write_text(args[0])
        return "fake mtz2sca\n"

    def shelxc(self, args, keywords, cwd):
        return "fake shelxc\n"

    def anode(self, args, keywords, cwd):
        mtz_file = Path((cwd / "anode.sca").read_text())
        d_min = gemmi.read_mtz_file(str(mtz_file)).resolution_high()
        write_pha(cwd / "anode.pha", cwd / "anode.pdb", mtz_wavelength(mtz_file), d_min)
        return "fake anode\n"
