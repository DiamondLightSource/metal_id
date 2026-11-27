"""
Tests of the whole workflow with the real CCP4 programs on real data: thermolysin with
data above and below the Ca K edge (three Ca sites with anomalous signal).

The input files are edited at test time to send metal_id down its different routes.

These need CCP4 and the data files, and are skipped without them. The data are read
from tests/data, or from the directory in METAL_ID_TEST_DATA:

    above.mtz   data above the Ca K edge (IMEAN, I(+), I(-))
    below.mtz   data below the Ca K edge
    model.pdb   model refined against these data (i.e. dimple's final.pdb)
    2TLX.pdb    the deposited thermolysin structure (its cell differs by 2-3%)

To run only these:

    module load ccp4
    pytest tests/test_real_data.py
"""

import os
from pathlib import Path

import gemmi
import pytest

from metal_id.metal_id import main
from metal_id.viewer import parse_peaks

pytestmark = pytest.mark.ccp4

DATA_DIR = Path(os.environ.get("METAL_ID_TEST_DATA", Path(__file__).parent / "data"))
DATA_FILES = ["above.mtz", "below.mtz", "model.pdb", "2TLX.pdb"]
# Number of Ca sites with anomalous signal (Ca323, Ca325, Ca326; Ca324 has none)
N_CA_SITES = 3
MIN_CA_PEAK_RMSD = 12.0


@pytest.fixture(scope="module")
def data() -> dict[str, Path]:
    missing = [name for name in DATA_FILES if not (DATA_DIR / name).is_file()]
    if missing:
        pytest.skip(
            f"test data not found in {DATA_DIR} ({', '.join(missing)}) - "
            "set METAL_ID_TEST_DATA to the data directory"
        )
    return {Path(name).stem: DATA_DIR / name for name in DATA_FILES}


def run(out: Path, above: Path, below: Path, model: Path, *options: str):
    main([str(above), str(below), str(model), "-o", str(out), *options])


def peaks_near(out: Path, element: str, max_distance: float = 1.0):
    """
    The peaks found by metal_id, each with the distance to the nearest atom of the
    element in the model metal_id used (so this works whatever frame that model is in)
    """
    structure = gemmi.read_structure(str(out / "final.pdb"))
    sites = [
        cra.atom.pos for cra in structure[0].all() if cra.atom.element.name == element
    ]
    results = []
    for peak in parse_peaks(out / "found_peaks.dat"):
        pos = gemmi.Position(peak.x, peak.y, peak.z)
        distance = min(
            structure.cell.find_nearest_image(site, pos, gemmi.Asu.Any).dist()
            for site in sites
        )
        results.append((peak, distance <= max_distance))
    return results


def assert_ca_sites_found(out: Path):
    """The largest peaks are at the Ca sites, and none is at Zn or a Met/DMSO S"""
    peaks = peaks_near(out, "Ca")
    top = peaks[:N_CA_SITES]
    assert len(top) == N_CA_SITES
    assert all(near_ca for _, near_ca in top), [p.height for p, _ in peaks]
    assert all(peak.height >= MIN_CA_PEAK_RMSD for peak, _ in top)
    for element in ["Zn", "S"]:
        assert not any(near for _, near in peaks_near(out, element, 1.5)), element


def set_model_cell(model: Path, out_file: Path, cell, spacegroup: str) -> Path:
    """Copy a model with a different cell and space group (as from another crystal form)"""
    structure = gemmi.read_structure(str(model))
    structure.cell = gemmi.UnitCell(*cell)
    structure.spacegroup_hm = spacegroup
    structure.write_pdb(str(out_file))
    return out_file


def scale_mtz_cell(mtz_file: Path, out_file: Path, factor: float) -> Path:
    """Copy an mtz file with the cell lengths multiplied by factor"""
    mtz = gemmi.read_mtz_file(str(mtz_file))
    a, b, c, alpha, beta, gamma = mtz.cell.parameters
    mtz.set_cell_for_all(
        gemmi.UnitCell(a * factor, b * factor, c * factor, alpha, beta, gamma)
    )
    mtz.write_to_file(str(out_file))
    return out_file


def truncate_model(model: Path, out_file: Path, fraction: float) -> Path:
    """Copy a model keeping only the first fraction of the residues in each chain"""
    structure = gemmi.read_structure(str(model))
    for chain in structure[0]:
        keep = int(len(chain) * fraction)
        for i in reversed(range(keep, len(chain))):
            del chain[i]
    structure.write_pdb(str(out_file))
    return out_file


def test_refmac(tmp_path, data):
    out = tmp_path / "out"
    run(out, data["above"], data["below"], data["model"])
    assert_ca_sites_found(out)
    assert not (out / "dimple").exists()


def test_anode(tmp_path, data):
    out = tmp_path / "out"
    run(out, data["above"], data["below"], data["model"], "--method", "anode")
    assert_ca_sites_found(out)


def test_deposited_model_runs_dimple(tmp_path, data):
    # The cell of 2TLX differs from the data's by 2-3%, so dimple refines it first
    out = tmp_path / "out"
    run(out, data["above"], data["below"], data["2TLX"])
    assert (out / "dimple" / "final.pdb").is_file()
    assert_ca_sites_found(out)


def test_other_crystal_form_runs_molecular_replacement(tmp_path, data):
    # A model from another crystal form doesn't match the data, and dimple goes
    # straight to MR. The results may be in a different frame from the input model.
    # (A model without CRYST1 would be the same test, but dimple in CCP4 9.0.017 fails
    # on one at its rwcontents step - a bug fixed in the next CCP4 release.)
    out = tmp_path / "out"
    model = set_model_cell(
        data["model"], tmp_path / "other_form.pdb", (60, 70, 120, 90, 90, 90), "P 1"
    )
    run(out, data["above"], data["below"], model)
    assert (out / "dimple" / "final.pdb").is_file()
    assert_ca_sites_found(out)


def test_swapped_datasets_find_no_ca(tmp_path, data):
    # With above and below swapped, the Ca signal is negative in the double difference
    out = tmp_path / "out"
    run(out, data["below"], data["above"], data["model"])
    if (out / "found_peaks.dat").exists():
        assert not any(near for _, near in peaks_near(out, "Ca"))


def test_incomplete_model_fails_r_factor_check(tmp_path, data):
    # A model with most of the protein missing fits the data badly
    model = truncate_model(data["model"], tmp_path / "partial.pdb", 0.3)
    with pytest.raises(SystemExit):
        run(tmp_path / "out", data["above"], data["below"], model)
    assert not (tmp_path / "out" / "found_peaks.dat").exists()


def test_incompatible_cells(tmp_path, data):
    # The below data with a cell 3% larger don't match the above data
    below = scale_mtz_cell(data["below"], tmp_path / "below_big.mtz", 1.03)
    with pytest.raises(SystemExit):
        run(tmp_path / "out", data["above"], below, data["model"])
    assert not (tmp_path / "out" / "above").exists()
