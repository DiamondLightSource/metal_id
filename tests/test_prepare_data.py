import gemmi
import pytest

from metal_id.prepare_data import (
    check_mtz,
    copy_model,
    find_matching_model,
    model_symmetry,
    mtz_symmetry,
    symmetry_difference,
)
from synthetic import (
    CELL,
    WAVELENGTH_ABOVE,
    make_structure,
    write_data_mtz,
    write_model,
)


def sym(cell, spacegroup):
    return gemmi.UnitCell(*cell), gemmi.SpaceGroup(spacegroup)


@pytest.mark.parametrize(
    "sym_1, sym_2",
    [
        # Identical
        (sym(CELL, "P 21 21 21"), sym(CELL, "P 21 21 21")),
        # Same point group, different screw axes (pointless sorts these out)
        (sym(CELL, "P 2 2 2"), sym(CELL, "P 21 21 21")),
        (
            sym((92.4, 92.4, 128.2, 90, 90, 120), "P 61 2 2"),
            sym((92.4, 92.4, 128.3, 90, 90, 120), "P 63 2 2"),
        ),
        # Within the 1% tolerance
        (sym(CELL, "P 21 21 21"), sym((40.3, 45.3, 50.4, 90, 90, 90), "P 21 21 21")),
        # Axes in a different order (alternative indexing)
        (
            sym((40, 45, 50, 90, 90, 90), "P 2 2 2"),
            sym((45, 50, 40, 90, 90, 90), "P 2 2 2"),
        ),
        # Monoclinic cell given with an obtuse or acute beta
        (
            sym((50, 60, 70, 90, 100, 90), "P 1 21 1"),
            sym((50, 60, 70, 90, 80, 90), "P 1 21 1"),
        ),
    ],
)
def test_compatible_symmetry(sym_1, sym_2):
    assert symmetry_difference(sym_1, sym_2, 0.01, 1.0) is None


@pytest.mark.parametrize(
    "sym_1, sym_2, reason",
    [
        (sym(CELL, "P 21 21 21"), sym(CELL, "C 2 2 21"), "point group or lattice"),
        (sym(CELL, "P 2 2 2"), sym((40, 45, 50, 90, 90, 90), "P 4 2 2"), "point group"),
        (
            sym(CELL, "P 21 21 21"),
            sym((40, 45, 51, 90, 90, 90), "P 21 21 21"),
            "unit cells differ",
        ),
        (
            sym((50, 60, 70, 90, 100, 90), "P 1 21 1"),
            sym((50, 60, 70, 90, 102, 90), "P 1 21 1"),
            "unit cells differ",
        ),
        ((gemmi.UnitCell(), None), sym(CELL, "P 21 21 21"), "missing"),
    ],
)
def test_incompatible_symmetry(sym_1, sym_2, reason):
    problem = symmetry_difference(sym_1, sym_2, 0.01, 1.0)
    assert problem is not None and reason in problem


def test_cell_tolerance():
    sym_1 = sym(CELL, "P 21 21 21")
    sym_2 = sym((41.0, 45.0, 50.0, 90, 90, 90), "P 21 21 21")
    assert symmetry_difference(sym_1, sym_2, 0.01, 1.0) is not None
    assert symmetry_difference(sym_1, sym_2, 0.03, 1.0) is None


def test_check_mtz_wavelength_from_anomalous_dataset(tmp_path):
    # The HKL_base dataset has a placeholder wavelength, which must not be used
    mtz_file = write_data_mtz(tmp_path / "data.mtz", WAVELENGTH_ABOVE)
    assert check_mtz(mtz_file) == pytest.approx(WAVELENGTH_ABOVE)
    mtz_file = write_data_mtz(
        tmp_path / "data0.mtz", WAVELENGTH_ABOVE, base_wavelength=0
    )
    assert check_mtz(mtz_file) == pytest.approx(WAVELENGTH_ABOVE)


def test_check_mtz_missing_column(tmp_path):
    mtz_file = write_data_mtz(
        tmp_path / "data.mtz", WAVELENGTH_ABOVE, labels=("IMEAN", "SIGIMEAN")
    )
    with pytest.raises(SystemExit):
        check_mtz(mtz_file)


def test_check_mtz_missing_wavelength(tmp_path):
    mtz_file = write_data_mtz(tmp_path / "data.mtz", 0.0)
    with pytest.raises(SystemExit):
        check_mtz(mtz_file)


def test_mtz_and_model_symmetry(tmp_path):
    mtz_cell, mtz_sg = mtz_symmetry(
        write_data_mtz(tmp_path / "d.mtz", WAVELENGTH_ABOVE)
    )
    model_cell, model_sg = model_symmetry(write_model(tmp_path / "m.pdb"))
    assert mtz_sg.hm == model_sg.hm == "P 21 21 21"
    assert mtz_cell.approx(model_cell, 1e-3)


def test_find_matching_model(tmp_path):
    data_sym = {"above": sym(CELL, "P 21 21 21"), "below": sym(CELL, "P 21 21 21")}
    wrong_cell = write_model(tmp_path / "big.pdb", cell=(42, 47, 52, 90, 90, 90))
    no_cell = tmp_path / "no_cell.pdb"
    no_cell.write_text(
        "".join(
            line + "\n"
            for line in make_structure().make_pdb_string().splitlines()
            if not line.startswith("CRYST1")
        )
    )
    matching = write_model(tmp_path / "match.pdb")

    assert find_matching_model([wrong_cell, no_cell], data_sym, 0.01, 1.0) is None
    assert (
        find_matching_model([wrong_cell, no_cell, matching], data_sym, 0.01, 1.0)
        == matching
    )
    assert find_matching_model([wrong_cell], data_sym, 0.06, 1.0) == wrong_cell


def test_copy_model(tmp_path):
    out_dir = tmp_path / "out"
    out_dir.mkdir()
    pdb = copy_model(write_model(tmp_path / "model.pdb"), out_dir)
    assert pdb == out_dir / "model.pdb" and pdb.is_file()
    # Copying a model that is already in the output directory leaves it in place
    assert copy_model(pdb, out_dir) == pdb

    # mmCIF models are converted to PDB format (ANODE needs PDB)
    cif = tmp_path / "other.cif"
    make_structure().make_mmcif_document().write_file(str(cif))
    pdb = copy_model(cif, out_dir)
    assert pdb == out_dir / "other.pdb"
    assert gemmi.read_structure(str(pdb)).find_spacegroup().hm == "P 21 21 21"
