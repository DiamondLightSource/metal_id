import logging

import gemmi
import numpy as np
import pytest

from metal_id.calc_map import (
    anode_maps,
    calc_double_diff_maps,
    find_peaks,
    hkl_keys,
    read_pha,
    refmac_anomalous_maps,
)
from metal_id.viewer import parse_peaks
from synthetic import (
    CA_SITE,
    S_SITE,
    WAVELENGTH_ABOVE,
    WAVELENGTH_BELOW,
    ZN_SITE,
    make_structure,
    write_model,
    write_pha,
    write_refmac_mtz,
)


def distance(structure, xyz, site):
    """Distance from xyz to the nearest symmetry copy of site"""
    return structure.cell.find_nearest_image(
        gemmi.Position(*site), gemmi.Position(*xyz), gemmi.Asu.Any
    ).dist()


def test_read_pha(tmp_path):
    # Large negative indices run into one another in the fixed width (3I4) columns
    pha = tmp_path / "test.pha"
    pha.write_text(
        "   1   2   3     10.000  0.5000   45.00   1.00\n"
        "-100-120-140      4.000  1.0000  -90.00   2.00\n"
        "\n"
    )
    data = read_pha(pha)
    np.testing.assert_allclose(
        data, [[1, 2, 3, 5.0, 45.0], [-100, -120, -140, 4.0, -90.0]]
    )


def test_hkl_keys_unique():
    hkl = np.stack(np.meshgrid(*[range(-3, 4)] * 3), axis=-1).reshape(-1, 3)
    keys = hkl_keys(np.column_stack([hkl, np.zeros(len(hkl))]))
    assert len(np.unique(keys)) == len(hkl)


@pytest.fixture
def model(tmp_path):
    return write_model(tmp_path / "model.pdb")


@pytest.fixture
def structure(model):
    return gemmi.read_structure(str(model))


def test_find_peaks(tmp_path, model, structure):
    # Anomalous map with Ca, Zn and S, in P212121 (4 symmetry copies of each site)
    mtz = write_refmac_mtz(tmp_path / "anom.mtz", model, WAVELENGTH_ABOVE, 2.0)
    grid, _ = refmac_anomalous_maps(mtz, mtz, structure)
    peaks = find_peaks(grid, structure, 5.0, 10)

    # One peak per site (no symmetry copies), largest first: f'' Ca > Zn > S
    assert len(peaks) == 3
    assert [p.rmsd for p in peaks] == sorted((p.rmsd for p in peaks), reverse=True)
    for peak, site in zip(peaks, [CA_SITE, ZN_SITE, S_SITE]):
        # Moved to the copy next to the model, not just a symmetry equivalent
        assert np.linalg.norm(np.subtract(peak.xyz, site)) < 0.3
        assert peak.density == pytest.approx(
            peak.rmsd * float(np.std(grid.array)), rel=0.01
        )


def test_find_peaks_max_peaks(tmp_path, model, structure, caplog):
    mtz = write_refmac_mtz(tmp_path / "anom.mtz", model, WAVELENGTH_ABOVE, 2.0)
    grid, _ = refmac_anomalous_maps(mtz, mtz, structure)
    with caplog.at_level(logging.WARNING):
        peaks = find_peaks(grid, structure, 5.0, 1)
    assert len(peaks) == 1
    assert distance(structure, peaks[0].xyz, CA_SITE) < 0.3
    assert "maximum of 1" in caplog.text


def test_refmac_maps_use_common_reflections(tmp_path, model, structure):
    # The above data extend to higher resolution than the below data. With the same
    # anomalous scattering in both, the maps should be identical once cut to the
    # common reflections.
    above = write_refmac_mtz(tmp_path / "above.mtz", model, WAVELENGTH_ABOVE, 1.8)
    mtz = gemmi.read_mtz_file(str(above))
    mtz.set_data(np.array(mtz)[mtz.make_d_array() >= 2.4])
    below = tmp_path / "below.mtz"
    mtz.write_to_file(str(below))
    map_above, map_below = refmac_anomalous_maps(above, below, structure)
    assert map_above.array.shape == map_below.array.shape
    np.testing.assert_allclose(map_above.array, map_below.array, atol=1e-5)


def test_refmac_maps_drop_missing_reflections(tmp_path, model, structure):
    above = write_refmac_mtz(tmp_path / "above.mtz", model, WAVELENGTH_ABOVE, 2.0)
    below = write_refmac_mtz(tmp_path / "below.mtz", model, WAVELENGTH_ABOVE, 2.0)
    # Remove FAN from some reflections in the below data
    mtz = gemmi.read_mtz_file(str(below))
    data = np.array(mtz, copy=True)
    data[::7, mtz.column_labels().index("FAN")] = np.nan
    mtz.set_data(data)
    mtz.write_to_file(str(below))

    map_above, map_below = refmac_anomalous_maps(above, below, structure)
    np.testing.assert_allclose(map_above.array, map_below.array, atol=1e-5)


def test_double_difference_finds_ca(tmp_path, model, structure):
    # Only Ca changes much across its edge, so it is the only peak in the difference
    above = write_refmac_mtz(tmp_path / "above.mtz", model, WAVELENGTH_ABOVE, 2.0)
    below = write_refmac_mtz(tmp_path / "below.mtz", model, WAVELENGTH_BELOW, 2.0)
    map_above, map_below = refmac_anomalous_maps(above, below, structure)
    out = tmp_path / "out"
    out.mkdir()
    files = calc_double_diff_maps(model, map_above, map_below, out, 5.0, 10)

    for name in [
        "final.pdb",
        "diff.map",
        "anom_above.map",
        "anom_below.map",
        "found_peaks.dat",
    ]:
        assert (out / name).is_file()
    peaks = parse_peaks(files["peaks"])
    assert len(peaks) == 1
    assert distance(structure, (peaks[0].x, peaks[0].y, peaks[0].z), CA_SITE) < 0.3


def test_double_difference_no_peaks(tmp_path, model, structure):
    # Identical data either side of the edge: no peaks, and no peak file
    mtz = write_refmac_mtz(tmp_path / "same.mtz", model, WAVELENGTH_ABOVE, 2.0)
    map_above, map_below = refmac_anomalous_maps(mtz, mtz, structure)
    map_below.array[:] += np.random.default_rng(0).normal(
        0, 1e-4, map_below.array.shape
    )
    files = calc_double_diff_maps(model, map_above, map_below, tmp_path, 5.0, 10)
    assert "peaks" not in files
    assert not (tmp_path / "found_peaks.dat").exists()


def test_anode_maps(tmp_path, model, structure):
    above = write_pha(tmp_path / "above.pha", model, WAVELENGTH_ABOVE, 2.0)
    below = write_pha(tmp_path / "below.pha", model, WAVELENGTH_BELOW, 2.0)
    map_above, map_below = anode_maps(above, below, make_structure())
    assert map_above.array.shape == map_below.array.shape
    diff = map_above.clone()
    diff.array[:] -= map_below.array
    peaks = find_peaks(diff, structure, 5.0, 10)
    assert len(peaks) == 1
    assert distance(structure, peaks[0].xyz, CA_SITE) < 0.3
