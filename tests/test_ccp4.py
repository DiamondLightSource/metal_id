"""
Tests of the whole workflow with the real CCP4 programs, on synthetic data with a Ca
site whose anomalous signal changes across the Ca K edge.

These are skipped unless CCP4 is loaded. To run only these:

    module load ccp4
    pytest -m ccp4
"""

import re

import gemmi
import pytest

from metal_id.metal_id import main
from metal_id.viewer import parse_peaks
from synthetic import (
    CA_SITE,
    WAVELENGTH_ABOVE,
    WAVELENGTH_BELOW,
    write_anomalous_data_mtz,
    write_model,
)

pytestmark = pytest.mark.ccp4


@pytest.fixture(scope="module")
def data(tmp_path_factory):
    """Above and below data, at different resolutions, and a matching model"""
    data_dir = tmp_path_factory.mktemp("data")
    return {
        "above": write_anomalous_data_mtz(
            data_dir / "above.mtz", WAVELENGTH_ABOVE, d_min=2.0, seed=1
        ),
        "below": write_anomalous_data_mtz(
            data_dir / "below.mtz", WAVELENGTH_BELOW, d_min=2.2, seed=2
        ),
        "model": write_model(data_dir / "model.pdb"),
    }


def run(data, out, *options, model=None):
    main(
        [
            str(data["above"]),
            str(data["below"]),
            str(model or data["model"]),
            "-o",
            str(out),
            *options,
        ]
    )


def assert_ca_is_top_peak(out):
    """The largest peak is at the Ca site, and much larger than any other"""
    structure = gemmi.read_structure(str(out / "final.pdb"))
    peaks = parse_peaks(out / "found_peaks.dat")
    top = gemmi.Position(peaks[0].x, peaks[0].y, peaks[0].z)
    ca = gemmi.Position(*CA_SITE)
    assert structure.cell.find_nearest_image(ca, top, gemmi.Asu.Any).dist() < 0.5
    assert all(peak.height < 0.5 * peaks[0].height for peak in peaks[1:])


def refmac_f_double_prime(log_file, element):
    """f'' that refmac used for an element, from its log"""
    match = re.search(
        rf"^\s*{element}\s+\d+\s+atoms, f'=.*f''=\s*([-\d.]+)",
        log_file.read_text(),
        re.M,
    )
    return float(match.group(1))


def test_refmac(tmp_path, data):
    out = tmp_path / "out"
    run(data, out)

    assert_ca_is_top_peak(out)
    # refmac used f'' for each dataset's own wavelength
    for name, wavelength in [("above", WAVELENGTH_ABOVE), ("below", WAVELENGTH_BELOW)]:
        expected = gemmi.cromer_liberman(z=20, energy=12398.42 / wavelength)[1]
        used = refmac_f_double_prime(out / name / "refmac.log", "CA")
        assert used == pytest.approx(expected, abs=0.05)


def test_anode(tmp_path, data):
    out = tmp_path / "out"
    run(data, out, "--method", "anode")

    assert_ca_is_top_peak(out)
    assert (out / "above" / "anode" / "anode.pha").is_file()


def test_dimple_when_model_does_not_match(tmp_path, data):
    # The model's cell is 2% larger than the data's, so dimple is run to get a model
    out = tmp_path / "out"
    model = write_model(tmp_path / "big.pdb", cell=(40.8, 45.9, 51.0, 90, 90, 90))
    run(data, out, model=model)

    assert (out / "dimple" / "final.pdb").is_file()
    assert_ca_is_top_peak(out)
