"""Tests of the whole metal_id workflow, with stand-ins for the CCP4 programs"""

import os
import sys
import subprocess
from pathlib import Path

import gemmi
import numpy as np
import pytest

import metal_id.metal_id
from metal_id.metal_id import main
from metal_id.viewer import parse_peaks
from synthetic import (
    CA_SITE,
    WAVELENGTH_ABOVE,
    WAVELENGTH_BELOW,
    write_data_mtz,
    write_model,
)

REPO_ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture
def data(tmp_path):
    """Above and below data, at different resolutions, and a matching model"""
    return {
        "above": write_data_mtz(tmp_path / "above.mtz", WAVELENGTH_ABOVE, d_min=2.0),
        "below": write_data_mtz(tmp_path / "below.mtz", WAVELENGTH_BELOW, d_min=2.2),
        "model": write_model(tmp_path / "model.pdb"),
    }


class FakeDimple:
    """Stand-in for run_dimple that writes a model in the cell of the data"""

    def __init__(self, succeed=True):
        self.succeed = succeed
        self.calls = []

    def __call__(self, mtz, pdb, dimple_dir):
        self.calls.append((mtz, pdb, dimple_dir))
        dimple_dir.mkdir(parents=True)
        if self.succeed:
            write_model(dimple_dir / "final.pdb")
        return subprocess.CompletedProcess(
            [], 0 if self.succeed else 1, "fake dimple\n"
        )


@pytest.fixture
def fake_dimple(monkeypatch):
    fake = FakeDimple()
    monkeypatch.setattr(metal_id.metal_id, "run_dimple", fake)
    return fake


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


def assert_single_ca_peak(out):
    peaks = parse_peaks(out / "found_peaks.dat")
    assert len(peaks) == 1
    assert (
        np.linalg.norm(np.subtract((peaks[0].x, peaks[0].y, peaks[0].z), CA_SITE)) < 0.3
    )


def test_refmac_workflow(tmp_path, data, fake_ccp4, fake_dimple):
    out = tmp_path / "out"
    run(data, out)

    assert_single_ca_peak(out)
    for name in [
        "metal_id.log",
        "final.pdb",
        "diff.map",
        "anom_above.map",
        "anom_below.map",
        "refmac_above.mtz",
        "refmac_below.mtz",
        "coot_view.py",
        "metal_id.html",
    ]:
        assert (out / name).is_file(), name
    assert not fake_dimple.calls
    assert not fake_ccp4.calls_to("anode")

    # Each dataset is reindexed against, and refined with, the same model
    models = {
        call[1][call[1].index("XYZIN") + 1] for call in fake_ccp4.calls_to("pointless")
    }
    models |= {
        call[1][call[1].index("XYZIN") + 1] for call in fake_ccp4.calls_to("refmac5")
    }
    assert models == {str(out / "model.pdb")}

    # refmac gets each dataset's own wavelength, not the HKL_base placeholder
    wavelengths = [
        next(k for k in keywords if k.startswith("ANOMALOUS WAVELENGTH"))
        for _, _, keywords in fake_ccp4.calls_to("refmac5")
    ]
    assert wavelengths == [
        f"ANOMALOUS WAVELENGTH {WAVELENGTH_ABOVE:.5f}",
        f"ANOMALOUS WAVELENGTH {WAVELENGTH_BELOW:.5f}",
    ]

    log = (out / "metal_id.log").read_text()
    assert log.startswith("metal_id version: ")
    assert "Using the" in log and "to 2.20 Å" in log


def test_anode_workflow(tmp_path, data, fake_ccp4, fake_dimple):
    out = tmp_path / "out"
    run(data, out, "--method", "anode")

    assert_single_ca_peak(out)
    assert len(fake_ccp4.calls_to("anode")) == 2
    # refmac still runs, for the R factor check and the 2Fo-Fc map
    assert len(fake_ccp4.calls_to("refmac5")) == 2


def test_dimple_when_model_does_not_match(tmp_path, data, fake_ccp4, fake_dimple):
    out = tmp_path / "out"
    big_model = write_model(tmp_path / "big.pdb", cell=(41.2, 46.4, 51.5, 90, 90, 90))
    run(data, out, model=big_model)

    ((mtz, pdb, dimple_dir),) = fake_dimple.calls
    assert mtz == data["above"].resolve()
    assert pdb == [str(big_model.resolve())]
    # The model from dimple is used for both datasets
    final = gemmi.read_structure(str(out / "final.pdb"))
    assert final.cell.a == pytest.approx(40.0)
    assert_single_ca_peak(out)


def test_cell_tolerance_option(tmp_path, data, fake_ccp4, fake_dimple):
    big_model = write_model(tmp_path / "big.pdb", cell=(41.2, 46.4, 51.5, 90, 90, 90))
    run(data, tmp_path / "out", "--cell-tolerance", "5", model=big_model)
    assert not fake_dimple.calls


def test_dimple_failure(tmp_path, data, fake_ccp4, monkeypatch):
    monkeypatch.setattr(metal_id.metal_id, "run_dimple", FakeDimple(succeed=False))
    big_model = write_model(tmp_path / "big.pdb", cell=(41.2, 46.4, 51.5, 90, 90, 90))
    with pytest.raises(SystemExit):
        run(data, tmp_path / "out", model=big_model)
    assert not fake_ccp4.calls


def test_incompatible_datasets(tmp_path, data, fake_ccp4, fake_dimple):
    data["below"] = write_data_mtz(
        tmp_path / "below_c2.mtz", WAVELENGTH_BELOW, spacegroup="C 2 2 21"
    )
    with pytest.raises(SystemExit):
        run(data, tmp_path / "out")
    assert not fake_ccp4.calls and not fake_dimple.calls


def test_r_factor_too_high(tmp_path, data, fake_ccp4, fake_dimple):
    with pytest.raises(SystemExit):
        run(data, tmp_path / "out", "--max-r-factor", "0.1")
    assert len(fake_ccp4.calls_to("refmac5")) == 1


def test_existing_output_directory(tmp_path, data, fake_ccp4, fake_dimple):
    (tmp_path / "out").mkdir()
    run(data, tmp_path / "out")
    assert (tmp_path / "out_1" / "found_peaks.dat").is_file()


def test_bin_wrapper_ignores_pythonpath(tmp_path):
    # Regression test for #19: dlstbx.python may reset PYTHONPATH, and a "metal_id"
    # output directory in the current directory must not be imported instead
    # The stand-in dlstbx.python replaces PYTHONPATH with a path to the dependencies,
    # and skips site processing (-S) so that an installed metal_id (e.g. an editable
    # install) isn't importable, as in the real dlstbx environment
    site_packages = Path(gemmi.__file__).resolve().parents[1]
    fake_bin = tmp_path / "bin"
    fake_bin.mkdir()
    dlstbx_python = fake_bin / "dlstbx.python"
    dlstbx_python.write_text(
        f'#!/bin/bash\nexport PYTHONPATH="{site_packages}"\nexec {sys.executable} -S "$@"\n'
    )
    dlstbx_python.chmod(0o755)
    (tmp_path / "metal_id").mkdir()

    result = subprocess.run(
        [str(REPO_ROOT / "bin" / "metal_id"), "--version"],
        cwd=tmp_path,
        env={**os.environ, "PATH": f"{fake_bin}:{os.environ['PATH']}"},
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stderr
    assert result.stdout.startswith("metal_id ")
