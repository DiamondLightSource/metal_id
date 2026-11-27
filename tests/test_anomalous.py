import pytest

from metal_id.anomalous import run_anode, run_ctruncate, run_refmac_anomalous
from metal_id.metal_id_helpers import run_ccp4
from synthetic import WAVELENGTH_ABOVE, write_data_mtz, write_model


def test_run_ccp4_passes_keywords_and_writes_log(tmp_path):
    log = tmp_path / "cat.log"
    output = run_ccp4("cat", [], log, keywords=["LINE ONE", "END"])
    assert output == "LINE ONE\nEND\n"
    assert log.read_text() == output


def test_run_ccp4_uses_cwd(tmp_path):
    output = run_ccp4("pwd", [], tmp_path / "pwd.log", cwd=tmp_path)
    assert output.strip() == str(tmp_path)


def test_run_ccp4_exits_on_failure(tmp_path):
    with pytest.raises(SystemExit):
        run_ccp4("false", [], tmp_path / "false.log")
    assert (tmp_path / "false.log").exists()


def test_run_ccp4_exits_if_program_missing(tmp_path):
    with pytest.raises(SystemExit):
        run_ccp4("not_a_ccp4_program", [], tmp_path / "missing.log")


def test_run_refmac_anomalous(tmp_path, fake_ccp4):
    model = write_model(tmp_path / "model.pdb")
    data = write_data_mtz(tmp_path / "data.mtz", WAVELENGTH_ABOVE)
    truncated = run_ctruncate(data, tmp_path)
    mtz, r_factor = run_refmac_anomalous(truncated, model, WAVELENGTH_ABOVE, tmp_path)

    assert mtz == tmp_path / "refmac.mtz" and mtz.is_file()
    assert r_factor == pytest.approx(0.2)
    ((_, _, keywords),) = fake_ccp4.calls_to("refmac5")
    assert f"ANOMALOUS WAVELENGTH {WAVELENGTH_ABOVE:.5f}" in keywords
    assert "NCYC 0" in keywords
    assert "MAKE NEWLIGAND CONTINUE" in keywords


def test_run_refmac_anomalous_exits_without_r_factor(tmp_path, fake_ccp4, monkeypatch):
    model = write_model(tmp_path / "model.pdb")
    data = write_data_mtz(tmp_path / "data.mtz", WAVELENGTH_ABOVE)
    monkeypatch.setattr(
        type(fake_ccp4), "refmac", lambda self, *args: "no R factor here\n"
    )
    with pytest.raises(SystemExit):
        run_refmac_anomalous(data, model, WAVELENGTH_ABOVE, tmp_path)


def test_run_anode(tmp_path, fake_ccp4):
    model = write_model(tmp_path / "model.pdb")
    data = write_data_mtz(tmp_path / "data.mtz", WAVELENGTH_ABOVE)
    pha = run_anode(data, model, tmp_path)

    assert pha == tmp_path / "anode" / "anode.pha" and pha.is_file()
    assert [call[0] for call in fake_ccp4.calls] == ["mtz2sca", "shelxc", "anode"]
    ((_, _, keywords),) = fake_ccp4.calls_to("shelxc")
    assert keywords == [
        "SAD anode.sca",
        "CELL 40.0000 45.0000 50.0000 90.0000 90.0000 90.0000",
        "SPAG P 21 21 21",
    ]
