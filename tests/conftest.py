import os
import shutil

import pytest

import metal_id.anomalous
import metal_id.prepare_data
from synthetic import FakeCCP4

# CCP4 programs needed by the tests marked "ccp4"
CCP4_PROGRAMS = [
    "pointless",
    "ctruncate",
    "refmac5",
    "mtz2sca",
    "shelxc",
    "anode",
    "dimple",
]


def pytest_collection_modifyitems(config, items):
    """
    Skip the tests that need CCP4 if it isn't loaded, unless METAL_ID_REQUIRE_CCP4 is
    set (e.g. on Jenkins), in which case stop with an error instead.
    """
    missing = [program for program in CCP4_PROGRAMS if shutil.which(program) is None]
    if not missing:
        return
    if os.environ.get("METAL_ID_REQUIRE_CCP4"):
        raise pytest.UsageError(
            f"METAL_ID_REQUIRE_CCP4 is set, but CCP4 isn't available ({', '.join(missing)} not found)"
        )
    skip = pytest.mark.skip(
        reason=f"needs CCP4 ({', '.join(missing)} not found) - load CCP4 to run"
    )
    for item in items:
        if "ccp4" in item.keywords:
            item.add_marker(skip)


@pytest.fixture
def fake_ccp4(monkeypatch) -> FakeCCP4:
    """Replace the CCP4 programs with stand-ins calculated from the model"""
    fake = FakeCCP4()
    monkeypatch.setattr(metal_id.prepare_data, "run_ccp4", fake)
    monkeypatch.setattr(metal_id.anomalous, "run_ccp4", fake)
    return fake
