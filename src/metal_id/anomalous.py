import logging
import re
import shutil
import sys
from pathlib import Path

import gemmi

from metal_id_helpers import run_ccp4


def run_ctruncate(mtz_file: Path, work_dir: Path) -> Path:
    """Calculate amplitudes, including F(+)/F(-), from the intensities with ctruncate"""
    hklout = work_dir / "truncated.mtz"
    run_ccp4(
        "ctruncate",
        [
            "-hklin",
            mtz_file,
            "-hklout",
            hklout,
            "-colin",
            "/*/*/[IMEAN,SIGIMEAN]",
            "-colano",
            "/*/*/[I(+),SIGI(+),I(-),SIGI(-)]",
        ],
        work_dir / "ctruncate.log",
    )
    return hklout


def run_refmac_anomalous(
    mtz_file: Path, model: Path, wavelength: float, work_dir: Path
) -> tuple[Path, float]:
    """
    Run refmac with the anomalous amplitudes, so that the data are scaled against the
    model and anomalous map coefficients (FAN/PHAN) are calculated. No refinement is
    done, so both datasets are phased with exactly the same model. f' and f'' are set
    from the wavelength. Returns the output mtz file and the overall R factor.
    """
    hklout = work_dir / "refmac.mtz"
    log_file = work_dir / "refmac.log"
    refmac_log = run_ccp4(
        "refmac5",
        [
            "HKLIN",
            mtz_file,
            "HKLOUT",
            hklout,
            "XYZIN",
            model,
            "XYZOUT",
            work_dir / "refmac.pdb",
        ],
        log_file,
        keywords=[
            "LABIN FP=F SIGFP=SIGF F+=F(+) SIGF+=SIGF(+) F-=F(-) SIGF-=SIGF(-)",
            f"ANOMALOUS WAVELENGTH {wavelength:.5f}",
            "NCYC 0",
            # refmac sets up restraints even with no refinement, and by default stops
            # if a residue isn't in the monomer library. The restraints aren't used for
            # the scaling or maps, so carry on with generated ones (as dimple does)
            "MAKE NEWLIGAND CONTINUE",
            "END",
        ],
    )
    r_factors = re.findall(r"Overall R factor\s*=\s*([\d.]+)", refmac_log)
    if not hklout.is_file() or not r_factors:
        logging.error(f"refmac did not complete, see {log_file}")
        sys.exit(1)
    return hklout, float(r_factors[-1])


def run_anode(mtz_file: Path, model: Path, work_dir: Path) -> Path:
    """
    Calculate an anomalous map with ANODE, as dimple does: mtz2sca, SHELXC, then ANODE
    using the model for the phases. Returns the path of the .pha file.
    """
    anode_dir = work_dir / "anode"
    anode_dir.mkdir()
    shutil.copyfile(model, anode_dir / "anode.pdb")
    mtz = gemmi.read_mtz_file(str(mtz_file), with_data=False)
    cell = " ".join(f"{p:.4f}" for p in mtz.cell.parameters)

    run_ccp4(
        "mtz2sca", [mtz_file, "anode.sca"], anode_dir / "mtz2sca.log", cwd=anode_dir
    )
    run_ccp4(
        "shelxc",
        ["anode"],
        anode_dir / "shelxc.log",
        keywords=["SAD anode.sca", f"CELL {cell}", f"SPAG {mtz.spacegroup.hm}"],
        cwd=anode_dir,
    )
    run_ccp4("anode", ["anode"], anode_dir / "anode.log", cwd=anode_dir)

    pha_file = anode_dir / "anode.pha"
    if not pha_file.is_file():
        logging.error(
            f"ANODE did not produce {pha_file}, see {anode_dir / 'anode.log'}"
        )
        sys.exit(1)
    return pha_file
