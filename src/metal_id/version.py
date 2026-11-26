import subprocess
from importlib import metadata
from pathlib import Path

# Root of the repository, or of a release made with git archive
# (src/metal_id/version.py)
_REPO_ROOT = Path(__file__).resolve().parents[2]


def _archive_version() -> str | None:
    """
    Return the version from .git_archival.txt, which git archive fills in with the
    output of git describe (see .gitattributes), or None if there isn't one.
    """
    archival_file = _REPO_ROOT / ".git_archival.txt"
    if not archival_file.is_file():
        return None
    for line in archival_file.read_text().splitlines():
        key, _, value = line.partition(":")
        value = value.strip()
        if key == "describe-name" and value and not value.startswith("$Format"):
            return value
    return None


def _git_version() -> str | None:
    """Return the output of git describe for a git checkout, or None"""
    if not (_REPO_ROOT / ".git").exists():
        return None
    try:
        result = subprocess.run(
            ["git", "describe", "--tags", "--always", "--dirty"],
            cwd=_REPO_ROOT,
            capture_output=True,
            text=True,
            check=True,
        )
    except (OSError, subprocess.CalledProcessError):
        return None
    return result.stdout.strip()


def get_version() -> str:
    """
    Return the metal_id version: the output of git describe for a git checkout (e.g.
    "v0.2.0-3-gabc1234-dirty"), the release tag for a release made with git archive
    (e.g. "v0.2.0"), or the package version for a pip install (e.g. "0.2.0").
    """
    version = _git_version() or _archive_version()
    if version:
        return version
    try:
        return metadata.version("metal_id")
    except metadata.PackageNotFoundError:
        return "unknown"
