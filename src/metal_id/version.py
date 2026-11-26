import subprocess
from pathlib import Path

# In a release made with git archive, this is replaced by the output of git describe
# (see .gitattributes), e.g. "v0.2.0"
_ARCHIVE_VERSION = "$Format:%(describe)$"


def get_version() -> str:
    """
    Return the metal_id version: the release tag for an installed release, or the
    output of git describe for a git checkout (e.g. "v0.2.0-3-gabc1234-dirty").
    """
    if not _ARCHIVE_VERSION.startswith("$Format"):
        # Blank if there is no annotated tag before the release
        return _ARCHIVE_VERSION or "unknown"
    try:
        result = subprocess.run(
            ["git", "describe", "--tags", "--always", "--dirty"],
            cwd=Path(__file__).parent,
            capture_output=True,
            text=True,
            check=True,
        )
    except (OSError, subprocess.CalledProcessError):
        return "unknown"
    return result.stdout.strip()
