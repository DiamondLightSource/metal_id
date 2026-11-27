import subprocess
from importlib import metadata

import pytest

from metal_id import version


@pytest.fixture
def repo_root(tmp_path, monkeypatch):
    """Pretend that metal_id is installed in tmp_path (not a git checkout)"""
    monkeypatch.setattr(version, "_REPO_ROOT", tmp_path)
    return tmp_path


def test_archive_version(repo_root):
    (repo_root / ".git_archival.txt").write_text(
        "node: 3417050abc\nnode-date: 2026-11-26T15:00:00+00:00\ndescribe-name: v0.2.1\n"
    )
    assert version.get_version() == "v0.2.1"


def test_unfilled_archive_falls_back_to_package_version(repo_root, monkeypatch):
    # .git_archival.txt as committed, before git archive fills it in
    (repo_root / ".git_archival.txt").write_text(
        "describe-name: $Format:%(describe:tags=true,match=*[0-9]*)$\n"
    )
    monkeypatch.setattr(metadata, "version", lambda name: "0.2.1")
    assert version.get_version() == "0.2.1"


def test_unknown_version(repo_root, monkeypatch):
    def not_installed(name):
        raise metadata.PackageNotFoundError(name)

    monkeypatch.setattr(metadata, "version", not_installed)
    assert version.get_version() == "unknown"


def test_git_checkout_version(repo_root):
    def git(*args):
        subprocess.run(
            ["git", "-c", "user.name=t", "-c", "user.email=t@t", *args],
            cwd=repo_root,
            check=True,
            capture_output=True,
        )

    git("init", "-q")
    (repo_root / "file").write_text("x")
    git("add", "file")
    git("commit", "-qm", "test")
    git("tag", "-a", "v1.2.3", "-m", "test release")
    assert version.get_version() == "v1.2.3"
    (repo_root / "file").write_text("y")
    assert version.get_version() == "v1.2.3-dirty"
