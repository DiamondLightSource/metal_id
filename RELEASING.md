# Releasing metal_id

Each release is a git tag, installed into its own directory in `/dls_sw/apps/metal_id`
with its own module file, so that earlier versions stay available.

## Version numbers

Versions are `MAJOR.MINOR.PATCH`, tagged as `vMAJOR.MINOR.PATCH` (e.g. `v0.2.0`):

- **PATCH** for bug fixes that don't change the results.
- **MINOR** for new features, or changes that alter the results (e.g. a change to the
  maps or peak finding).
- **MAJOR** for changes to the command line that break existing usage.

## Making a release

1. Make sure everything for the release is merged into `main`, and update your
   checkout:

   ```
   git checkout main
   git pull
   ```

2. Tag the release with an annotated tag (`-a`), with a short summary of the changes
   as the message, and push the tag:

   ```
   git tag -a v0.2.0 -m "Anomalous maps from refmac against a single model"
   git push origin v0.2.0
   ```

   Optionally, create a GitHub release from the tag with fuller release notes.

3. Install the release with `git archive`, which writes exactly the tagged files,
   without `.git`, and fills in the version (see below):

   ```
   mkdir /dls_sw/apps/metal_id/0.2.0
   git archive v0.2.0 | tar -x -C /dls_sw/apps/metal_id/0.2.0
   ```

4. Add a module file for the new version, by copying the previous version's module
   file and changing the version and path to the new directory. Test it:

   ```
   module load metal_id/0.2.0
   metal_id --version
   ```

5. Make the new version the default, e.g. by updating the `.version` file (or
   `default` symlink) in the metal_id module file directory.

To roll back, set the default module version back to the previous release.

## Checking the version

`metal_id --version` prints the version, and it is written at the top of every
`metal_id.log`. An installed release reports its tag (e.g. `v0.2.0`). A git checkout
reports the output of `git describe` (e.g. `v0.2.0-3-gabc1234-dirty`: 3 commits after
v0.2.0, with uncommitted changes).

The version comes from the git tag, in one of three ways:

- In a git checkout, from `git describe`.
- In a release made with `git archive`, from `.git_archival.txt`, which `git archive`
  fills in with the tag (see `.gitattributes`).
- When installed with pip, from the package version, which `setuptools-scm` takes from
  the tag when the package is built (e.g. `0.2.0`, without the `v`).

So don't install by copying a clone with `.git` deleted: the version would then be
reported as `unknown`.

## Installing with pip instead

metal_id can also be installed with pip into its own environment for each release,
e.g.:

```
python -m venv /dls_sw/apps/metal_id/0.2.0
/dls_sw/apps/metal_id/0.2.0/bin/pip install git+https://github.com/DiamondLightSource/metal_id@v0.2.0
```

The module file then puts that environment's `bin` directory on `PATH` (and loads
CCP4). This doesn't need `dlstbx.python`, but installs its own copies of gemmi, numpy
and molviewspec.
