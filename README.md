# metal_id
Pipeline for analysing metal (or heavy element) identification experiments in macromolecular crystallography experiments.

## Usage
metal_id needs Python 3.11 or later with `gemmi`, `numpy` and `molviewspec`, and the CCP4 programs (`pointless`, `ctruncate`, `refmac5`, `dimple`, and for `--method anode` also `mtz2sca`, `shelxc` and `anode`). At Diamond, these are loaded with the module:
```
module load metal_id

metal_id mtz_above mtz_below pdb [additional pdbs] [-o output_dir] [options]
```

Run `metal_id --help` for the options and `metal_id --version` for the version. If no output directory is given, the pipeline outputs to a directory named `metal_id`.

### Installing
metal_id is a Python package and can be installed with pip, which provides the `metal_id` command (CCP4 still needs to be loaded separately):
```
pip install git+https://github.com/DiamondLightSource/metal_id@v0.2.0
```

Or, from a copy of this repository, it can be run without installing, using the dials/dlstbx Python environment:
```
bin/metal_id mtz_above mtz_below pdb
```
or with any Python environment that has the dependencies:
```
PYTHONPATH=src python -m metal_id mtz_above mtz_below pdb
```

For development, install it in editable mode with the test dependencies, and run the tests with pytest:
```
pip install -e ".[test]"
pytest
```
Most of the tests don't need CCP4: the CCP4 programs are replaced by stand-ins that calculate their outputs from a small synthetic structure (see `tests/synthetic.py`). The tests marked `ccp4` (in `tests/test_ccp4.py`) run the whole workflow with the real CCP4 programs, on synthetic anomalous data with a Ca site. They are skipped unless CCP4 is loaded:
```
module load ccp4
pytest              # all tests
pytest -m ccp4      # only the tests with the real CCP4 programs
```

## Pipeline overview
Takes two mtz files and one or more models (pdb files and/or 4 character PDB codes) as input. The mtz files should correspond to data collected above and below an element's absorption edge and must contain `IMEAN`, `I(+)` and `I(-)` columns and the wavelength. The pipeline calculates an anomalous map for each dataset and subtracts the 'below' map from the 'above' map to give an anomalous 'double' difference map, `diff.map`, whose peaks give the location(s) of the element in question.

1. **Check the data.** The two datasets must have the same point group and lattice and similar unit cells (Niggli reduced cells within `--cell-tolerance` percent, default 1%), otherwise metal_id stops.
2. **Find a model that matches the data.** PDB codes are downloaded, and the first model with the same point group and a similar unit cell to both datasets is used. If no model matches, dimple is run on the 'above' data (with molecular replacement if needed) and its refined model is used for both datasets. Note that after molecular replacement the model may be in a different origin/frame from the input model.
3. **Reindex.** `pointless` puts each dataset into the space group and indexing of the model.
4. **Scale against the model.** `ctruncate` calculates F(+)/F(-), then `refmac5` (0 cycles, so both datasets are phased with exactly the same model) scales each dataset against the model and calculates anomalous map coefficients (`FAN`/`PHAN`), with f′ and f″ set from the wavelength. metal_id stops if the R factor against either dataset is above `--max-r-factor` (default 0.35).
5. **Anomalous maps** (`--method`):
   - `refmac` (default): maps from the refmac `FAN`/`PHAN` coefficients, using only the reflections common to both datasets, to the resolution of the lower resolution dataset. As the data are scaled against the model, the maps are on an approximately absolute scale (e/Å^3), so the signal from atoms whose f″ hardly changes across the edge (e.g. S) cancels in the double difference.
   - `anode`: maps from ANODE (`mtz2sca`, `SHELXC`, `ANODE`, as run by dimple) phased with the model. ANODE normalises the data, so these maps are in units of σ.
6. **Double difference and peaks.** The 'below' map is subtracted from the 'above' map and peaks above `--peak-threshold` (in rmsd, default 5) are written to `found_peaks.dat`, along with a coot script and a Mol* viewer page (`metal_id.html`).

Main outputs: `diff.map`, `anom_above.map`, `anom_below.map`, `found_peaks.dat`, `final.pdb` (the model used), `refmac_above.mtz` and `refmac_below.mtz` (including `FWT`/`PHWT` and `FAN`/`PHAN`), `metal_id.html` and `coot_view.py`. Intermediate files and logs are in the `above` and `below` subdirectories (and `dimple`, if it was run).
