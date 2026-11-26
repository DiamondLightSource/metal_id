# metal_id
Pipeline for analysing metal (or heavy element) identification experiments in macromolecular crystallography experiments.

## Usage
This pipeline is dependent on the dials python environment (`cctbx.python`) and ccp4 packages. In the dls_sw installation, these are loaded in as a dependenc and thus the pipeline can be used as follows:
```
module load metal_id

metal_id mtz_above mtz_below pdb output_dir
```

Otherwise, load dials and ccp4, then from this repository call: 

```
bin/metal_id mtz_above mtz_below pdb output_dir
```

If no output directory name is specified, the pipeline will output to a directory named `metal_id`

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
