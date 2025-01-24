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
Takes two mtz files and a pdb file as input. The mtz files should correspond to data collected above and below an element's absorption edge and the pdb file should contain the protein strucuture. The pipeline then combines and scales the two data files to put them on a matching coarse scale, runs dimple on both datasets with the `--anode` option set to produce anomalous Fourier maps, then subtracts the 'below' map from the 'above' map to generate an anomalous Fourier 'double' difference that gives the location(s) of the element in question.

Prior to scaling, `pointless` is run on the mtz files to ensure that they are of compatible symmetry and to index them into matching space groups if they are not already. The second run of `dimple` uses the output `final.pdb` file from the first run of dimple, to ensure that the structure is defined in a similar way if any molecular replacement methods are needed in the model generation. A check is then made to ensure that the output `final.pdb` files from both runs of dimple are sufficently similar, before calculating the 'double difference' map `diff.map`. 