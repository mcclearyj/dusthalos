# dusthalos
Dusty, dusty halos! README still under construction. 

This repository takes galaxy catalogs from sources like COSMOS-Web or SuperBIT (mainly redMaGiC) and uses their data to calculate dust reddening corrections for the Milky Way, cross-correlations between dust extinction and galaxy source density, and dereddening for weak lensing. 

## Installation 
```bash
git clone https://github.com/mcclearyj/dusthalos
cd dusthalos
pip install -e .
```

Dependencies:
- numpy
- astropy
- healpy
- dustmaps 
- matplotlib
- scipy
- fitsio

Plots look best with TeX fonts; you can follow Unix installation instructions for TeX live here: https://www.tug.org/texlive/quickinstall.html
Once installed, add something like this to your path:
`export PATH='.':$PATH:'/path/to/your/texlive-bin/x86_64-linux'`

## Process Catalogs


`python runner_scripts/cat_prep_runner.py --config configs/prep_catalog_config.yaml`

`cat_prep_runner.py` is the first of two main runner scripts. This script takes the inputted galaxy catalogs and processes them to be used throughout the code and calculations. 


## Run Reddening Calculation
`python runner_scripts/dust_calc_runner.py --config configs/dust_calc_config.yaml`

`dust_calc_runner.py` is the second main runner script that calculates the dust reddening and extinction of the galaxy catalogs. 

## Configuration
The `dust_calc_config.yaml` configuration file points to catalog-level configuration files in `configs/catalog_configs`. 
These do useful things like specify the WCS coordinate system of the catalog and the lon/lat column names, the name of the redshift column, and the path to the catalog itself. 

The parameters of the extinction model to use is also specified in `dust_calc_config.yaml` [here]([url](https://github.com/mcclearyj/dusthalos/blob/4371d0fc41bbef86c4ce33c1e07c525801027854/configs/dust_calc_config.yaml#L38))

## Operation Pipeline
The flow is as follows:
1. Prepare catalog: process, reads, and filters galaxy catalogs  through `catalog.py` and `cat_prep_runner.py`
```bash
python runner_scripts/cat_prep_runner.py -config configs/prep_catalog_config.yaml
```
2. HEALPix coordinates: converts the galaxy catalog RA/Dec sky coordinates into HEALPix indices for dust mapping. Matches catalog areas based on HEALPix map. 
3. Optional: MW dust extinction correction, if not already de-reddened.
4. Outputs dereddened, masked catalog. 
5. Actual reddening calculation and correlation: 

```bash
python runner_scripts/dust_calc_runner.py -config configs/dust_calc_config.yaml
```

Note: Steps 3 and 4 happen internally within ``prep_cat_runner.py `` script. 

## Output
TBD

## Citation 
If you use this code in your research, please cite:

McCleary, J. E., Huff, E. M., Bartlett, J. G., & Hensley, B. S. (2026). A Detection of Circumgalactic Dust at Megaparsec Scales with Maximum Likelihood Estimation. The Astrophysical Journal, 1000(2), 313, IOP. https://doi.org/10.3847/1538-4357/ae4c3d
