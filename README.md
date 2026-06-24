# dusthalos
Dusty, dusty halos!

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
2. HEALPix coordinates: converts the galaxy catalog RA/Dec sky coordinates into HEALPix indices for dust mapping, converts from RA/Dec to galactic (l,b) coords, and works to exclude bad regions of the catalogs

3. Query dust map: looks up each E(B-V) reddening value at specific points in the sky to create a full dust map

4. Apply extinction model: converts the E(B-V) value using the extinction model to corrections

```bash
python runner_scripts/dust_calc_runner.py -config configs/dust_calc_config.yaml
```

5. Output corrected catalog: creates the deredenned values and outputs them into a new catalog

Note: Steps 2 and 3 happen internally within ``prep_cat_runner.py `` script and only applies to galaxy cross correlation analyses.

## Output
The pipeline above creates new corrected catalogs to the `output/` directory. This directory has the original values as well as:
- Extinction Value: `A_V` 
- Reddening Value: `E(B-V)`
- Corrected values

## Citation 
If you use this code in your research, please cite:

McCleary et al. 2025, "A Detection of Circumgalactic Dust at Megaparsec Scales with Maximum Likelihood Estimation" arXiv:2503.04098
