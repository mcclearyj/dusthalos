# dusthalos
Dusty, dusty halos

## Process catalogs
`python runner_scripts/cat_prep_runner.py --config configs/prep_catalog_config.yaml`

## Run reddening calculation
`python runner_scripts/dust_calc_runner.py --config configs/dust_calc_config.yaml`

The `dust_calc_config.yaml` configuration file points to catalog-level configuration files in `configs/catalog_configs`. 
These do useful things like specify the WCS coordinate system of the catalog and the lon/lat column names, the name
of the redshift column, and the path to the catalog itself. 

The parameters of the extinction model to use is also specified in `dust_calc_config.yaml` [here]([url](https://github.com/mcclearyj/dusthalos/blob/4371d0fc41bbef86c4ce33c1e07c525801027854/configs/dust_calc_config.yaml#L38)https://github.com/mcclearyj/dusthalos/blob/4371d0fc41bbef86c4ce33c1e07c525801027854/configs/dust_calc_config.yaml#L38)



