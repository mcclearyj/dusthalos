# Local imports
import src.utils as utils
import os
import src.calculate_source_density as calculate_source_density
import argparse
import numpy as np


def main(args):
    # Read in configuration file
    config_file = args.config
    correl_config = utils.read_yaml(config_file)

    # Create output directory if it doesn't exist
    if not os.path.isdir(correl_config['output_path']):
        os.makedirs(correl_config['output_path'])

    # Set theta_edges (angular bins for source density calculation)
    theta_edges = np.array([0.5, 1, 2, 5, 10, 20, 50, 100, 200])

    # Run the source density calculation
    # NOTE: this function ignores the angular bins set in correl_config 
    # since the calculation is a special-purpose one
    calculate_source_density.run_all(correl_config, theta_edges=theta_edges)

if __name__ == "__main__":   
    parser = argparse.ArgumentParser(description="Runner script for calculating source density. Output is technically a density-ratio profile aka 'boost factor' profile")
    parser.add_argument('-config', '-c', type=str, required=True, help='Path to configuration file')
    args = parser.parse_args()
    main(args)