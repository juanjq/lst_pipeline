# Importing necessary libraries
import matplotlib.pyplot as plt
import numpy as np
import astropy.units as u
from datetime import datetime
import pickle, json, sys, os, glob
import pandas as pd
from astropy.coordinates import SkyCoord

# Importing Lstchain packages
from traitlets.config.loader import Config
from lstchain.io.config import get_standard_config
from ctapipe.io import read_table

import logging
logger = logging.getLogger(__name__)
logger.addHandler(logging.StreamHandler())
logger.setLevel(logging.INFO)

# Importing custom utility functions
import utils


"""
Standard paths for the data in the IT cluster
"""
# DL1 data root
dl1_root = "/fefs/aswg/data/real/DL1/*/v0.*/tailcut84/"
# RFs root main directory
rfs_root = "/fefs/aswg/data/models/AllSky/20230901_v0.10.4_allsky_base_prod/"
# MCs dl2 main directory
mcs_root = "/fefs/aswg/data/mc/DL2/AllSky/20230901_v0.10.4_allsky_base_prod/TestingDataset/"


def main(run_number, source_name, root_data, subruns_num=None):
    """
    Main function for the analysis of the data.
    
    Parameters
    ----------
    run_number : int
        Run number.
    source_name : str
        Source name.
    root_data : str
        Path to the root data.
    subruns_num : int, optional
        Number of subruns. The default is None, that means all the subruns
    """

    # Root path of this script
    root = os.getcwd() + "/"
    # Path to store the configuration file we are going to use
    config_file = root + "config/standard_config.json"
    # Data main directory
    root_data = os.path.join(root_data, source_name)

    run_numbers = [run_number]

    # Creating the directories in case they don't exist
    # Directories for the data
    dir_dl1b = root_data + "dl1b/"
    dir_dl1m = root_data + "dl1_merged/"
    dir_dl2  = root_data + "dl2/"
    dir_dl3  = root_data + "dl3/"
    dir_irfs = root_data + "irfs/"
    for path in [os.path.dirname(config_file), dir_dl1b, dir_dl2, dir_dl1m, dir_dl3, dir_irfs]:
        if not os.path.exists(path):
            os.makedirs(os.path.join(path), exist_ok=True)

    

    # Creating the config file
    config_dict = get_standard_config()
    # Changes in the configuration should be done here
    # We select the heuristic flatfield option in the standard configuration
    config_dict["source_config"]["LSTEventSource"]["use_flatfield_heuristic"] = True
    with open(config_file, 'w') as json_file:
        json.dump(config_dict, json_file)

    # Creating a directory dictionary
    # Getting coordinates of source
    source_coords = SkyCoord.from_name(source_name)
    dict_source = {
        "name"   : source_name,
        "coords" : source_coords,
        "ra"     : source_coords.ra.deg  * u.deg, # ra in degrees
        "dec"    : source_coords.dec.deg * u.deg, # dec in degrees
    }
    # We create a empty dictionary to store all the information needed inside
    DICT = {}
    for run in run_numbers:
        DICT[run] = {
            "run_num" : run
        }
    DICT = utils.add_dl1_paths_to_dict(DICT, dl1_root)
    DICT = utils.add_dl1_paths_to_dict(DICT, dl1_root, dchecking=True)

    # Adding info from the datachecks
    for run in run_numbers:
        tab = read_table(DICT[run]["dchecks"]["runwise"], "/dl1datacheck/cosmics")
        # reading the variables
        _zd, _az = 90 - np.rad2deg(np.array(tab["mean_alt_tel"])), np.rad2deg(np.array(tab["mean_az_tel"]))
        _t_start, _t_elapsed = tab["dragon_time"][0][0], np.array(tab["elapsed_time"])
        DICT[run]["time"] = {
            "tstart"   : _t_start,            # datetime object
            "telapsed" : np.sum(_t_elapsed),  # s
            "srunwise" : {
                "telapsed" : _t_elapsed,      # s      
            },
        }
        DICT[run]["pointing"] = {
            "zd" : np.mean(_zd),  # deg
            "az" : np.mean(_az),  # deg
            "srunwise" : {
                "zd" : _zd,       # deg
                "az" : _az,       # deg
            },
        }
    # then we also select the RFs and MC files looking at the nodes available
    DICT, dict_nodes = utils.add_mc_and_rfs_nodes(DICT, rfs_root, mcs_root, dict_source)


    """
    dl1a_to_dl1b
    """
    for ir, run in enumerate(DICT.keys()):
    
        dir_run = dir_dl1b + f"{run:05}" + "/"
        sruns = [int(path.split(".")[-2]) for path in DICT[run]["dl1a"]["srunwise"]]
        DICT[run]["dl1b"] = {"srunwise" : []}
        
        # Create a folder for each run
        if not os.path.exists(dir_run):
            os.makedirs(os.path.join(dir_run), exist_ok=True)
    
        for i, srun in enumerate(sruns[:subruns_num]):
    
            input_fname  = DICT[run]["dl1a"]["srunwise"][srun]
            output_fname = dir_run + f"dl1_LST-1.Run{run:05}.{srun:04}.h5"
    
            logger.info(f"\nComputing dl1b Run {run:5} Subrun {srun:04} - {i/len(sruns)*100:3.1f}% sruns {ir+1}/{len(DICT.keys())} runs")
            logger.info(f"--> {output_fname}\n")
    
            !lstchain_dl1ab \
              --input-file $input_fname \
              --output-file $output_fname \
              --config $config_file \
              --no-image
    
            DICT[run]["dl1b"]["srunwise"].append(output_fname)


    """
    dl1_merging
    """
    for ir, run in enumerate(DICT.keys()):
    
        dir_run = dir_dl1b + f"{run:05}" + "/"
        output_fname = dir_dl1m + f"dl1_LST-1.Run{run:05}.h5"
        
        !lstchain_merge_hdf5_files \
          --input-dir $dir_run \
          --output-file $output_fname \
          --run-number $run \
          --no-image
        
        DICT[run]["dl1b"]["runwise"] = output_fname


    """
    dl1b_to_dl2
    """
    for ir, run in enumerate(DICT.keys()):
    
        input_fname  = DICT[run]["dl1b"]["runwise"]
        output_fname = dir_dl2 + input_fname.split("/")[-1].replace("dl1", "dl2", 1)
        rf_node      = DICT[run]["simulations"]["rf"]
    
        # Check if the file exists and delete if exists (may be empty or half filled)
        if os.path.exists(output_fname):
            logger.warning(f"File already exists, deleting and re-computing:\n-->{output_fname}")
            os.remove(output_fname)
    
        logger.info(f"\nComputing dl2 from dl1b Run {run:5} - {ir+1}/{len(DICT.keys())} runs")
        logger.info(f"--> {output_fname}\n")
    
        !lstchain_dl1_to_dl2 \
          --input-files $input_fname \
          --path-models $rf_node \
          --output-dir $dir_dl2 \
          --config $config_file
            
        DICT[run]["dl2"] = output_fname

    """
    mcs_to_irfs
    """
    # Already computed IRFs
    computed_irfs = glob.glob(dir_irfs + "*")
    
    for ir, run in enumerate(DICT.keys()):
        
        input_mc = DICT[run]["simulations"]["mc"]
    
        output_irf = dir_irfs + "irf_{}_{}.fits.gz".format(input_mc.split("/")[-3], input_mc.split("/")[-2])
    
        # we don't compute the IRF if it has been already done
        if output_irf not in computed_irfs:
            
            logger.info(f"\nComputing IRF for Run {run:5}, {ir+1}/{len(DICT.keys())} runs")
            logger.info(f"--> {output_irf}\n")
            
            !lstchain_create_irf_files \
              --input-gamma-dl2 $input_mc \
              --output-irf-file $output_irf \
              --point-like \
              --energy-dependent-gh \
              --energy-dependent-theta \
              --overwrite
        else:
            print("\nIRF {}_{} already computed\n".format(input_mc.split("/")[-3], input_mc.split("/")[-2]))
        DICT[run]["irf"] = output_irf

    """
    dl2_to_dl3
    """
    ra_str  = "{}".format(dict_source["ra"]).replace(" ", "")
    dec_str = "{}".format(dict_source["dec"]).replace(" ", "")
    
    for ir, run in enumerate(DICT.keys()):
    
        # dir_run = dir_dl3 + f"{run:05}" + "/"    
        dl2_fname = DICT[run]["dl2"]
    
        output_dl3 = dir_dl3 + f"dl3_LST-1.Run{run:05}.fits"
        
        logger.info(f"\nConverting dl2 for {run:5}, {ir+1}/{len(DICT.keys())} runs")
        logger.info(f"--> {output_dl3}\n")
        
        !lstchain_create_dl3_file \
          --input-dl2 $dl2_fname \
          --input-irf-path $dir_irfs \
          --output-dl3-path $dir_dl3 \
          --source-name $source_name \
          --source-ra $ra_str \
          --source-dec $dec_str \
          --config $config_file \
          --overwrite
    
        DICT[run]["dl3"] = output_dl3

    """
    index files
    """
    logger.info(f"All dl3 files created 100%\n\n\nCreating index files...")
    
    # Creating the index file
    !lstchain_create_dl3_index_files \
      --input-dl3-dir $dir_dl3 \
      --file-pattern 'dl3*.fits' \
      --overwrite
    
    logger.info(f"\nFinished with the dl3 process")


if __name__ == "__main__":
    run_num     = sys.argv[1]
    source_name = sys.argv[2]
    root_data   = sys.argv[3]
    main(run_num, source_name, root_data)