import numpy as np
import glob
import os
import sys

def add_dl1_paths_to_dict(DICT, dl1_root, dchecking=False):
    """
    Function to add a DICT with run numbers, the respective paths associated to dl1
    for run and subrun
    Input:
    - DICT      --> The dictionary
    - dl1_root  --> The path string where all dl1 data is inside (with possibility of using *)
                    format example: "/fefs/aswg/data/real/DL1/*/v*/tailcut84/"
    - dchecking --> Boolean for extracting the paths for datachecks or for normal files
    """

    str_dchecks = "" if not dchecking else "datacheck_"
    log_errors  = ""
    print(f"\nAdding dl1 {str_dchecks[:-1]} data to dictionary...")

    main_name = f"{str_dchecks}dl1_LST-1.Run?????"
    total_dl1a_runwise    = glob.glob(dl1_root + "*/" + f"{main_name}.h5")      + glob.glob(dl1_root + f"{main_name}.h5")
    total_dl1a_subrunwise = glob.glob(dl1_root + "*/" + f"{main_name}.????.h5") + glob.glob(dl1_root + f"{main_name}.????.h5")
    # print(f"DL1 files:  Found {len(total_dl1a_runwise):4} run-wise and {len(total_dl1a_subrunwise):6} subrun-wise")
    # print(f"Datachecks: Found {len(total_dcheck_runwise):4} run-wise and {len(total_dcheck_subrunwise):6}  subrun_wise\n")

    for run in DICT.keys():
        # checking for files of this certain run
        runfiles = []
        for rf in total_dl1a_runwise:
            if f"{run:05}" in rf:
                runfiles.append(rf)

        # checking runs we have, not, or we have duplicated
        if len(runfiles) > 1:
            warning = f"WARNING: Run {run:5} presented {len(runfiles)} files:"
            log_errors = log_errors + "\n" + warning 
            print(warning)
            for i, runfile in enumerate(runfiles):
                selected = "(SELECTED)" if i == 0 else ""
                warning = f"--> {runfile} {selected}"
                log_errors = log_errors + "\n" + warning 
                print(warning)
        if len(runfiles) == 0:
            error = f"ERROR: Run {run:5} not found in {dl1_root}"
            log_errors = log_errors + "\n" + error
            print(error)
            run_path = None
        else:
            run_path  = runfiles[0]        

        # checking subruns
        subrunfiles = []
        subrun_num  = []
        for f in total_dl1a_subrunwise:
            if str(run) in f:
                subrunfiles.append(f)
                subrun_num.append(int(f.split(".")[-2]))
        subrun_num  = list(np.sort(subrun_num))
        subrunfiles = np.sort(subrunfiles)

        # first we make sure there is at least one subrun
        if len(subrunfiles) > 0:

            # checking we have all subruns, or if some is repeated
            subrun_paths = []
            for i in range(subrun_num[-1]+1):

                try:
                    fi = subrun_num.index(i)
                    n_coincidences = np.sum(np.array(subrun_num) == i)
                    if n_coincidences > 1:
                        
                        warning = f"WARNING: Subrun {i:04} presented {n_coincidences} files:"
                        log_errors = log_errors + "\n" + warning
                        print(warning)
                        
                        for i, subrunfile in enumerate(subrunfiles[np.array(subrun_num) == i]):
                            selected = "(SELECTED)" if i == 0 else ""
                            warning = f"--> {subrunfile} {selected}"
                            log_errors = log_errors + "\n" + warning
                            print(warning)        

                    subrun_paths.append(subrunfiles[fi])
                except ValueError:
                    error = f"ERROR: Subrun {i:04} not found in {dl1_root}"
                    log_errors = log_errors + "\n" + error
                    print(error)
                    subrun_paths.append(None)

        else:
            subrun_paths = []

        # checking if the branch of dict exists
        str_dchecks_dict = "dl1a" if not dchecking else "dchecks"
        try: 
            DICT[run][str_dchecks_dict]["runwise"]    = run_path
            DICT[run][str_dchecks_dict]["srunwise"] = subrun_paths
        except KeyError:
            DICT[run][str_dchecks_dict] = {
                "runwise"  : run_path, 
                "srunwise" : subrun_paths,
            }

        # adding a log error if do not exist in the case that there has been errors
        if log_errors != "":
            try:
                DICT[run]["errors"] = DICT[run]["errors"] + f"\n# error logs dl1 {str_dchecks[:-1]}" + log_errors + "\n"
            except KeyError:
                DICT[run]["errors"] =  f"# error logs dl1 {str_dchecks[:-1]}" + log_errors + "\n"
        else:
            DICT[run]["errors"] = ""
                
    print(f"...Finished adding dl1 data to dictionary")
    return DICT


def angular_dist(az1, az2):
    """
    Angular distance for azimuths where az=0 is the same as az=360
    """
    angular_distance_abs = abs(az1 - az2)
    return min(angular_distance_abs, 360 - angular_distance_abs)

def add_mc_and_rfs_nodes(DICT, rfs_root, mcs_root, dict_source):

    """
    
    """
    
    # finding the RF nodes in DEC
    rfs_decs = os.listdir(rfs_root)
    rf_nodes_dec = np.array([float(d.split("_")[-1][:-2] + "." + d.split("_")[-1][-2:]) for d in rfs_decs])

    dist_rfs = np.abs(rf_nodes_dec - dict_source["dec"].value)
    closest_rf_node = rfs_decs[np.argmin(dist_rfs)]

    # and the MC nodes in AZ ZD
    mcs_decs = os.listdir(mcs_root)
    mc_nodes_dec = [float(d.split("_")[-1][:-2] + "." + d.split("_")[-1][-2:]) for d in mcs_decs]

    dist_mc_dec = np.abs(mc_nodes_dec - dict_source["dec"].value)
    closest_mc_dec_node = rfs_decs[np.argmin(dist_mc_dec)]

    nodes = np.array(os.listdir(mcs_root + mcs_decs[0]))
    zds = np.array([float(n.split("_")[2]) for n in nodes])
    azs = np.array([float(n.split("_")[4]) for n in nodes])

    for run in DICT.keys():
        _zd = DICT[run]["pointing"]["zd"]
        _az = DICT[run]["pointing"]["az"]

        dist_mcs_zd = np.abs(zds - _zd)

        zd_closest_node = zds[np.argmin(dist_mcs_zd)]

        mask_zd  = (zds == zd_closest_node)
        nodes_zd = nodes[mask_zd]
        zds_zd   = zds[mask_zd]
        azs_zd   = azs[mask_zd]

        dist_mcs_az = np.array([angular_dist(azs_zd[i], _az) for i in range(len(azs_zd))])

        closest_node = nodes[np.argmin(dist_mcs_az)]

        # now looking inside the folder to find the MC .h5 file
        mc_fnames = glob.glob(mcs_root + closest_mc_dec_node + "/" + closest_node + "/*.h5")
        if len(mc_fnames) == 0:
            sys.exit("ERROR: no MC files found inside {}".format(mcs_root + closest_mc_dec_node + "/" + closest_node))
        elif len(mc_fnames) > 1:
            print("WARNING: MC path {} presented {} .h5 files:".format(closest_mc_dec_node + "/" + closest_node, len(mc_fnames)))
            for f in mc_fnames:
                selected = "(SELECTED)" if i == 0 else ""
                print(f"--> {f} {selected}")

        mc_fname = mc_fnames[0]

        DICT[run]["simulations"] = {
            "mc" : mc_fname,
            "rf" : rfs_root + closest_rf_node,
        }

    dict_nodes = {
        "dec" : mc_nodes_dec,
        "pointing" : {
            "az" : azs,
            "zd" : zds,
        },
    }

    return DICT, dict_nodes
