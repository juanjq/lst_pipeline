import numpy as np
import glob



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
    log_errors  = f"# error logs dl1 {str_dchecks[:-1]}"
    print(f"\nAdding dl1 {str_dchecks[:-1]} data to dictionary...")
    
    total_dl1a_runwise      = glob.glob(dl1_root + f"{str_dchecks}dl1_LST-1.Run?????.h5")
    total_dl1a_subrunwise   = glob.glob(dl1_root + f"{str_dchecks}dl1_LST-1.Run?????.????.h5")
    # print(f"DL1 files:  Found {len(total_dl1a_runwise):4} run-wise and {len(total_dl1a_subrunwise):6} subrun-wise")
    # print(f"Datachecks: Found {len(total_dcheck_runwise):4} run-wise and {len(total_dcheck_subrunwise):6}  subrun_wise\n")

    for run in DICT.keys():
        # checking for files of this certain run
        runfiles = []
        for rf in total_dl1a_runwise:
            if str(run) in rf:
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

        DICT[run][f"dl1a_{str_dchecks}run"]    = run_path
        DICT[run][f"dl1a_{str_dchecks}subrun"] = subrun_paths

        # adding a log error if do not exist in the case that there has been errors
        if log_errors != f"# error logs dl1 {str_dchecks[:-1]}":
            try:
                DICT[run]["errors"] = DICT[run]["errors"] + "\n" + log_errors + "\n"
            except KeyError:
                DICT[run]["errors"] = log_errors + "\n"
                
    print(f"...Finished adding dl1 data to dictionary")
    return DICT