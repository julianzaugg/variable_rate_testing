"""
Given an input data directory and input parameters, this script will run FastML.
Requires FastML installed to path, i.e., fastml.

This can done with brew:
 brew install homebrew/science/fastml

-qf = FASTA

-g Assume among site rate variation model (Gamma) [By default the program   |
 |   will assume an homogenous model. very fast, but less accurate!]
    > -g and no -g


"fastml -s %s -t %s -g -qf" %(alignment[i], tree[i])
"fastml -s %s -t %s -qf" %(alignment[i], tree[i])
"""

import sys, os, argparse
import subprocess
from pathlib import Path
import glob
import shutil
import re
import time

from sequence import *



def _process_arguments(myparser, myargs):
    # Collect directory paths for analysis
    analysis_dirs = []
    init_dir = Path(os.curdir).resolve().absolute()

    input_path = Path(myargs.input).resolve()
    for dir in glob.glob(os.path.join(str(input_path.absolute()), "**"), recursive=True):
        if os.path.isdir(dir):
            base_name = os.path.basename(os.path.normpath(dir))
            if base_name.startswith(("group")) and "size" not in base_name:
                analysis_dirs.append(dir)

    # Sort analysis directories so we iterate from lowest to highest group
    analysis_dirs = sorted(analysis_dirs, key=lambda x: (int(x.split("/")[-2].split("_")[-1]),
                                                         int(x.split("/")[-1].split("_")[-1])))

    for ad in analysis_dirs:
        # Create the results directory
        results_path = Path(os.path.join(ad, myargs.id))
        results_path.mkdir(parents=True, exist_ok=True)
        results_path_abs_str = str(results_path.absolute())

        #Get input sequence filename
        input_sequence_filename = Path(glob.glob(ad + "/*.fasta")[0]).absolute()
        input_tree_filename = Path([fn for fn in glob.glob(ad + "/*.newick") if "_phylip" not in fn][0]).absolute()
        try:
            start_time = time.time()
            print("Running FastML on " + edited_control_filename)
            # echo_str = "echo %s | " %  edited_control_filename

            # cwd = os.getcwd()
            base_command_str = "fastml -mj " \
            "-s ../CYP_gs_40_r_1_g_1.fasta " \
            "-t ../CYP_gs_40_r_1_g_1.newick " \
            "-x ./temp_tree_newick.txt  " \
            "-j ./joint_seqs.txt " \
            "-k ./marginal_seqs.txt " \
            "-qf"
            if myargs.gamma

            subprocess.call("fastml -mj JTT")
            subprocess.call("{}codeml".format(echo_str), shell=True,  stdout=subprocess.PIPE)
            end_time = time.time() - start_time
        except Exception as e:
            print(e)
            pass
        sys.exit()



if __name__ == "__main__":
    #FASTML_joint_fixed, FASTML_marginal_fixed, FASTML_joint_variable, FASTML_marginal_variable

    parser = argparse.ArgumentParser(description='Run FastML to evaluate variable rates')
    parser.add_argument('-id', '--id', help='ID tag for run. Names the result folder.', default = "PAML")
    parser.add_argument('-i', '--input', help='Base directory for input data.', required=True)
    parser.add_argument('-g', '--gamma', help='Apply site rate variation model.', action = 'store_true')

    myargs = ["-id", "fastml_marginal",
              "-i", "/Users/julianzaugg/Documents/University/Phd/Projects/GRASP/Data/test_CYP/test_out/CYP/group_size_40"]

    args = parser.parse_args(myargs)
    _process_arguments(parser, args)
