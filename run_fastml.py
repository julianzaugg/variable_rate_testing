"""
Given an input data directory and input parameters, this script will run FastML.
Requires FastML installed to path, i.e., fastml.

This can done with brew:
 brew install homebrew/science/fastml
"""

import sys
import os
import argparse
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
            os.chdir(str(results_path.absolute()))
            start_time = time.time()
            print("Running FastML on " + results_path_abs_str)
            # cwd = os.getcwd()
            base_command_str = "fastml -mj " \
            "-s {} " \
            "-t {} " \
            "-x {} " \
            "-j {} " \
            "-k {} " \
            "-y {} " \
            "-e {} " \
            "-d {} " \
            "-qf".format(input_sequence_filename, # -s
                         input_tree_filename, # -t
                         Path(os.path.join(results_path_abs_str, "fastml_reconstruction_tree.newick")).absolute(), # -x
                         Path(os.path.join(results_path_abs_str, "fastml_joint_seqs.fasta")).absolute(), # -j
                         Path(os.path.join(results_path_abs_str, "fastml_marginal_seqs.fasta")).absolute(), # k
                         Path(os.path.join(results_path_abs_str,
                                           "fastml_reconstruction_tree_ancestor_format.txt")).absolute(), # -y
                         Path(os.path.join(results_path_abs_str, "marginal_probs.txt")).absolute(), # -e
                         Path(os.path.join(results_path_abs_str, "joint_probs.txt")).absolute()  # -d
                         )

            if myargs.gamma:
                subprocess.call(base_command_str + " -g", shell=True)
            else:
                subprocess.call(base_command_str, shell=True)
            total_time = time.time() - start_time

            # Pull out just the root
            joint_seqs = readFastaFile(os.path.join(results_path_abs_str, "fastml_joint_seqs.fasta"))
            marginal_seqs = readFastaFile(os.path.join(results_path_abs_str, "fastml_marginal_seqs.fasta"))
            for s in joint_seqs:
                if s.name == "N1":
                    writeFastaFile(os.path.join(results_path_abs_str, "fastml_joint_root_reconstruction.fasta"),
                                   [s])
                    break
            for s in marginal_seqs:
                if s.name == "N1":
                    writeFastaFile(os.path.join(results_path_abs_str, "fastml_marginal_root_reconstruction.fasta"),
                                   [s])
                    break

            # FastML generates its own log file. Rename it to make way for our log file.
            os.rename(os.path.join(results_path_abs_str, "log.txt"),
                      os.path.join(results_path_abs_str, "fastml_log.txt"))
            # Remove the marginal probs file. Individually they are not large, but if analysing many groups across
            # multiple repeats, this can add up.
            os.remove(os.path.join(results_path_abs_str, "marginal_probs.txt"))

        except Exception as e:
            print(e)
            pass

        # Write log file for the run, i.e time to run program
        with open(os.path.join(results_path_abs_str, "log.txt"), 'w') as fh:
            print("Time to run : {} seconds".format(str(round(total_time, 3))), file = fh)
        # Change directory back to initial directory. This `should' address issues with
        # user providing relative file paths as input parameters.
        os.chdir(str(init_dir))
        # sys.exit()



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Run FastML to evaluate variable rates.')
    parser.add_argument('-id', '--id', help='ID tag for run. Names the result folder. Must include name of tool '
                                            'and either "variable" or "fixed"',required=True)
    parser.add_argument('-i', '--input', help='Base directory for input data.', required=True)
    parser.add_argument('-g', '--gamma', help='Apply site rate variation gamma model.', action = 'store_true')

    myargs = ["-id", "fastml_fixed",
              "-i", "/Users/julianzaugg/Documents/University/Phd/Projects/GRASP/Data/test_CYP/test_out/CYP/group_size_40"]

    # myargs = ["-id", "fastml_variable",
    #           "-i", "/Users/julianzaugg/Documents/University/Phd/Projects/GRASP/Data/test_CYP/test_out/CYP/group_size_40",
    #           "-g"]

    myargs = ["-id", "fastml_fixed",
        "-i", "/Users/julianzaugg/Documents/University/Phd/Projects/GRASP/Data/test_KARI/test_out/KARI/group_size_40"]
    # myargs = ["-id", "fastml_variable",
    #     "-i", "/Users/julianzaugg/Documents/University/Phd/Projects/GRASP/Data/test_KARI/test_out/KARI/group_size_40",
    #           '-g']
    # args = parser.parse_args(myargs)
    args = parser.parse_args()
    _process_arguments(parser, args)
