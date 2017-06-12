"""
Given an input data directory and input parameters, this script will run GRASP.
Requires GRASP.jar available in same directory as this script.

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
    grasp_jar_path = Path(myargs.grasp).absolute()

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

        # Get input sequence filename
        input_sequence_filename = Path(glob.glob(ad + "/*.fasta")[0]).absolute()
        input_tree_filename = Path([fn for fn in glob.glob(ad + "/*.newick") if "_phylip" not in fn][0]).absolute()
        try:
            os.chdir(str(results_path.absolute()))
            start_time = time.time()
            # < tree_file > < aln_file > < Inference > < ID >
            print("Running GRASP on " + results_path_abs_str)
            if myargs.marginal:
                subprocess.call("java -jar {} {} {} {} {} {}".format(str(grasp_jar_path),
                                                         input_tree_filename,
                                                         input_sequence_filename,
                                                         "Marginal",
                                                         "N0",
                                                         "grasp_marginal"), shell=True)
                total_time = time.time() - start_time

                # Since marginal reconstruction does not output sequences, we need to turn the
                # amino acid distribution matrix into a single consensus sequence
                with open(os.path.join(results_path_abs_str, "grasp_marginal_distribution.txt"), 'r') as fh:
                    marginal_seq_content = []
                    data = [line.strip().split() for line in fh.readlines()]
                    # There is a small bug in the marginal matrix, the number of columns is short by 1 in the header
                    # Alternatively we could just use the length of the frequencies in the line below the header
                    n_cols = int(data[0][-1]) + 2 # +2 to compensate for the bug and the key column
                    keys = [line[0] for line in data[1:-1]] # 0 = header, -1 = gap
                    for col in range(1,n_cols): # loop over each column
                        frequency = map(float, [line[col] for line in data[1:-1]])
                        matched_freqs_dict = dict(zip(keys, frequency))
                        marginal_seq_content.append(max(matched_freqs_dict, key=lambda key: matched_freqs_dict[key]))
                    marginal_seq = Sequence(name ="N0", sequence = ''.join(marginal_seq_content))
                    writeFastaFile(os.path.join(results_path_abs_str,
                                                "grasp_marginal_root_reconstruction.fasta"), [marginal_seq])
            else:
                subprocess.call("java -jar {} {} {} {} {}".format(str(grasp_jar_path),
                                                         input_tree_filename,
                                                         input_sequence_filename,
                                                         "Joint",
                                                         "grasp_joint"), shell=True)
                total_time = time.time() - start_time
                joint_seqs = readFastaFile(os.path.join(results_path_abs_str, "grasp_joint_aln_full.fa"))
                for s in joint_seqs:
                    if s.name == "N0":
                        writeFastaFile(os.path.join(results_path_abs_str, "grasp_joint_root_reconstruction.fasta"),
                                       [s])
                        break

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
    #GRASP_joint_fixed, GRASP_marginal_fixed
    parser = argparse.ArgumentParser(description='Run GRASP to evaluate variable rates')
    parser.add_argument('-i', '--input', help='Base directory for input data.', required=True)
    parser.add_argument('-id', '--id', help='ID tag for run', default="GRASP")
    parser.add_argument('-grasp', help='JAR file for GRASP', required = True)
    parser.add_argument('-m', '--marginal', help='Apply marginal reconstruction', action = 'store_true')

    myargs = ["-id", "grasp_joint",
              "-i", "/Users/julianzaugg/Documents/University/Phd/Projects/GRASP/Data/test_CYP/test_out/CYP/group_size_40",
              "-grasp", "/Users/julianzaugg/Documents/University/Phd/Projects/GRASP/Code/ASR.jar"]

    myargs = ["-id", "grasp_marginal",
              "-i", "/Users/julianzaugg/Documents/University/Phd/Projects/GRASP/Data/test_CYP/test_out/CYP/group_size_40",
              "-grasp", "/Users/julianzaugg/Documents/University/Phd/Projects/GRASP/Code/ASR.jar",
              "-m"]

    myargs = ["-id", "grasp_marginal",
              "-i", "/Users/julianzaugg/Documents/University/Phd/Projects/GRASP/Data/test_KARI/test_out/KARI/group_size_40",
              "-grasp", "/Users/julianzaugg/Documents/University/Phd/Projects/GRASP/Code/ASR.jar",
              "-m"]

    # args = parser.parse_args(myargs)
    args = parser.parse_args()
    _process_arguments(parser, args)