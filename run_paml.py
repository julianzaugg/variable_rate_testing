"""
Given an input data directory and a template control file (.ctl), this script will run PAML.
Requires PAML installed to path, i.e., codeml.

This can done with brew:
 brew install homebrew/science/paml

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
    # TODO: apply hard limit option to stop trying to process jobs with too few/many sequences?
    # Collect directory paths for analysis
    analysis_dirs = []
    init_dir = Path(os.curdir).resolve().absolute()
    if "paml" not in myargs.id.lower():
        raise RuntimeError("Job ID must include the word 'paml'")
    if "variable" not in myargs.id.lower() and "fixed" not in myargs.id.lower():
        raise RuntimeError("Job ID must include the word 'fixed' or 'variable'")

    input_path = Path(myargs.input).resolve()
    for dir in glob.glob(os.path.join(str(input_path.absolute()), "**"), recursive=True):
        if os.path.isdir(dir):
            base_name = os.path.basename(os.path.normpath(dir))
            if base_name.startswith(("group")) and "size" not in base_name:
                analysis_dirs.append(dir)

    # Sort analysis directories so we iterate from lowest to highest group
    analysis_dirs = sorted(analysis_dirs, key = lambda x: (int(x.split("/")[-2].split("_")[-1]),
                                                 int(x.split("/")[-1].split("_")[-1])
                                                 ))
    for ad in analysis_dirs:

        # Create the results directory
        results_path = Path(os.path.join(ad, myargs.id))
        results_path.mkdir(parents=True, exist_ok=True)

        results_path_abs_str = str(results_path.absolute())

        #Get input sequence filename
        input_sequence_filename = Path(glob.glob(ad + "/*.phylip")[0]).absolute()
        input_tree_filename = Path(glob.glob(ad + "/*_phylip.newick")[0]).absolute()
        #open template control file
        new_lines = []
        with open(str(Path(myargs.control).absolute()), 'r') as fh:
            data = [line.strip() for line in fh.readlines()]
            for line in data:
                if "seqfile" in line:
                    splitted_line = line.split()
                    to_replace = splitted_line[2]
                    line = line.replace(to_replace, str(input_sequence_filename))
                elif "treefile" in line:
                    splitted_line = line.split()
                    to_replace = splitted_line[2]
                    line = line.replace(to_replace, str(input_tree_filename))
                elif "aaRatefile" in line:
                    splitted_line = line.split()
                    to_replace = splitted_line[2]
                    line = line.replace(to_replace, str(Path(myargs.model).absolute()))
                new_lines.append(line)
        # write edited control file
        edited_control_filename = os.path.basename(
                                    os.path.normpath(
                                    str(input_sequence_filename))
                                    ).split(".")[0] + "_paml.ctl"
        edited_control_filename = os.path.join(results_path_abs_str, edited_control_filename)
        with open(edited_control_filename, 'w') as fh:
            print("\n".join(new_lines), file = fh)
        try:
            start_time = time.time()
            # Run PAML
            print("Running PAML on " + edited_control_filename)
            echo_str = "echo %s | " %  edited_control_filename

            # Change the current working directory to the results location. This is necessary as PAML generates
            # files wherever it is run from. If running multiple jobs with this script, generated files would
            # be overwritten if we generated them at the location of this script.
            os.chdir(str(results_path.absolute()))
            # cwd = os.getcwd()
            # subprocess.call("{}codeml".format(echo_str), shell=True,  stdout=subprocess.PIPE)
            subprocess.call("{}codeml".format(echo_str), shell=True)

            end_time = time.time()
            total_time = end_time - start_time

            # Load the rst file and pull out the marginal and, if applicable, joint reconstruction for the root node
            # and write to a new separate file.
            # The internal nodes are lines starting with `node', followed by #N where N is the node number. The lowest
            # is the root node. The first one encountered in the file will be marginal and the second the joint.

            # Joint and marginal reconstructions are found here
            rst_filename = os.path.join(results_path_abs_str, "rst")
            marginal_reconstructions = []
            joint_reconstructions = []
            names_seen = set()
            with open(rst_filename, 'r') as fh:
                data = [line.strip() for line in fh.readlines()]
                for line in data:
                    if line.startswith("node"):
                        # print(line)
                        match = re.compile("\s{5,}").split(line)
                        name = "_".join(match[0].split())
                        content = "".join(match[1].split())
                        if name in names_seen:
                            joint_reconstructions.append(Sequence(name = name, sequence=content))
                        else:
                            marginal_reconstructions.append(Sequence(name=name, sequence=content))
                        names_seen.add(name)

            # Get the root reconstructions and write to file
            ## internal tag is simply the details encoding the group size, repeat number and current group
            ## TODO: apply this tag?
            internal_tag = re.search("gs_\d*_r_\d*_g_\d*", edited_control_filename).group()

            # sort the reconstructions as the lowest number node is the root node
            root_marginal_reconstruction = sorted(marginal_reconstructions, key=lambda x: int(x.name.split("#")[-1]))[0]
            writeFastaFile(os.path.join(results_path_abs_str, "marginal_root_reconstruction.txt"),
                           [root_marginal_reconstruction])
            writeFastaFile(os.path.join(str(results_path.absolute()), "marginal_reconstructions.txt"), marginal_reconstructions)

            # Joint is not always generated. If this fails the program will continue due to the error handling
            root_joint_reconstruction = sorted(joint_reconstructions, key=lambda x: int(x.name.split("#")[-1]))[0]
            writeFastaFile(os.path.join(results_path_abs_str, "joint_root_reconstruction.txt"),
                           [root_joint_reconstruction])
            writeFastaFile(os.path.join(results_path_abs_str, "joint_reconstructions.txt"), joint_reconstructions)

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
    # PAML_joint_fixed, PAML_marginal_fixed, PAML_marginal_variable

    parser = argparse.ArgumentParser(description='Run PAML to evaluate variable rates')
    parser.add_argument('-id', '--id', help='ID tag for run. Names the result folder.  Must include name of tool '
                                            'and either "variable" or "fixed"', required=True)
    parser.add_argument('-m', '--model', help='Location of model file (aaRatefile) for PAML', required=True)
    parser.add_argument('-c', '--control', help='Template control (.ctl) file. Only input'
                                                'sequence and tree entries are changed.', required=True)
    parser.add_argument('-i', '--input', help='Base directory for input data.', required=True)

    # myargs = ["-id", "paml_marginal_joint_fixed",
    #           "-m", "/Users/julianzaugg/Documents/University/Phd/Projects/GRASP/Data/jones.dat",
    #           "-c", "/Users/julianzaugg/Documents/University/Phd/Projects/GRASP/Data/marginal_joint_fixed_codeml_ctl.txt",
    #           "-i", "/Users/julianzaugg/Documents/University/Phd/Projects/GRASP/Data/test_CYP/test_out/CYP/group_size_40"]

    # myargs = ["-id", "paml_marginal_variable",
    #           "-m", "/Users/julianzaugg/Documents/University/Phd/Projects/GRASP/Data/jones.dat",
    #           "-c", "/Users/julianzaugg/Documents/University/Phd/Projects/GRASP/Data/marginal_variable_codeml_ctl.txt",
    #           "-i", "/Users/julianzaugg/Documents/University/Phd/Projects/GRASP/Data/test_CYP/test_out/CYP/group_size_3"]

    # myargs = ["-id", "paml_marginal_joint_fixed",
    #           "-m", "/Users/julianzaugg/Documents/University/Phd/Projects/GRASP/Data/jones.dat",
    #           "-c", "/Users/julianzaugg/Documents/University/Phd/Projects/GRASP/Data/marginal_joint_fixed_codeml_ctl.txt",
    #           "-i", "/Users/julianzaugg/Documents/University/Phd/Projects/GRASP/Data/test_KARI/test_out/KARI/group_size_10"]

    # args = parser.parse_args(myargs)
    args = parser.parse_args()
    _process_arguments(parser, args)