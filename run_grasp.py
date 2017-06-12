"""
Given an input data directory and input parameters, this script will run GRASP.
Requires GRASP.jar available in same directory as this script.


output
    groupsize
        repeat
            GRASP
                marginal
                    logfile (time it took to run, number of sequences, length of sequences)
                    ancestor sequences
                joint
                    ancestor sequences
            PAML
                marginal
                    variable_rates
                        logfile
                    not_variable_rates
                        logfile
                        ancestor sequences
            FASTML
                logfile
"""

import sys, os, argparse



def _process_arguments(myparser, myargs):
    pass





if __name__ == "__main__":
    #GRASP_joint_fixed, GRASP_marginal_fixed
    parser = argparse.ArgumentParser(description='Run GRASP to evaluate variable rates')
    parser.add_argument('-id', '--id', help='ID tag for run', default="GRASP")

    myargs = []

    args = parser.parse_args(myargs)
    _process_arguments(parser, args)