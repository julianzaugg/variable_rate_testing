"""
Given an input data directory and input parameters, this script will run FastML.
Requires PAML installed to path, i.e., codeml.


-qf = FASTA

-g Assume among site rate variation model (Gamma) [By default the program   |
 |   will assume an homogenous model. very fast, but less accurate!]
    > -g and no -g


"fastml -s %s -t %s -g -qf" %(alignment[i], tree[i])
"fastml -s %s -t %s -qf" %(alignment[i], tree[i])

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
                        ancestor sequences
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
    #FASTML_joint_fixed, FASTML_marginal_fixed, FASTML_joint_variable, FASTML_marginal_variable


    parser = argparse.ArgumentParser(description='Run FastML to evaluate variable rates')
    parser.add_argument('-id', '--id', help='ID tag for run', default="FastML")

    myargs = []

    args = parser.parse_args(myargs)
    _process_arguments(parser, args)
