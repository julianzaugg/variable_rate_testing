"""
Generate processed results from the output results for the `run_...' scripts.

Output from this script includes:
    - Table summarising the run time for each tool-sequence type(joint/marginal)-jobtype(fixed/variable)
    - Table with run times for each groupsize-repeat-group-seqtype-jobtype
    - Distance matrix for sequences in each group size
    - Table with mean and variability in evolutionary distances for sequences in each group size
    - Distance matrix for all sequences
    - Fasta file with all the sequences



The phylip package should be installed to PATH. Either the binaries for the different packages need to be added to PATH
or the package can be installed using :
    brew install homebrew/science/phylip

The specific program that we use in the phylip package is `protdist'. This program takes sequences in the phylip
format and will calculate the distance between them. We use the JTT rate matrix in this script.

"""


import sys, os
import random
import copy
import time
import argparse
from pathlib import Path
import glob
from sequence import *
from collections import defaultdict
import numpy as np
from itertools import groupby
import subprocess
import shutil


def _generateRandomStringID(length=5, alpha = None, seed  = None):
    """
    @param length: length of string
    @return: string
    """
    if seed:
        random.seed(seed)
    if not alpha:
        alpha = string.ascii_uppercase + string.digits
    return  ''.join(random.choice(alpha) for x in range(length))

class AncestorResult():
    """
    An AncestorResult object stores the results for an individual ancestral sequence.
    """

    def __init__(self, ID, repeat_id, group_size, seq_type, job_type, calculation_time, tool, number_of_seqs_used,
                 number_of_columns_used, sequence):
        self.ID = ID # Folder name that results was in.
        self.repeat = repeat_id # Specific repeat
        self.group_size = group_size # Group division size
        self.seq_type = seq_type # marginal or joint
        self.job_type = job_type # variable or fixed
        self.calculation_time = calculation_time # Run time
        self.tool = tool # GRASP, PAML or FastML
        self.number_of_seqs_used = number_of_seqs_used
        self.number_of_columns_used = number_of_columns_used
        self.sequence = sequence
        self.seq_object = Sequence(name = _generateRandomStringID(10, Protein_Alphabet),
                              sequence = self.sequence, info=self.get_fasta_info())


    def __str__(self):
        return "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(
                                            self.ID,
                                            self.repeat,
                                            self.group_size,
                                            self.seq_type,
                                            self.job_type,
                                            self.calculation_time,
                                            self.tool,
                                            self.number_of_seqs_used,
                                            self.number_of_columns_used,
                                            self.sequence)

    def get_fasta_info(self):
        return ' '.join([self.ID, self.repeat, self.group_size, self.seq_type, self.job_type, self.tool])


def _phylip_padding(name, width=8):
    return name + "_" * (width - len(name))

def _write_phylip_file(filename, seqs):
    """
    Write sequences as phylip format that is readable by the PHYLIP programs.
    :param filename:
    :param seqs:
    :return:
    """
    with open(filename, 'w') as fh:

        print("   {}   {}".format(len(seqs), len(seqs[0])), file = fh)
        for seq in seqs:
            window = 0
            chunks = []
            while window < len(seq):
                chunks.append(seq[window:window + 10])
                window += 10
            formatted_seq = " ".join(chunks)
            print ("{:10}{}".format(_phylip_padding(seq.name, 8), formatted_seq), file = fh, sep = "\n")



def _process_arguments(myparser, myargs):
    """
    Calculate:
        ID  Sequence    Type(Marginal/Joint)    TimeToRun
        ID  GroupSize   Type(Marginal/Joint)    MeanDist    StdevDist

            ID1 ID2 ID3
        ID1
        ID2
        ID3
    """
    input_path = Path(myargs.input).resolve()
    out_dir = Path(myargs.output).resolve().absolute()
    result_leaf_dirs = []
    for root, dirs, files in os.walk(str(input_path)):
        if not dirs:
            result_leaf_dirs.append(root)

    # Collate and write details for each ancestral sequence
    ar_results = []
    for result_dir in result_leaf_dirs:
        splitted_result_dir = result_dir.split("/")
        group_id = splitted_result_dir[-2]
        repeat_id = splitted_result_dir[-3]
        group_size = splitted_result_dir[-4].split("_")[-1]

        # For this to work, job ID must have the name of the tool (GRASP, PAML or FastML)
        # and type of job (variable or fixed) included in the directory name.
        # This should be enforced in the `run_...' scripts.
        job_ID = splitted_result_dir[-1]

        base_log_fn = result_dir + "/../" +\
              [base_log_fn for base_log_fn in os.listdir(result_dir + "/../") if "_log" in base_log_fn][0]
        number_of_seqs, number_of_columns_used = None, None
        with open(base_log_fn, 'r') as fh:
            data = [line.strip() for line in fh.readlines()]
            for line in data:
                if "Nseqs" in line:
                    number_of_seqs = line.split(":")[-1].strip()
                elif "Alignlen_without_gaps" in line:
                    number_of_columns_used = line.split(":")[-1].strip()

        tool = None
        if "grasp" in job_ID.lower():
            tool = "GRASP"
        elif "paml" in job_ID.lower():
            tool = "PAML"
        elif "fastml" in job_ID.lower():
            tool = "FastML"
        job_type = "fixed" # default to fixed. If "variable" in directory (jobid) name, we change
        if "variable" in job_ID.lower():
            job_type = "variable"
        joint_seq, marginal_seq = None, None
        for result_file in os.listdir(result_dir):
            full_result_file_path = os.path.join(result_dir, result_file)
            if "joint_root_reconstruction" in result_file:
                joint_seq = readFastaFile(full_result_file_path, Protein_Alphabet)[0]
            elif "marginal_root_reconstruction" in result_file:
                marginal_seq = readFastaFile(full_result_file_path, Protein_Alphabet)[0]
        with open(glob.glob(result_dir + "/log.txt")[0], 'r') as fh:
            time = fh.readline().split()[-2]


        if joint_seq:
            joint_seq_result = \
                AncestorResult(
                    job_ID,
                    repeat_id,
                    group_size,
                    "joint",
                    job_type,
                    time,
                    tool,
                    number_of_seqs_used = number_of_seqs,
                    number_of_columns_used = number_of_columns_used,
                    sequence="".join(joint_seq.sequence))
            ar_results.append(joint_seq_result)
        if marginal_seq:
            marginal_seq_result = \
                AncestorResult(
                    job_ID,
                    repeat_id,
                    group_size,
                    "marginal",
                    job_type,
                    time,
                    tool,
                    number_of_seqs_used = number_of_seqs,
                    number_of_columns_used = number_of_columns_used,
                    sequence="".join(marginal_seq.sequence))
            ar_results.append(marginal_seq_result)

    ################################################################################
    # Write full result table for time differences
    with open(os.path.join(str(out_dir), "asr_results_table.txt"), 'w') as fh:
        print("\t".join(["ID","RepeatID","GroupSize","SeqType",
                     "JobType","Time","Tool","SeqNumber","ColumnNumber","Sequence"]), file = fh)
        for ar in ar_results:
            print(ar, file = fh)

    ######################################
    # Write summary results for time differences
    summary_results_time = defaultdict(lambda : defaultdict(lambda : defaultdict(lambda : defaultdict(list))))
    summary_results_seq_range = defaultdict(lambda : defaultdict(lambda : defaultdict(lambda : defaultdict(list))))
    for ar in ar_results:
        summary_results_time[ar.tool][ar.seq_type][ar.job_type][ar.group_size].append(float(ar.calculation_time))
        summary_results_seq_range[ar.tool][ar.seq_type][ar.job_type][ar.group_size].append(float(ar.number_of_seqs_used))
    with open(os.path.join(str(out_dir), "asr_results_table_summary.txt"), 'w') as fh:
        print("Tool\tSeqType\tJobType\tGroupSize\tSeqRange\tTime_Mean\tTime_SD", file = fh)
        for tool, seq_type_data in summary_results_time.items():
            for seq_type, job_type_data in seq_type_data.items():
                for job_type, group_size_data in job_type_data.items():
                    for group_size, times in group_size_data.items():
                        min_seq = int(min(summary_results_seq_range[tool][seq_type][job_type][group_size]))
                        max_seq = int(max(summary_results_seq_range[tool][seq_type][job_type][group_size]))
                        seq_range_str = "%i-%i" % (min_seq, max_seq)
                        print ("{}\t{}\t{}\t{}\t{}\t{}\t{}".format(
                            tool, seq_type, job_type, group_size,
                            seq_range_str, round(np.mean(times),3), round(np.std(times)),3), file = fh)
    ################################################################################
    # Calculate evolutionary distances between ancestral sequences:
    #   - Within the same group
    #   - All vs All
    ########################################
    # Generate combined sequence files from ASR results
    # write out fasta format
    writeFastaFile(os.path.join(str(out_dir), "asr_seqs.fasta"), [ar.seq_object for ar in ar_results])
    # write out phylip format
    _write_phylip_file(os.path.join(str(out_dir), "asr_seqs.phylip"), [ar.seq_object for ar in ar_results])

    ########################################
    # Group ASR results by attributes to calculate evolutionary distance within groups.
    # For groupby to work correctly, the data needs to be sorted by the attribute/s first

    # Group by tool, group size, joint/marginal, fixed/variable.
    sorted_ar_results = sorted(ar_results, key = lambda x: (x.tool, x.group_size, x.seq_type, x.job_type))
    ar_groups = groupby(sorted_ar_results, lambda x:
                    (x.tool, x.group_size, x.seq_type, x.job_type))
    # Calculate the distance between sequences in each group. This requires an input file in phylip format.
    # "The output file contains on its first line the number of species. The distance matrix is then printed in standard
    # form, with each species starting on a new line with the species name, followed by the distances to the species in
    # order. These continue onto a new line after every nine distances. The distance matrix is square with zero
    # distances on the diagonal. In general the format of the distance matrix is such that it can serve as input
    # to any of the distance matrix programs."
    for ar_g, ar_g_members in ar_groups:
        ar_g_member_seqs = [member.seq_object for member in ar_g_members]
        # Write out the phylip formatted input file for each group
        group_filename = os.path.join(str(out_dir), "_".join(ar_g) + ".phylip")
        try:
            os.remove(os.path.join(str(out_dir), "outfile"))
        except:
            pass
        _write_phylip_file(group_filename, ar_g_member_seqs)
        _write_phylip_file(os.path.join(str(out_dir), "infile"), ar_g_member_seqs)
        os.chdir(str(out_dir))
        subprocess.call("echo Y | protdist", shell = True)
        shutil.move(os.path.join(str(out_dir), "outfile"), group_filename.replace(".phylip", ".matrix"))
        # Have to save the phylip file as `infile' as well
        # Call protdist with `echo Y | protdist'
        # Have to be in same directory when calling
        # rename outfile as the distance matrix file
        # clean up at end

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Analyse results generated by ancestral reconstruction tools')
    parser.add_argument('-i', '--input', help='Base directory to search for results', required=True)
    parser.add_argument('-o', '--output', help='Output directory for results', required=True)

    result_dir = "/Users/julianzaugg/Documents/University/Phd/Projects/GRASP/Results/CYP2"
    out_dir = "/Users/julianzaugg/Documents/University/Phd/Projects/GRASP/Processed_results/CYP2"
    myargs = ["-i", result_dir,
              "-o", out_dir]
    args = parser.parse_args(myargs)
    # args = parser.parse_args()
    _process_arguments(parser, args)
