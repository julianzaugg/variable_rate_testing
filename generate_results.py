"""
Generate processed results from the output results for the `run_...' scripts.

Output from this script includes:
    - Table summarising the run time for each tool-sequence type(joint/marginal)-jobtype(fixed/variable)
    - Table with run times for each groupsize-repeat-group-seqtype-jobtype
    - Distance matrix for sequences in each group size
    - Table with mean and variability in evolutionary distances for sequences in each group size
    - Distance matrix for all sequences
    - Fasta file with all the sequences
    - Various other distances matrices and sequences file for different groupings of tools/groupsizes/seqtype/jobtype


The phylip package must be installed to PATH. Either the binaries for the different packages need to be added to PATH
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


def _summarise_distances_from_matrix(filename):
    """
    Read a symmetrical distance matrix from protdist output and calculate the mean and variability in distances
    :param filename: distance matrix filename
    :return: name order, distance matrix, mean distance, standard deviation
    """
    names = [] # the name order
    with open(filename, 'r') as fh:
        data = [line.strip() for line in fh.readlines()]
        nseqs = int(data[0])
        dist_matrix = np.zeros((nseqs, nseqs))
        row_idx = -1
        column_idx = 0
        for line in data[1:]:
            splitted_line = line.split()
            if splitted_line[0].isalpha():
                names.append(splitted_line[0])
                column_idx = 0
                row_idx += 1
                for entry in splitted_line[1:]:
                    dist_matrix[row_idx][column_idx] = float(entry)
                    column_idx += 1
            else:
                for entry in splitted_line:
                    dist_matrix[row_idx][column_idx] = float(entry)
                    column_idx += 1
    # Calculate mean and variability
    lower_matrix_greater_0 = np.tril(dist_matrix)[np.where(np.tril(dist_matrix) > 0)]
    mean_dist = np.mean(lower_matrix_greater_0)
    std_dist = np.std(lower_matrix_greater_0)
    return names, dist_matrix, mean_dist, std_dist


def _process_arguments(myparser, myargs):
    input_path = Path(myargs.input).resolve()
    out_dir = Path(myargs.output).resolve().absolute()

    # Find all result directories
    result_leaf_dirs = []
    for root, dirs, files in os.walk(str(input_path)):
        if not dirs:
            result_leaf_dirs.append(root)

    # Collate and write details for each ancestral sequence
    ar_results = []
    for result_dir in result_leaf_dirs: # Loop through each result directory
        splitted_result_dir = result_dir.split("/")
        group_id = splitted_result_dir[-2]
        repeat_id = splitted_result_dir[-3]
        group_size = splitted_result_dir[-4].split("_")[-1]

        # For this to work, job ID must have the name of the tool (GRASP, PAML or FastML)
        # and type of job (variable or fixed) included in the directory name.
        # This should be enforced in the `run_...' scripts.
        job_ID = splitted_result_dir[-1]

        # Get the log filename for the group folder
        base_log_fn = result_dir + "/../" +\
              [base_log_fn for base_log_fn in os.listdir(result_dir + "/../") if "_log" in base_log_fn][0]
        # Determine the number of sequences and columns used for the input
        number_of_seqs, number_of_columns_used = None, None
        with open(base_log_fn, 'r') as fh:
            data = [line.strip() for line in fh.readlines()]
            for line in data:
                if "Nseqs" in line:
                    number_of_seqs = line.split(":")[-1].strip()
                elif "Alignlen_without_gaps" in line:
                    number_of_columns_used = line.split(":")[-1].strip()

        # Figure out which tool the current results belong to based on the job ID in the folder name
        tool = None
        if "grasp" in job_ID.lower():
            tool = "GRASP"
        elif "paml" in job_ID.lower():
            tool = "PAML"
        elif "fastml" in job_ID.lower():
            tool = "FastML"
        # Figure out the job type
        job_type = "fixed" # default to fixed. If "variable" in directory (jobid) name, we change
        if "variable" in job_ID.lower():
            job_type = "variable"

        # Load the marginal and joint ancestral sequences if they are in the result folder
        joint_seq, marginal_seq = None, None
        for result_file in os.listdir(result_dir):
            full_result_file_path = os.path.join(result_dir, result_file)
            if "joint_root_reconstruction" in result_file:
                joint_seq = readFastaFile(full_result_file_path, Protein_Alphabet)[0]
            elif "marginal_root_reconstruction" in result_file:
                marginal_seq = readFastaFile(full_result_file_path, Protein_Alphabet)[0]

        # Try getting the run time from the result log file. If this fails, it probably
        # means the job never completed and so we need to skip this result
        try:
            with open(glob.glob(result_dir + "/log.txt")[0], 'r') as fh:
                time = fh.readline().split()[-2]
        except:
            continue

        # Store the ancestral sequences with the appropriate attribute settings
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
    with open(os.path.join(str(out_dir), myargs.id + "_asr_results_table.txt"), 'w') as fh:
        print("\t".join(["ID","RepeatID","GroupSize","SeqType",
                     "JobType","Time","Tool","SeqNumber","ColumnNumber","Sequence"]), file = fh)
        for ar in ar_results:
            print(ar, file = fh)

    ################################################################################
    # Calculate time differences between calculation of ancestral sequences with different tools and methods

    summary_results_time = defaultdict(lambda : defaultdict(lambda : defaultdict(lambda : defaultdict(list))))
    summary_results_seq_range = defaultdict(lambda : defaultdict(lambda : defaultdict(lambda : defaultdict(list))))

    # Loop through all the ancestral results and collect their run times and the number of sequences
    # used to calculate them
    for ar in ar_results:
        summary_results_time[ar.tool][ar.seq_type][ar.job_type][ar.group_size].append(float(ar.calculation_time))
        summary_results_seq_range[ar.tool][ar.seq_type][ar.job_type][ar.group_size].append(float(ar.number_of_seqs_used))
    # Write the time mean and variance time differances
    with open(os.path.join(str(out_dir), myargs.id + "_asr_results_table_summary.txt"), 'w') as fh:
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
    ########################################
    # Generate combined sequence files from ASR results
    # write out fasta format
    writeFastaFile(os.path.join(str(out_dir), myargs.id + "_asr_seqs.fasta"), [ar.seq_object for ar in ar_results])
    # write out phylip format
    _write_phylip_file(os.path.join(str(out_dir), myargs.id + "_asr_seqs.phylip"), [ar.seq_object for ar in ar_results])

    ########################################
    # Group ASR results by different attributes to calculate evolutionary distance within groups.
    # For groupby to work correctly, the data needs to be sorted by the attribute/s first

    # Group by tool, group size, joint/marginal, fixed/variable.
    # This is a INTRA tool and parameter comparison
    sorted_ar_results = sorted(ar_results, key = lambda x: (x.tool, x.group_size, x.seq_type, x.job_type))
    ar_groups = groupby(sorted_ar_results, lambda x:
                    (x.tool, x.group_size, x.seq_type, x.job_type))
    """
    Calculate the distance between sequences in each group. This requires an input file in phylip format.
    'The output file contains on its first line the number of species. The distance matrix is then printed in standard
    form, with each species starting on a new line with the species name, followed by the distances to the species in
    order. These continue onto a new line after every nine distances. The distance matrix is square with zero
    distances on the diagonal. In general the format of the distance matrix is such that it can serve as input
    to any of the distance matrix programs.'
    """

    # Create the output comparison strings
    group_dist_summary_strs = []
    for ar_g, ar_g_members in ar_groups:
        ar_g_member_seqs = [member.seq_object for member in ar_g_members]
        # Write out the phylip formatted input file for each group
        group_filename = os.path.join(str(out_dir), myargs.id + "_" + "_".join(ar_g) + ".phylip")
        try:
            os.remove(os.path.join(str(out_dir), "outfile"))
        except:
            pass
        _write_phylip_file(group_filename, ar_g_member_seqs)
        _write_phylip_file(os.path.join(str(out_dir), "infile"), ar_g_member_seqs) # write seqs to infile as well
        os.chdir(str(out_dir)) # Protdist looks for infile in the current directory
        subprocess.call("echo Y | protdist", shell = True, stdout=subprocess.PIPE)
        matrix_filename = group_filename.replace(".phylip", ".matrix")
        shutil.move(os.path.join(str(out_dir), "outfile"), matrix_filename) # move the outfile matrix to new file
        name_order, dist_matrix, group_mean_dist, group_std_dist = _summarise_distances_from_matrix(matrix_filename)
        tool = ar_g[0]
        group_size = ar_g[1]
        seq_type = ar_g[2]
        job_type = ar_g[3]
        group_dist_summary_strs.append(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}".format(
                "_".join(ar_g),
                tool,
                group_size,
                seq_type,
                job_type,
                round(group_mean_dist,3),
                round(group_std_dist,3)))
    # Write the distance summaries to file
    with open(os.path.join(str(out_dir), myargs.id + "_group_distance_summaries.txt"), 'w') as fh:
        print("ID\tTool\tGroupSize\tSeqType\tJobType\tMean_distance\tStdev_distance", file = fh)
        for out_str in group_dist_summary_strs:
            print(out_str, file=fh)


    ########################################
    """
    Below, groups contain results for mixtures of tools. We need to know mean distances between sequences from 
    one tool compared to another. The distances between sequences from different tools are stored in the same distance
    matrix. We know the name order of the distance matrix, so we can get the names and map them 
    to a particular tool. Once we have this mapping, we can collect the distances for each tool/method and calculate
    the mean distance.
    """

    # Now group by just groupsize, seq type and job type to compare tools against each other
    # This is an INTER tool comparison
    group_group_dist_summary_strs = []
    sorted_ar_results = sorted(ar_results, key=lambda x: (x.group_size, x.seq_type, x.job_type))
    ar_groups = groupby(sorted_ar_results, lambda x: (x.group_size, x.seq_type, x.job_type))
    for ar_g, ar_g_members in ar_groups:
        ar_g_member_objects = [member for member in ar_g_members]
        ar_g_member_seqs = [member.seq_object for member in ar_g_member_objects]
        # Map sequence names to tool_groupsize_seqtype_jobtype
        mapped_seq_ar = dict(zip(map(lambda x: x.name, ar_g_member_seqs),
                                 map(lambda x: "_".join([x.tool,x.group_size,x.seq_type,x.job_type]), ar_g_member_objects)))

        # Write out the phylip formatted input file for each group
        group_filename = os.path.join(str(out_dir), myargs.id + "_tool_comparison_" + "_".join(ar_g) + ".phylip")
        try:
            os.remove(os.path.join(str(out_dir), "outfile"))
        except:
            pass
        _write_phylip_file(group_filename, ar_g_member_seqs)
        _write_phylip_file(os.path.join(str(out_dir), "infile"), ar_g_member_seqs) # write seqs to infile as well
        os.chdir(str(out_dir)) # Protdist looks for infile in the current directory
        subprocess.call("echo Y | protdist", shell = True, stdout=subprocess.PIPE)
        matrix_filename = group_filename.replace(".phylip", ".matrix")
        shutil.move(os.path.join(str(out_dir), "outfile"), matrix_filename) # move the outfile matrix to new file
        name_order, dist_matrix, group_mean_dist, group_std_dist = _summarise_distances_from_matrix(matrix_filename)

        tool_to_tool_distances = defaultdict(lambda: defaultdict(list))
        lower_matrix = np.tril(dist_matrix) # lower half of the distance matrix
        for row in range(len(lower_matrix)):
            row_seq_name = name_order[row]
            row_seq_group = mapped_seq_ar[row_seq_name]
            for col in range(len(lower_matrix[row])):
                if lower_matrix[row][col] == 0.0: continue
                col_seq_name = name_order[col]
                col_seq_group = mapped_seq_ar[col_seq_name]
                tool_to_tool_distances[row_seq_group][col_seq_group].append(lower_matrix[row][col])

        # Summarise the tool distances
        for g1, v1 in tool_to_tool_distances.items():
            splitted_g1 = g1.split("_")
            tool1 = splitted_g1[0]
            group_size1 = splitted_g1[1]
            seq_type1 = splitted_g1[2]
            job_type1 = splitted_g1[3]
            for g2, v2 in v1.items():
                splitted_g2 = g2.split("_")
                tool2 = splitted_g2[0]
                # Skip self comparison
                if tool1 == tool2: continue
                group_size2 = splitted_g2[1]
                seq_type2 = splitted_g2[2]
                job_type2 = splitted_g2[3]
                group_group_mean_dist = round(np.mean(v2), 3)
                group_group_stdev_dist = round(np.std(v2), 3)
                group_group_dist_summary_strs.append("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(
                    tool1, group_size1, seq_type1, job_type1,
                    tool2, group_size2, seq_type2, job_type2,
                    group_group_mean_dist, group_group_stdev_dist
                    ))

    # Write the tool comparison results to file
    with open(os.path.join(str(out_dir), myargs.id + "_tool_group_distance_comparison.txt"), 'w') as fh:
        print("Tool1\tGroupSize1\tSeqType1\tJobType1\t" + \
              "Tool2\tGroupSize2\tSeqType2\tJobType2\tMean_distance\tStdev_distance", file=fh)
        for ggdss in group_group_dist_summary_strs:
            print (ggdss, file = fh)

    ########################################
    # Now group by just groupsize, to compare tools and methods against each other
    # This is a INTER tools and method comparison
    group_group_dist_summary_strs = []
    sorted_ar_results = sorted(ar_results, key=lambda x: (x.group_size))
    ar_groups = groupby(sorted_ar_results, lambda x: (x.group_size))
    for ar_g, ar_g_members in ar_groups:
        ar_g_member_objects = [member for member in ar_g_members]
        ar_g_member_seqs = [member.seq_object for member in ar_g_member_objects]
        # Map sequence names to tool_groupsize_seqtype_jobtype
        mapped_seq_ar = dict(zip(map(lambda x: x.name, ar_g_member_seqs),
                                 map(lambda x: "_".join([x.tool, x.group_size, x.seq_type, x.job_type]),
                                     ar_g_member_objects)))
        # possible_combinations = set(map(lambda x: "_".join([x.tool, x.group_size, x.seq_type, x.job_type]),
        #                                 ar_g_member_objects))

        # mapped_seq_ar_reverse = dict([(v,[]) for k, v in mapped_seq_ar.items()])

        # Write out the phylip formatted input file for each group
        group_filename = os.path.join(str(out_dir), myargs.id + "_tool_and_method_comparison_group" + ar_g + ".phylip")
        try:
            os.remove(os.path.join(str(out_dir), "outfile"))
        except:
            pass
        _write_phylip_file(group_filename, ar_g_member_seqs)
        _write_phylip_file(os.path.join(str(out_dir), "infile"), ar_g_member_seqs)  # write seqs to infile as well
        os.chdir(str(out_dir))  # Protdist looks for infile in the current directory
        subprocess.call("echo Y | protdist", shell=True, stdout=subprocess.PIPE)
        matrix_filename = group_filename.replace(".phylip", ".matrix")
        shutil.move(os.path.join(str(out_dir), "outfile"), matrix_filename)  # move the outfile matrix to new file
        name_order, dist_matrix, group_mean_dist, group_std_dist = _summarise_distances_from_matrix(matrix_filename)

        tool_to_tool_distances = defaultdict(lambda: defaultdict(list))
        lower_matrix = np.tril(dist_matrix)  # lower half of the distance matrix
        for row in range(len(lower_matrix)):
            row_seq_name = name_order[row]
            row_seq_group = mapped_seq_ar[row_seq_name]
            for col in range(len(lower_matrix[row])):
                if lower_matrix[row][col] == 0.0: continue
                col_seq_name = name_order[col]
                col_seq_group = mapped_seq_ar[col_seq_name]
                tool_to_tool_distances[row_seq_group][col_seq_group].append(lower_matrix[row][col])

        # Summarise the tool distances
        for g1, v1 in tool_to_tool_distances.items():
            splitted_g1 = g1.split("_")
            tool1 = splitted_g1[0]
            group_size1 = splitted_g1[1]
            seq_type1 = splitted_g1[2]
            job_type1 = splitted_g1[3]
            for g2, v2 in v1.items():
                splitted_g2 = g2.split("_")
                tool2 = splitted_g2[0]
                group_size2 = splitted_g2[1]
                seq_type2 = splitted_g2[2]
                job_type2 = splitted_g2[3]
                # Skip self comparison
                if tool1 == tool2 and seq_type1 == seq_type2 and job_type1 == job_type2: continue
                group_group_mean_dist = round(np.mean(v2), 3)
                group_group_stdev_dist = round(np.std(v2), 3)
                group_group_dist_summary_strs.append("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(
                    tool1, group_size1, seq_type1, job_type1,
                    tool2, group_size2, seq_type2, job_type2,
                    group_group_mean_dist, group_group_stdev_dist
                ))
    with open(os.path.join(str(out_dir), myargs.id + "_tool_group_and_method_distance_comparison.txt"), 'w') as fh:
        print("Tool1\tGroupSize1\tSeqType1\tJobType1\t" + \
              "Tool2\tGroupSize2\tSeqType2\tJobType2\tMean_distance\tStdev_distance", file=fh)
        for ggdss in group_group_dist_summary_strs:
            print(ggdss, file=fh)
    try:
        os.remove(os.path.join(str(out_dir), "infile"))
    except:
        pass



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Analyse results generated by ancestral reconstruction tools')
    parser.add_argument('-i', '--input', help='Base directory to search for results', required=True)
    parser.add_argument('-o', '--output', help='Output directory for results', required=True)
    parser.add_argument('-id', help='ID tag for results', required=True)

    result_dir = "/Users/julianzaugg/Documents/University/Phd/Projects/GRASP/Results/CYP2"
    out_dir = "/Users/julianzaugg/Documents/University/Phd/Projects/GRASP/Processed_results/CYP2"
    result_dir = "/Users/julianzaugg/Documents/University/Phd/Projects/GRASP/Results/KARI"
    out_dir = "/Users/julianzaugg/Documents/University/Phd/Projects/GRASP/Processed_results/KARI"
    myargs = ["-i", result_dir,
              "-o", out_dir,
              "-id", "KARI"]
    # args = parser.parse_args(myargs)
    args = parser.parse_args()
    _process_arguments(parser, args)
