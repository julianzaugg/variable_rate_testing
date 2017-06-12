"""
This script will generate the input data required for run_paml.py, run_fastml.py and run_grasp.py.
"""

import sys, os
import random
import copy
import time
import argparse
from pathlib import Path


import Bio
from Bio import Phylo
from Bio import SeqIO

from sequence import *



##################################################
#For PAML we require sequences to be in phylip-sequential format (with 2 spaces between 10-letter name and
#  sequence content)
#The corresponding tree also needs the same 10-letter names


def _correct_phylip_file(filename):
    """
    Take a sequence file formatted in phylip sequential format, as
    written by Biopython's Bio.SeqIO.convert, and correct it so it works with PAML.
    i.e. after formatting should look like:
    114158638_  MVASGILLVALLTCLTVMVLM...
    511984835_  MLASGLLLVALLACLTIMVLM...
    585192830_  MLASGLLLVALLTCLTTMVLM...
    where each sequence takes a single line and there are 2 spaces between the 10-letter name and sequence content)

    :param filename: sequence file in phylip sequential format
    :return: None
    """
    out = []
    with open(filename, 'r') as fh:
        lines = [line.strip() for line in fh.readlines()]
        out.append(lines[0]) # header line we keep
        for line in lines[1:]:
            first = line[:10]
            second = line[10:]
            corrected_line = first + "  " + second
            out.append(corrected_line)
    with open(filename, 'w') as fh:
        print('\n'.join(out), sep= "", end="", file=fh)


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


def _partition(lst, n):
    """
    Take a list, randomly shuffle and partition into n groups.
    :return : nested list
    """
    random.shuffle(lst)
    division = len(lst) / n # number of sequences per group
    out = []
    for i in range(n):
        out.append(lst[round(division * i) : round(division * (i + 1))])
    return out


def _clean_aligned_seqs(seqs):
    """
    Return sequences with gappy columns removed
    """
    gappy_columns = set()
    align_len = len(seqs[0])
    for seq in seqs:
        for i in range(align_len):
            if i not in gappy_columns and seq[i] == "-":
                gappy_columns.add(i)
    new_seqs = []
    for seq in seqs:
        content = "".join([seq[i] for i in range(align_len) if i not in gappy_columns])
        new_seqs.append(Sequence(sequence = content, alphabet = seq.alphabet, name = seq.name))
    return new_seqs, sorted(list(gappy_columns))

def _process_arguments(myparser, myargs):
    # Load base sequences and tree
    base_seqs = readFastaFile(myargs.sequences, Protein_Alphabet, gappy=True)
    # FIXME: truncation can lead to identical names. Fix this by assigning a random alphanumeric name and/or
    # FIXME : storing the matching original name in a textfile. I don't think we need correct names for analysis.

    # Assign random names to avoid duplicate names
    new_names = dict()
    for seq in base_seqs:
        new_names[seq.name] = _generateRandomStringID(10,alpha="ABCDEFGHIJKLMOPQRSTUVWXYZ")
        seq.name = new_names[seq.name]
    base_tree = Phylo.read(myargs.tree, "newick")
    for node in base_tree.get_terminals():
        node.name = new_names[node.name]
        node.confidence = None
    for node in base_tree.get_nonterminals():
        node.name = ''
        node.confidence = None

    # Make results directory for current group size
    group_size_path = Path(os.path.join(myargs.output, myargs.id, "group_size_" + str(myargs.groups)))
    group_size_path.mkdir(parents = True, exist_ok=True)

    cleaned_seqs, gappy_columns = _clean_aligned_seqs(base_seqs) # **** Comment out to not clean base seqs


    # For each repeat
    for repeat in range(myargs.repeats):
        # Make folder for repeat
        repeat_path = Path(os.path.join(str(group_size_path.absolute()), "repeat_" + str(repeat +1)))
        repeat_path.mkdir(parents = True, exist_ok=True)

        # Subset sequences
        # seqs_copy = copy.deepcopy(base_seqs) # **** replace line below with this to clean each subset
        seqs_copy = copy.deepcopy(cleaned_seqs)
        seq_subsets = _partition(seqs_copy, myargs.groups)

        # Generate the input sequences and tree for each subset. Write to file.
        group_cnt = 1
        for sb in seq_subsets:
            # Identify and purge gappy columns from the subset of aligned sequences
            # cleaned_seqs, gappy_columns = _clean_aligned_seqs(sb)
            group_path = Path(os.path.join(str(repeat_path.absolute()), "group_" + str(group_cnt)))
            group_path.mkdir(parents=True, exist_ok=True)
            # Copy the original tree
            copy_tree = copy.deepcopy(base_tree)
            sb_names = [s.name for s in sb] # names of sequences in subset
            to_prune = [s.name for s in seqs_copy if s.name not in sb_names] # Nodes to prune
            # Prune and collapse nodes. Maintains bifurication. If collapsing occurs, adds branch lengths together.
            for tp in to_prune:
                copy_tree.prune(tp)

            # Write subset newick tree
            with open(os.path.join(str(group_path.absolute()), "{}_gs_{}_r_{}_g_{}.newick".format(myargs.id, myargs.groups,
                                                                               repeat+1,
                                                                               group_cnt)),'w') as fh:
                # print(copy_tree.format("newick"), file = fh)
                # Biopython attaches a branch length of 0.0000 to the end of the newick string.
                # This needs to be removed so GRASP can read the newick string
                print(':'.join(copy_tree.format('newick').split(":")[:-1]) + ";", file=fh)

            # Write subset sequences as fasta file
            subset_seq_filename = os.path.join(str(group_path.absolute()), "{}_gs_{}_r_{}_g_{}.fasta".format(myargs.id,
                                                                                             myargs.groups,
                                                                                             repeat+1,
                                                                                             group_cnt))
            # writeFastaFile(subset_seq_filename, cleaned_seqs) # **** replace line below with this to used cleaned subset
            writeFastaFile(subset_seq_filename, sb)

            # Convert fasta file to phylip-sequential format
            subset_seq_corrected_filename = os.path.join(str(group_path.absolute()), "{}_gs_{}_r_{}_g_{}.phylip".format(myargs.id,
                                                                                                        myargs.groups,
                                                                                                        repeat+1,
                                                                                                        group_cnt))
            ## First convert the exiting sequence file to phylip format
            SeqIO.convert(subset_seq_filename, "fasta", subset_seq_corrected_filename, "phylip-sequential")
            ## And then correct the file so there are spaces after each sequence name
            _correct_phylip_file(subset_seq_corrected_filename)

            # Now write another copy of the tree but with matching 10-letter names
            terminals = copy_tree.get_terminals()
            for node in terminals:
                node.name = node.name[:10]
            with open(os.path.join(str(group_path.absolute()), "{}_gs_{}_r_{}_g_{}_phylip.newick".format(myargs.id,myargs.groups,
                                                                               repeat+1,
                                                                               group_cnt)),'w') as fh:
                print(copy_tree.format("newick"), file = fh)

            # Write log file
            with open(os.path.join(str(group_path.absolute()), "{}_gs_{}_r_{}_g_{}_log.txt".format(myargs.id,myargs.groups,
                                                                               repeat+1,
                                                                               group_cnt)),'w') as fh:

                print("gappy columns: " + " ".join(map(str, gappy_columns)), file = fh)
                print("Nseqs: " + str(len(sb)), file = fh)
                print("Alignlen_with_gaps: " + str(len(base_seqs[0])), file=fh)
                print("Alignlen_without_gaps: " + str(len(cleaned_seqs[0])), file=fh)


            group_cnt += 1




if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Generate input data and directory structures for evaluating ASR tools')
    parser.add_argument('-s', '--sequences', help='Input alignment file formatted as fasta', required=True)
    parser.add_argument('-t', '--tree', help='Input tree file formatted as newick', required=True)
    parser.add_argument('-o', '--output', help='Output directory', required=True)
    parser.add_argument('-r', '--repeats', help='Number of repeats to perform', type=int, default=1)
    parser.add_argument('-g', '--groups', help='Number of groups to randomly separate sequences into',
                        type=int, default=10)
    parser.add_argument('-id', '--id', help='ID tag for data', required=True)


    myargs = ["-s","../../Data/CYP2/CYP2_input.fa",
            "-t", "../../Data/CYP2/CYP2_input.nwk",
            "-o", "../../Data/test_CYP/test_out/",
            "-r", "1",
            "-g", "40",
            "-id", "CYP"]


    myargs = ["-s","../../Data/KARI/KARI_EC_mafft3.fasta",
            "-t", "../../Data/KARI/KARI_EC_mafft3.nwk",
            "-o", "../../Data/test_KARI/test_out/",
            "-r", "1",
            "-g", "40",
            "-id", "KARI"]

    # args = parser.parse_args(myargs)
    args = parser.parse_args()
    _process_arguments(parser, args)

