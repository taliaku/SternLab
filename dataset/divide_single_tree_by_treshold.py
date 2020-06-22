#! /powerapps/share/python-anaconda-3.6/bin/python

import argparse
from file_utilities import check_dirname, check_filename
import glob

from seqFileTools import count_seqs_in_fasta, convert_fasta_to_phylip
from seqFileAnalyzer import get_consensus_percentage
import pandas as pd
from PAML_utilities import write_ctl_file
from Bio import SeqIO
from Bio import Phylo
import re
import os
from phylogenetic_utilities import branch_over_threshold
import time

def main(args):
    tree = check_filename(args.tree)
    fasta = check_filename(args.fasta)
    treshold = args.treshold
    alignment_type = args.alignment_type
    tree_type = args.tree_type
    count = 1

    file_prefix = fasta.split(".fasta")[0] + "_divided"
    groups, added = divide_tree_by_threshold_wrapper(tree, treshold)

    print(f"found {len(groups)} groups")
    for i in range(len(groups)):
        print(f"group number {i+1}: {len(groups[i])}")

    fasta_file_text = open(fasta, "rb").read()
    num_of_seq = count_seqs_in_fasta(fasta)

    if len(groups) == 1 and len(groups[0]) == num_of_seq: #all seqs are in the same group
        print("all seqs are considered to be in the same group")
    else:
        print("running on groups:")
        for group_members in groups:
            print(f"group number {count}")
            if len(group_members) >= 8:
                count = run_on_group(group_members, fasta, file_prefix, count, treshold, alignment_type, tree_type)
            else:
                print("the group is too small (bellow 8)")


def divide_tree_by_threshold_wrapper(treefile, threshold):
    # a wrapper to the divide_tree_by_threshold function. just opens the tree file.
    tree = Phylo.read(treefile, "newick")
    groups, added = divide_tree_by_treshold(tree, threshold, [], [])
    return groups, added


def divide_tree_by_treshold(tree, threshold, groups=[], added=[]):
    # runs on all branches of the tree and divides it if there are branches over the threshold
    # runs on a recursive manner
    clades = list(tree.find_clades())
    pruned = []
    had_branches_over_threshold = False
    for clade in clades:
        if clade.branch_length==None:
            continue
        if clade.branch_length >= threshold:
            inside_tree = Phylo.BaseTree.Tree.from_clade(clade)
            inside_tree.root.branch_length = 0

            if len(clade.get_terminals()) == 1:
                try:
                    tree.prune(clade)
                except:
                    pass
                continue

            for c in clade.get_terminals():
                if not c.name in pruned:
                    tree.prune(c)
                    pruned.append(c.name)

            groups, added = divide_tree_by_treshold(tree, threshold, groups, added)
            groups, added = divide_tree_by_treshold(inside_tree, threshold, groups, added)
            had_branches_over_threshold = True

    if had_branches_over_threshold == False:
        terminals = [str(terminal.name) for terminal in tree.get_terminals()]
        if len(terminals) > 0:
            new_terminals = []
            for terminal in terminals:
                if terminal not in added:
                    new_terminals.append(terminal)
            if new_terminals != []:
                groups.append(new_terminals)

            for terminal in new_terminals:
                added.append(terminal)
    return groups, added


def run_on_group(group_members, fasta_file, file_prefix, count, treshold, alignment_type="mafft", tree_type="njTree"):
    print(f"group contenis {len(group_members)} members")
    print(group_members)
    with open(fasta_file, "r") as handle:
        fasta_file_text = handle.read()
    new_fasta_text = ""
    for id in group_members:
        id_for_search = id.replace("|", "\|")
        id_for_search = id_for_search.replace(".", "\.")
        pattern = re.compile("\>%s[^>]*" % id_for_search, re.DOTALL)
        seqs = re.findall(pattern, fasta_file_text)
        found = False
        for s in seqs:
            #print(s.split("\n"))
            if s.split("\n")[0].split(">")[1] == id:
                found = True
                break
        if found:
            new_fasta_text += s
        else:
            raise(TypeError("didn't find seq!!!"))

    new_fasta_file = file_prefix + f'_{count}.fasta'

    with open(new_fasta_file, "w") as handle:
        handle.write(new_fasta_text)
    if alignment_type=="mafft":
        alignment_file = file_prefix + f"_{count}.aln"
        cmd = f"/sternadi/home/volume1/taliakustin/software/mafft-7.300-with-extensions/scripts/mafft --quiet  --retree 1 --maxiterate 0 {new_fasta_file} > {alignment_file}"
    elif alignment_type=="prank":
        alignment_file_for_cmd = file_prefix + f"_{count}.aln"
        alignment_file = file_prefix + f"_{count}.aln.best.fas"
        cmd = f"/powerapps/share/bin/prank -d={new_fasta_file}  -F -o={alignment_file_for_cmd}"
    if not os.path.isfile(alignment_file):
        print(f"running {alignment_type} alignment")
        print(cmd)
        start_time = time.time()
        os.system(cmd)
        end_time = time.time()
        print(f"{alignment_type} took {end_time - start_time} seconds")
    else:
        print(f"alignment file {alignment_file} exist")
    if tree_type == "njTree":
        tree_file = file_prefix + f"_{count}.tree"
        cmd = f"/sternadi/home/volume1/shared/tools/phylogenyCode/programs/treeUtil/njTreeJCdist -i {alignment_file} -o {tree_file} -an"
    elif tree_type == "phyml":
        phylip_file = convert_fasta_to_phylip(alignment_file)
        tree_file = phylip_file + "_phyml_tree.txt"
        cmd = f"/sternadi/home/volume1/shared/tools/PhyML/PhyML_3.0_linux64 -i {phylip_file} -b 0"
    if not os.path.isfile(tree_file):
        print(f"running {tree_type} alignment")
        print(cmd)
        start_time = time.time()
        os.system(cmd)
        end_time = time.time()
        print(f"{tree_type} took {end_time - start_time} seconds")
    else:
        print(f"tree file {tree_file} exist")
    count += 1
    if not branch_over_threshold(tree_file, treshold):
        return count
    groups, added = divide_tree_by_threshold_wrapper(tree_file, treshold)
    if len(groups) == 1:
        new_group_members = groups[0]
        new_group_members.sort()
        if group_members == new_group_members:
            return count

    print(f"    found {len(groups)} groups")
    for group_members in groups:
        if len(group_members) >= 8:
            count = run_on_group(group_members, new_fasta_file, file_prefix, count, treshold, alignment_type, tree_type)
        else:
            print("     the group is too small (bellow 8)")
    return count




if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-t", "--tree", type=str,
                        help="tree", required=True)
    parser.add_argument("--treshold", type=float,
                        help="branch treshold to cut files", required=True, default=0.5)
    parser.add_argument("-f", "--fasta", type=str,
                        help="fasta file", required=True)
    parser.add_argument("-i", "--index", type=int, help="index", default=1)
    parser.add_argument("--alignment_type", type=str, help="mafft or prank (default:mafft)", default="mafft")
    parser.add_argument("--tree_type", type=str, help="njTree or phyml (default:njTree)", default="njTree")

    args = parser.parse_args()
    main(args)
