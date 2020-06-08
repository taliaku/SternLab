#! /usr/local/python_anaconda/bin/python3.4

from optparse import OptionParser
import glob
from file_utilities import check_dirname, check_filename
import pandas as pd
import re
from scipy.stats import chi2_contingency, fisher_exact
from Bio.Seq import Seq
from Bio.Seq import Seq
from os import path
from collections import Counter
from itertools import permutations
from general_utilities import merge_two_dicts
from fastq_utilities import get_sequence_and_ancestry_data, get_position_to_remove



def main():
    parser = OptionParser("usage: %prog [options]")
    parser.add_option("-d", "--dir", dest="dir", help="directory with fastml output files")
    parser.add_option("-o", "--output", dest="output", help="output file perfix")
    parser.add_option("-m", "--marginal", dest="marginal", help="prob_marginal_file")


    (options, args) = parser.parse_args()
    dir = options.dir
    if dir != None:
        dir = check_dirname(dir) + "/"
    output = options.output
    marginal = options.marginal
    marginal = check_filename(marginal)
    basenames = [marginal.split(".prob.marginal.txt")[0]]
    if output == None and marginal != None:
        output   = "%s_fastml_analysis.csv" % basenames[0]
    elif output == None:
        output = "fastml_anlysis_new.csv"


    print (output)


    print(dir)
    if marginal == None:
        files = glob.glob(dir  + "*prob.marginal.txt")
        print(files)
        basenames = [f.split(".prob.marginal.txt")[0] for f in files]
    df = pd.DataFrame(columns=["Basename", "Mutation", "Branch", "Context", "Mutation_type",
                               "Codon_position", "APOBEC_context_GA","APOBEC_context_CT", "type_of_pos"])

    for basename in basenames:
        print(basename)
        prob_marginal = basename + ".prob.marginal.txt"
        seq_marginal = basename + ".seq.marginal.txt"
        tree_ancestor = basename + ".tree.ancestor.txt"

        basename = basename.split("/")[-1]

        positions_to_remove = get_position_to_remove(prob_marginal)
        ancestor_info, seqs = get_sequence_and_ancestry_data(tree_ancestor, seq_marginal)
        df = go_over_positions(ancestor_info, seqs, positions_to_remove, basename, df)

    print(df)

    if(dir==None):
        df.to_csv(output, index=False)
    else:
        df.to_csv(dir + "/" + output, index=False)


def go_over_positions_old(ancestor_info, seqs, positions_to_remove, basename, df, mutations_to_check=["GA"]):
    same_positions = 0
    bad_positions = 0
    positions_with_one_difference = 0
    positions_with_more_than_one_difference = 0

    for son in ancestor_info:
        father = ancestor_info[son]
        if father == "root!":
            continue
        son_seq = seqs[son]
        father_seq = seqs[father]
        pattern = re.compile("^N\d*") #checks if the node is an internal node
        if pattern.findall(son) == []:
            branch = "external"
        else:
            branch = "internal"


        for i in range(0, len(father_seq) - 1):
            if i in positions_to_remove:
                bad_positions += 1
                type_of_pos = "diverse"
            else:
                type_of_pos = "not-diverse"
            context = father_seq[i-1:i+2]

            father_nuc = father_seq[i]
            son_nuc = son_seq[i]

            if son_nuc == father_nuc:
                same_positions += 1

                df = df.append({"Basename": basename, "Mutation": son_nuc+son_nuc, "Branch": branch, "Context": context,
                                "Mutation_type": None, "Codon_position": (i % 3) + 1,
                                "APOBEC_context_GA": None,
                                "APOBEC_context_CT": None, "type_of_pos": type_of_pos}, ignore_index=True)
                continue
            if son_nuc not in ["A", "C", "T", "G"]:
                continue
            elif "-" in son_nuc:
                continue
            if i == 0 or i == len(father_seq) - 2:
                continue
            codon_position = i % 3
            if codon_position == 0:
                father_codon = father_seq[i:i+3]
                son_codon = son_seq[i:i+3]
            elif codon_position == 1:
                father_codon = father_seq[i-1:i+2]
                son_codon = son_seq[i-1:i+2]
            else:
                father_codon = father_seq[i-2:i+1]
                son_codon = son_seq[i-2:i+1]
            if "-" in son_codon:
                continue
            changes = [j for j in range(len(father_codon)) if father_codon[j] != son_codon[j]]
            number_of_changes = len(changes)
            mutation = father_nuc + son_nuc
            father_aa = str(Seq(father_codon).translate())
            son_aa = str(Seq(son_codon).translate())
            if father_aa == "*":
                continue
            if number_of_changes == 1:
                positions_with_one_difference += 1

                if father_aa == son_aa:
                    mutation_type = "synonymous"
                    if i%3 == 1:
                        print(i, father_codon, son_codon)
                else:
                    mutation_type = "non-synonymous"

            elif number_of_changes >= 2:
                positions_with_more_than_one_difference +=1
                if i%3 == 0:
                    mutation_type = "non-synonymous"
                elif i%3 == 1:
                    mutation_type = "non-synonymous"
                elif i%3 == 2:
                    mutation_type = "synonymous"

            #print(i, context)
            if  context[2] in ["G", "A"]:
                APOBEC_context_GA = True
            else:

                APOBEC_context_GA =  False
            if context[0] in ["C", "T"]:
                APOBEC_context_CT = True
            else:
                APOBEC_context_CT = False




            df = df.append({"Basename": basename, "Mutation": mutation, "Branch": branch, "Context": context,
                                "Mutation_type": mutation_type, "Codon_position": (i % 3) + 1,
                            "APOBEC_context_GA":APOBEC_context_GA,
                            "APOBEC_context_CT":APOBEC_context_CT, "type_of_pos":type_of_pos}, ignore_index=True)

    return(df)


def go_over_positions(ancestor_info, seqs, positions_to_remove, basename, df, mutations_to_check=["GA"]):
    same_positions = 0
    bad_positions = 0

    df=pd.DataFrame()
    options = ["".join(i) for i in list(permutations(["A", "C", "G", "T"], 2))] + ["AA", "GG", "CC", "TT"]
    external_counts = {i:0 for i in options}
    internal_counts = {i:0 for i in options}
    overall_counts = {i:0 for i in options}
    print(external_counts, internal_counts, overall_counts)

    for son in (ancestor_info):
        father = ancestor_info[son]
        if father == "root!":
            continue
        son_seq = seqs[son]
        father_seq = seqs[father]
        pattern = re.compile("^N\d*") #checks if the node is an internal node
        if pattern.findall(son) == []:
            branch = "external"
        else:
            branch = "internal"

        non_diverse_seq_father = ""
        non_diverse_seq_son = ""

        for i in range(len(father_seq)):
            if i not in positions_to_remove:
                non_diverse_seq_father += father_seq[i]
                non_diverse_seq_son += son_seq[i]

        diffs = [non_diverse_seq_father[i]+non_diverse_seq_son[i] for i in range(len(non_diverse_seq_father))]
        counter = Counter(diffs)

        for i in options:
            overall_counts[i] += counter[i]
            if branch == "external":
                external_counts[i] += counter[i]
            else:
                internal_counts[i] += counter[i]

    print(external_counts, internal_counts, overall_counts)

    for i in options:
        sum_all = 0
        sum_all_external = 0
        sum_all_internal = 0
        for k in ["A", "C", "G", "T"]:
            sum_all_external +=  external_counts[i[0]+k]
            sum_all_internal +=  internal_counts[i[0]+k]
            sum_all +=  overall_counts[i[0]+k]

        external_counts[i + "_ratio"] = external_counts[i] / sum_all_external
        internal_counts[i + "_ratio"] = internal_counts[i] / sum_all_internal
        overall_counts[i + "_ratio"] = overall_counts[i] / sum_all


    print(external_counts, internal_counts, overall_counts)
    external_counts["basename"] = basename
    external_counts["branch"] = "external"
    internal_counts["basename"] = basename
    internal_counts["branch"] = "internal"
    overall_counts["basename"] = basename
    overall_counts["branch"] ="all"

    df = df.append(external_counts, ignore_index=True)
    df = df.append(internal_counts, ignore_index=True)
    df = df.append(overall_counts, ignore_index=True)

    return(df)


if __name__ == "__main__":
    main()
