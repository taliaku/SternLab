#! /powerapps/share/python-anaconda-3.6/bin/python

import os
import glob
from file_utilities import check_dirname, check_filename
import pandas as pd
import re
from PAML_utilities import *
from phyVirus.get_baltimore import get_baltimore_classifiaction
from scipy.stats.distributions import chi2
from statsmodels.stats.multitest import multipletests
import numpy as np
from io import StringIO
from Bio import AlignIO, SeqIO
from Bio.Seq import Seq
from Bio import Phylo
import tqdm
#from phyVirus.phyVirus_utilities import get_basename
import statsmodels.api as sm

def get_mutation_df(basename, cutoff = 0.75, onlyDiffPoss = False, infoDictionary={}, phyVirus=True, rooted=True, overwrite = False):
    phyVirus_analysis_dir = "/sternadi/home/volume3/taliakustin/phyVirus_analysis/"
    fastml_dir = "fastml/"
    if not rooted:
        fastml_dir = "fastml_unrooted_tree/"
    output = phyVirus_analysis_dir + fastml_dir + basename + "mutation.map"

    if os.path.isfile(output) and not overwrite:
        print("output file already exists - returns. %s" % output)
        return output
    prob_marginal = glob.glob(phyVirus_analysis_dir + fastml_dir + basename + "*prob.marginal.txt")[0]
    seq_marginal = glob.glob(phyVirus_analysis_dir + fastml_dir +basename + "*seq.marginal.txt")[0]
    tree_ancestor = glob.glob(phyVirus_analysis_dir + fastml_dir +basename + "*tree.ancestor.txt")[0]

    #get all positions to remove from analysis because their marginal probability is lower than the cutoff
    positions_to_remove = {}
    pm = open(prob_marginal, "r").read()
    poss = pm.split("\n\nmarginal probabilities at ")[1:]
    p = re.compile(": (.+): p\([GATC]\)=([01].\d*)") #gets probabilities under 0.7
    for pos in poss:
        r = p.findall(pos)
        for i in r:
            if float(i[1]) < cutoff:
                pos_num = int(pos.split("\n")[0].split("position: ")[1]) - 1
                node = i[0]
                if node in positions_to_remove:
                    positions_to_remove[node].append(pos_num)
                else:
                    positions_to_remove[node] = [pos_num]

    alignment_len = int(pos_num) - 1


    #create ancestor info - for each son - who is the father
    pattern = re.compile(r'\s+')
    ancestors =  open(tree_ancestor, "r").read().split("\n")[2:-2]
    ancestor_info = {}
    for line in ancestors:
        line = re.sub(pattern, '$', line)
        line = line.split("$")
        son = line[0]
        father = line[1]
        ancestor_info[son] = father

    # create seq dictionary - name and sequance
    seqs = open(seq_marginal, "r").read()
    seqs = seqs.split("\n>")[1:]
    seqs = {seq.split("\n")[0]: seq.split("\n")[1] for seq in seqs}

    res = pd.DataFrame()

    for son in ancestor_info:
        father = ancestor_info[son]
        if father == "root!":
            continue
        son_seq = list(seqs[son])
        father_seq = list(seqs[father])
        pattern = re.compile("^N\d*")  # checks if the node is an internal node
        if pattern.findall(son) == []:
            node_type = "external"
        else:
            node_type = "internal"

        if onlyDiffPoss:
            posToIterate = [i for i in range(len(son_seq)) if son_seq[i] != father_seq[i]]
        else:
            posToIterate = range(len(son_seq))

        to_remove = []
        if son in positions_to_remove:
            for i in positions_to_remove[son]:
                to_remove.append(i)
        if father in positions_to_remove:
            for i in positions_to_remove[father]:
                to_remove.append(i)

        final_to_remove = []
        [final_to_remove.append(int(x)) for x in to_remove if x not in final_to_remove]
        for p in sorted(final_to_remove, reverse=True):
            del son_seq[p]
            del father_seq[p]


        noGap_son_seq = []
        noGap_father_seq = []
        for i in range(len(son_seq)):
            if son_seq[i] != "-" and father_seq[i] != "-":
                noGap_son_seq.append(son_seq[i])
                noGap_father_seq.append(father_seq[i])

        if len(noGap_son_seq) != len(noGap_father_seq):
            raise Error("not the same length")
        muts = ["%s%s" % (noGap_father_seq[i], noGap_son_seq[i]) for i in range(len(noGap_son_seq))]
        temp_res = pd.DataFrame(index=range(len(muts)), columns=["father_node_name", "son_node_name", "letter1", "letter2", "mutation", "branch"])
        temp_res.mutation = muts
        temp_res.father_node_name = father
        temp_res.son_node_name = son
        temp_res.branch = node_type
        temp_res.letter1 = noGap_father_seq
        temp_res.letter2 = noGap_son_seq


        res = res.append(temp_res, ignore_index=True)

    for key in infoDictionary.keys():
        res[key] = infoDictionary[key]

    res.to_csv(output)
    return(output)


def summmerize_mutations(file, ratios_output=None, summary_output=None, overwrite=False, cutoff="no_cutoff"):
    print(file)
    file = check_filename(file)
    nucs = ["A", "G", "T", "C"]
    if ratios_output == None:
        ratios_output = file.split(".mutation.map")[0] + ".to_nucs_ratio"
    else:
        ratios_output = check_filename(ratios_output, Truefile=False)
    if summary_output == None:
        summary_output = file.split(".mutation.map")[0] + ".mutation_summary"
    else:
        summary_output = check_filename(summary_output, Truefile=False)
    subs_ratios_output =  file.split(".mutation.map")[0] + ".to_substitutions_ratio"

    if os.path.isfile(ratios_output) and not overwrite:
        print("output file %s exists" % ratios_output)
        print(summary_output)
        return(summary_output)
    base = file.split("/")[-1].split(".mutation.map")[0].split(".joint")[0]
    family = base.split("_")[0]
    baltimore = get_baltimore_classifiaction(family)

    df = pd.read_csv(file, index_col=False)

    if df.shape[0] == 0:
        print("%s - has not mutation mapping" % file)
        return
    df["mut"] = df.apply(mut_or_not, axis=1)
    df["count"] = 1

    summary = df.groupby(["letter1", "letter2", str(cutoff), "mut"], as_index=False).agg({"count": "count"})
    summary["branch"] = summary[str(cutoff)]
    summary["base"] = base
    summary["family"] = family
    summary["baltimore"] = baltimore
    summary["cutoff"] = cutoff
    summary = summary.drop([str(cutoff)], axis=1)
    mutation_count = sum(summary.loc[(summary.branch == "external") & (summary.mut == "mut"), "count"])
    summary.to_csv(summary_output)
    ratios_to_nucs = pd.DataFrame()

    for branch in ["external", "internal"]:
        for nuc in nucs:
            if (sum(summary.loc[(summary["mut"]=="mut") & (summary["branch"] == branch), "count"].values)) != 0:
                #rate of mutation to this nuc from all mutations to all nucs
                to_nuc_rate = sum(summary.loc[(summary["letter2"] == nuc) & (summary["letter1"] != nuc) & (
                summary["branch"] == branch), "count"].values) \
                              / sum(summary.loc[(summary["mut"]=="mut") & (summary["branch"] == branch), "count"].values)
            else:
                to_nuc_rate = 0
            #num of position that can go to nuc from all positions available to mutate
            if sum(summary.loc[(summary["letter1"] != nuc) & (summary["branch"] == branch), "count"].values) == 0:
                potential = 0
            else:
                potential = sum(summary.loc[(summary["letter1"] != nuc) & (summary["branch"] == branch), "count"].values)\
                        / sum(summary.loc[(summary["branch"] == branch), "count"].values)
            if potential == 0:
                ratios_to_nucs = ratios_to_nucs.append(
                    {"basename": base, "family": family, "baltimore": baltimore, "branch": branch,
                     "nuc": nuc, "rate": 0, "cutoff": cutoff}, ignore_index=True)
            else:
                ratios_to_nucs = ratios_to_nucs.append(
                {"basename": base, "family": family, "baltimore": baltimore, "branch": branch,
                 "nuc": nuc, "rate": to_nuc_rate/potential, "cutoff":cutoff}, ignore_index=True)
    ratios_to_nucs["mutation_count"] = mutation_count
    ratios_to_nucs.to_csv(ratios_output)

    ratio_to_subs = pd.DataFrame()
    for branch in ["external", "internal"]:
        for first_nuc in nucs:
            for second_nuc in nucs:
                if first_nuc == second_nuc:
                    continue
                if (sum(summary.loc[(summary["mut"]=="mut") & (summary["branch"] == branch), "count"].values)) != 0:
                    #rate of mutation to this nuc from all mutations to all nucs
                    to_nuc_rate = sum(summary.loc[(summary["letter2"] == second_nuc) & (summary["letter1"] == first_nuc) & (
                    summary["branch"] == branch), "count"].values) \
                                  / sum(summary.loc[(summary["mut"]=="mut") & (summary["branch"] == branch), "count"].values)
                else:
                    to_nuc_rate = 0
                #num of position that can go to nuc from all positions available to mutate
                if sum(summary.loc[(summary["letter1"] != nuc) & (summary["branch"] == branch), "count"].values) == 0:
                    potential = 0
                else:
                    potential = sum(summary.loc[(summary["letter1"] == first_nuc) & (summary["branch"] == branch), "count"].values)\
                            / sum(summary.loc[(summary["branch"] == branch), "count"].values)

                if potential == 0:
                    ratio_to_subs = ratio_to_subs.append(
                        {"basename": base, "family": family, "baltimore": baltimore, "branch": branch,
                         "sub": first_nuc + second_nuc, "rate": 0, "cutoff": cutoff},
                        ignore_index=True)
                else:
                    ratio_to_subs = ratio_to_subs.append(
                        {"basename": base, "family": family, "baltimore": baltimore, "branch": branch,
                         "sub": first_nuc+second_nuc, "rate": to_nuc_rate/potential,  "cutoff":cutoff}, ignore_index=True)
    ratio_to_subs["mutation_count"] = mutation_count
    ratio_to_subs.to_csv(subs_ratios_output)

    return(summary_output)


def mut_or_not(row):
    if row["letter1"] == row["letter2"]:
        return "no_mut"
    return "mut"


def get_letter1(row):
    return row["mutation"][0]

def get_letter2(row):
    return row["mutation"][1]

def get_mutation_joint_df(basename, onlyDiffPoss = False, infoDictionary={}, rooted=True, overwrite = False):
    phyVirus_analysis_dir = "/sternadi/home/volume3/taliakustin/phyVirus_analysis/"
    fastml_dir = "fastml_midpoint_tree/"
    if not rooted:
        fastml_dir = "fastml_unrooted_tree/"
    output = phyVirus_analysis_dir + fastml_dir + basename + ".joint.mutation.map"
    branch_file = glob.glob(phyVirus_analysis_dir + "/branch_length_info/" + basename + ".*")[0]

    if os.path.isfile(output) and not overwrite:
        print("output file already exists - returns. %s" % output)
        return output
    print(phyVirus_analysis_dir + fastml_dir +basename + ".seq.joint.txt")
    seq_joint = glob.glob(phyVirus_analysis_dir + fastml_dir +basename + ".seq.joint.txt")[0]
    tree_ancestor = glob.glob(phyVirus_analysis_dir + fastml_dir +basename + ".tree.ancestor.txt")[0]


    #create ancestor info - for each son - who is the father
    pattern = re.compile(r'\s+')
    ancestors =  open(tree_ancestor, "r").read().split("\n")[2:-2]
    ancestor_info = {}
    for line in ancestors:
        line = re.sub(pattern, '$', line)
        line = line.split("$")
        son = line[0]
        father = line[1]
        ancestor_info[son] = father

    # create seq dictionary - name and sequance
    seqs = open(seq_joint, "r").read()
    seqs = seqs.split("\n>")[1:]
    seqs = {seq.split("\n")[0]: seq.split("\n")[1] for seq in seqs}
    res = pd.DataFrame()
    for son in ancestor_info:
        father = ancestor_info[son]
        if father == "root!":
            continue
        son_seq = list(seqs[son])
        father_seq = list(seqs[father])
        pattern = re.compile("^N\d*")  # checks if the node is an internal node
        if pattern.findall(son) == []:
            node_type = "external"
        else:
            node_type = "internal"

        if onlyDiffPoss:
            posToIterate = [i for i in range(len(son_seq)) if son_seq[i] != father_seq[i]]
        else:
            posToIterate = range(len(son_seq))

        noGap_son_seq = []
        noGap_father_seq = []
        for i in range(len(son_seq)):
            if son_seq[i] != "-" and father_seq[i] != "-":
                noGap_son_seq.append(son_seq[i])
                noGap_father_seq.append(father_seq[i])

        if len(noGap_son_seq) != len(noGap_father_seq):
            raise Error("not the same length")
        muts = ["%s%s" % (noGap_father_seq[i], noGap_son_seq[i]) for i in range(len(noGap_son_seq))]
        temp_res = pd.DataFrame(index=range(len(muts)), columns=["father_node_name", "son_node_name", "letter1", "letter2", "mutation", "no_cutoff"])
        temp_res.mutation = muts
        temp_res.father_node_name = father
        temp_res.son_node_name = son
        temp_res.no_cutoff = node_type
        temp_res.letter1 = noGap_father_seq
        temp_res.letter2 = noGap_son_seq


        res = res.append(temp_res, ignore_index=True)



    for key in infoDictionary.keys():
        res[key] = infoDictionary[key]



    res["node_name"] = res["son_node_name"]
    branch = pd.read_csv(branch_file)
    res_merged = res.merge(branch, how='outer', on=['node_name', "base", "family", "baltimore"])

    for c in [0.05, 0.04, 0.03, 0.02, 0.01, 0.005, 0.001]:
        res_merged[str(c)] = "internal"
        res_merged.loc[res_merged.branch_len < c, str(c)] = "external"

    res_merged.to_csv(output)
    return (summmerize_mutations)

def different_cutoffs(mut_file, branch_file, cutoffs=[0.05, 0.01, 0.005, 0.001]):
    print(mut_file)
    base = mut_file.split("/")[-1].split(".mutation.map")[0].split(".joint.mutation.map")[0]
    family = base.split("_")[0]
    baltimore = get_baltimore_classifiaction(family)

    mut = pd.read_csv(mut_file)
    mut = mut.rename(columns={"son_node_name": "node_name"})
    branch = pd.read_csv(branch_file)
    summary_output = mut_file.split(".joint.mutation.map")[0].split(".mutation.map")[0] + ".mutation_summary_cutoffs"
    ratios_output = mut_file.split(".joint.mutation.map")[0].split(".mutation.map")[0] + ".ratios_to_nucs_cutoffs"
    nucs = ["A", "G", "T", "C"]
    mut_merged = mut.merge(branch, how='outer', on=['node_name', "base", "family", "baltimore"])
    colnames = ["branch"]
    for c in cutoffs:
        col_name = "%s_cutoff" % str(c)
        colnames.append(col_name)
        mut_merged[col_name] = "internal"
        mut_merged.loc[mut_merged.branch_len  < c, col_name] = "external"
    mut_merged["mut"] = mut_merged.apply(mut_or_not, axis=1)
    mut_merged["count"] = 1
    all_summaries = pd.DataFrame()
    all_ratios_to_nucs = pd.DataFrame()
    for col_name in colnames:
        ratios_to_nucs = pd.DataFrame()
        summary = mut_merged.groupby(["letter1", "letter2", col_name, "mut", "base", "family", "baltimore"], as_index=False).agg({"count": "count"})    
        cutoff = col_name
        summary = summary.rename(columns={col_name:"branch"})
        if cutoff == "branch":
            cutoff = "no_cutoff"
        summary["cutoff"] = cutoff
        for branch in ["external", "internal"]:
            for nuc in nucs:
                if (sum(summary.loc[(summary["mut"] == "mut") & (summary["branch"] == branch), "count"].values)) != 0:
                    # rate of mutation to this nuc from all mutations to all nucs
                    to_nuc_rate = sum(summary.loc[(summary["letter2"] == nuc) & (summary["letter1"] != nuc) & (
                            summary["branch"] == branch), "count"].values) \
                                  / sum(
                        summary.loc[(summary["mut"] == "mut") & (summary["branch"] == branch), "count"].values)
                else:
                    to_nuc_rate = 0

                # num of position that can go to nuc from all positions available to mutate
                potential = sum(
                    summary.loc[(summary["letter1"] != nuc) & (summary["branch"] == branch), "count"].values) \
                            / sum(summary.loc[(summary["branch"] == branch), "count"].values)

                ratios_to_nucs = ratios_to_nucs.append(
                    {"base": base, "family": family, "baltimore": baltimore, "branch": branch,
                     "nuc": nuc, "rate": to_nuc_rate / potential, "cutoff":cutoff}, ignore_index=True)

        all_ratios_to_nucs = all_ratios_to_nucs.append(ratios_to_nucs, ignore_index=True)
        all_summaries = all_summaries.append(summary, ignore_index = True)

    all_summaries.to_csv(summary_output)
    all_ratios_to_nucs.to_csv(ratios_output)


def get_mutation_joint_df_codon(basename, infoDictionary={},  overwrite = False):
        phyVirus_analysis_dir = "/sternadi/home/volume3/taliakustin/phyVirus_analysis/"
        fastml_dir = "fastml_codon_midpoint_root/"
        output = phyVirus_analysis_dir + fastml_dir + basename + ".joint.mutation.map"
        branch_file = glob.glob(phyVirus_analysis_dir + "/branch_length_info/" + basename + ".*")[0]
        if os.path.isfile(output) and not overwrite:
            print("output file already exists - returns. %s" % output)
            return output
        seq_joint = glob.glob(phyVirus_analysis_dir + fastml_dir + basename + "*seq.joint.txt")[0]
        tree_ancestor = glob.glob(phyVirus_analysis_dir + fastml_dir + basename + "*tree.ancestor.txt")[0]

        # create ancestor info - for each son - who is the father
        pattern = re.compile(r'\s+')
        ancestors = open(tree_ancestor, "r").read().split("\n")[2:-2]
        ancestor_info = {}
        for line in ancestors:
            line = re.sub(pattern, '$', line)
            line = line.split("$")
            son = line[0]
            father = line[1]
            ancestor_info[son] = father

        # create seq dictionary - name and sequance
        seqs = open(seq_joint, "r").read()
        seqs = seqs.split("\n>")[1:]
        seqs = {seq.split("\n")[0]: seq.split("\n")[1] for seq in seqs}

        res = pd.DataFrame()

        for son in ancestor_info:
            pattern = re.compile("^N\d*")  # checks if the node is an internal node
            if pattern.findall(son) == []:
                node_type = "external"
            else:
                node_type = "internal"

            father = ancestor_info[son]
            if father == "root!":
                continue
            son_seq = (seqs[son])
            father_seq = (seqs[father])

            son_triplets = [son_seq[i:i + 3] for i in range(0, len(son_seq), 3)]
            father_triplets = [father_seq[i:i + 3] for i in range(0, len(father_seq), 3)]
            no_diffs = ""
            diffs=[]
            for i in range(len(son_triplets)):
                if son_triplets[i] != father_triplets[i] and "-" not in son_triplets[i] and "-" not in father_triplets[i]:
                    diffs.append([father_triplets[i], son_triplets[i]])
                elif son_triplets[i] == father_triplets[i] and "-" not in son_triplets[i] and "-" not in father_triplets[i]:
                    no_diffs += father_triplets[i]

            muts = ["%s%s" % (no_diffs[i], no_diffs[i]) for i in range(len(no_diffs))]

            temp_res = pd.DataFrame(index=range(len(muts)),
                                    columns=["father_node_name", "son_node_name", "mutation",
                                             "branch"])
            temp_res.mutation = muts
            temp_res.father_node_name = father
            temp_res.son_node_name = son
            temp_res.branch = node_type

            for i in diffs:
                if Seq(i[0]).translate() == Seq(i[1]).translate():
                    mut = "syn"
                else:
                    mut = "non_syn"
                changes = [[j, f"{i[0][j]}{i[1][j]}"] for j in range(len(i[1])) if i[0][j] != i[1][j]]
                if len(changes) > 1: #more than 1 change - can't know if it's synonymous or not
                    mut = "unknown"

                for j in range(len(i[1])):
                    if i[0][j] == i[1][j]:
                        temp_res = temp_res.append({"father_node_name":father, "son_node_name":son,
                                                    "mutation":f"{i[0][j]}{i[1][j]}", "branch":node_type},  ignore_index=True)
                    else:
                        temp_res = temp_res.append({"father_node_name": father, "son_node_name": son,
                                                    "mutation": f"{i[0][j]}{i[1][j]}", "branch":node_type,
                                                    "mut_type":mut, "position":j}, ignore_index=True)
            res = res.append(temp_res, ignore_index=True)

        for key in infoDictionary.keys():
            res[key] = infoDictionary[key]

        res["node_name"] = res["son_node_name"]
        branch = pd.read_csv(branch_file)
        res_merged = res.merge(branch, how='outer', on=['node_name', "base", "family", "baltimore"])

        res_merged.to_csv(output)
        return (output)


def get_mutation_context(basename, overwrite=False, cutoff="no_cutoff"):
    phyVirus_analysis_dir = "/sternadi/home/volume3/taliakustin/phyVirus_analysis/"
    family = basename.split("_")[0]
    baltimore = get_baltimore_classifiaction(family)
    fastml_dir = "fastml_midpoint_tree/"
    output_f0 = phyVirus_analysis_dir + fastml_dir + basename + ".joint.context_f0.csv"
    output_f2 = phyVirus_analysis_dir + fastml_dir + basename + ".joint.context_f2.csv"
    res_f0 = pd.DataFrame()
    res_f2 = pd.DataFrame()
    branch_file = glob.glob(phyVirus_analysis_dir + "/branch_length_info/" + basename + ".*")[0]
    #if os.path.isfile(output_CU) and not overwrite:
    #    print("output file already exists - returns. %s" % output)
    #    #return pd.DataFrame()
    seq_joint = glob.glob(phyVirus_analysis_dir + fastml_dir + basename + ".seq.joint.txt")[0]
    tree_ancestor = glob.glob(phyVirus_analysis_dir + fastml_dir + basename + ".tree.ancestor.txt")[0]

    # create ancestor info - for each son - who is the father
    pattern = re.compile(r'\s+')
    ancestors = open(tree_ancestor, "r").read().split("\n")[2:-2]
    ancestor_info = {}
    for line in ancestors:
        line = re.sub(pattern, '$', line)
        line = line.split("$")
        son = line[0]
        father = line[1]
        ancestor_info[son] = father
    # create seq dictionary - name and sequance
    seqs = open(seq_joint, "r").read()
    seqs = seqs.split("\n>")[1:]
    seqs = {seq.split("\n")[0]: seq.split("\n")[1] for seq in seqs}

    res = pd.DataFrame()
    branch = pd.read_csv(branch_file)
    nucs = {"A": 0, "C": 0, "G": 0, "T": 0}
    dinucs = {"GA":0, "AG":0, "CT":0, "TC":0, "CA":0, "AC":0, "AT":0, "TA":0, "CG":0, "GC":0, "TG":0, "GT":0, "GG":0, "AA":0, "CC":0, "TT":0}

    for son in ancestor_info:
        father = ancestor_info[son]
        if father == "root!":
            continue

        son_seq = seqs[son]
        father_seq = seqs[father]
        pattern = re.compile("^N\d*") #checks if the node is an internal node

        if pattern.findall(son) == []:
            node_type = "external"
        else:
            node_type = "internal"
            continue
        if cutoff != "no_cutoff":
            if float(branch.loc[branch.node_name == son, "branch_len"]) > float(cutoff):
                continue

        #print(son)
        noGap_son_seq = []
        noGap_father_seq = []
        for i in range(len(son_seq)):
            if son_seq[i] != "-": #and father_seq[i] != "-":
                noGap_son_seq.append(son_seq[i])
                noGap_father_seq.append(father_seq[i])
        if len(noGap_son_seq) != len(noGap_father_seq):
            raise Error("not the same length")
        for i in range(len(noGap_father_seq)-1):
            if noGap_father_seq[i] != noGap_son_seq[i]:
                if noGap_father_seq[i-1] != noGap_son_seq[i-1] or noGap_father_seq[i+1] != noGap_son_seq[i+1]:
                    #print(son, noGap_father_seq[i-1:i+2], noGap_son_seq[i-1:i+2])
                    continue
                res =  res.append({"node_name":son, "mut":f"{noGap_father_seq[i]}{noGap_son_seq[i]}", "letter1":noGap_father_seq[i], "letter2":noGap_son_seq[i],
                                  "f0":noGap_father_seq[i-1], "f2":noGap_father_seq[i+1], "father_node_name":father, "son_node_name":son, "branch":"external",
                                   "base":basename, "family":family, "baltimore":baltimore}, ignore_index=True)


        str_father_seq = "".join(noGap_father_seq)
        for n in nucs.keys():
            nucs[n] += str_father_seq.count(n)
        for s in dinucs.keys():
            dinucs[s] += str_father_seq.count(s)
    print(dinucs)

    if res.empty:
        print("EMPTY!!!!")
        return pd.DataFrame()

    nucs_sum = sum(nucs.values())
    dinucs_sum = sum(dinucs.values())

    nuc_freqs = {}
    for n in nucs.keys():
        nuc_freqs[n] = nucs[n] / nucs_sum
    dinucs_freqs = {}
    for d in dinucs.keys():
        dinucs_freqs[d] = dinucs[d] / dinucs_sum
    print(dinucs_freqs)

    subs = ["GA", "AG", "CT", "TC", "CA", "AC", "AT", "TA", "CG", "GC", "TG", "GT"]

    for s in subs:
        s_context = res.loc[res.mut == s]
        if s_context.empty:
            continue
        s_context.loc[:, 'count'] = 1
        s_context_f0 = s_context.groupby(["branch", "f0", "base", "family", "baltimore", "mut"], as_index=False).agg({"count": "count"})
        s_context_f2 = s_context.groupby(["branch", "f2", "base", "family", "baltimore", "mut"], as_index=False).agg({"count": "count"})

        s_count_f0 = s_context_f0.loc[:, "count"].sum()
        s_count_f2 = s_context_f2.loc[:, "count"].sum()

        s_context_f0["freq"] = s_context_f0["count"] / s_count_f0
        s_context_f2["freq"] = s_context_f2["count"] / s_count_f2



        from_nuc = s[0]
        to_nuc = s[1]


        for n in nucs:
            #s_context_f0.loc[s_context_f0.f0 == n, "ratio"] = s_context_f0.loc[s_context_f0.f0 == n, "freq"] / (
            #            nucs[n] / nucs_sum)
            #s_context_f2.loc[s_context_f2.f2 == n, "ratio"] = s_context_f2.loc[s_context_f2.f2 == n, "freq"] / (
            #        nucs[n] / nucs_sum)
            s_context_f0.loc[s_context_f0.f0 == n, "ratio"] = s_context_f0.loc[s_context_f0.f0 == n, "freq"] / (
                        dinucs_freqs[n+from_nuc] / (nuc_freqs[n]*nuc_freqs[from_nuc]))
            s_context_f2.loc[s_context_f2.f2 == n, "ratio"] = s_context_f2.loc[s_context_f2.f2 == n, "freq"] / (
                    dinucs_freqs[from_nuc + n] / (nuc_freqs[n] * nuc_freqs[from_nuc]))

        sum_ratios_f0 = s_context_f0.loc[:, "ratio"].sum()
        sum_ratios_f2 = s_context_f2.loc[:, "ratio"].sum()

        s_context_f0.loc[:, "norm_ratio"] = s_context_f0.loc[:, "ratio"] / sum_ratios_f0
        s_context_f2.loc[:, "norm_ratio"] = s_context_f2.loc[:, "ratio"] / sum_ratios_f2
        if s_count_f0  > 20:
            res_f0 = res_f0.append(s_context_f0, ignore_index = True)
        if s_count_f2 > 20:
            res_f2 = res_f2.append(s_context_f2, ignore_index=True)
    res_f0.to_csv(output_f0)
    res_f2.to_csv(output_f2)





def get_mutation_context_apobec_contingency(basename, overwrite=False, cutoff="no_cutoff"):
    phyVirus_analysis_dir = "/sternadi/home/volume3/taliakustin/phyVirus_analysis/"
    fastml_dir = "fastml_midpoint_tree/"
    output = phyVirus_analysis_dir + fastml_dir + basename + ".joint.apobec_context.csv"
    branch_file = glob.glob(phyVirus_analysis_dir + "/branch_length_info/" + basename + ".*")[0]
    if os.path.isfile(output) and not overwrite:
        print("output file already exists - returns. %s" % output)
        return output
    seq_joint = glob.glob(phyVirus_analysis_dir + fastml_dir + basename + ".seq.joint.txt")[0]
    tree_ancestor = glob.glob(phyVirus_analysis_dir + fastml_dir + basename + ".tree.ancestor.txt")[0]

    # create ancestor info - for each son - who is the father
    pattern = re.compile(r'\s+')
    ancestors = open(tree_ancestor, "r").read().split("\n")[2:-2]
    ancestor_info = {}
    for line in ancestors:
        line = re.sub(pattern, '$', line)
        line = line.split("$")
        son = line[0]
        father = line[1]
        ancestor_info[son] = father
    # create seq dictionary - name and sequance
    seqs = open(seq_joint, "r").read()
    seqs = seqs.split("\n>")[1:]
    seqs = {seq.split("\n")[0]: seq.split("\n")[1] for seq in seqs}
    alignment_len = (len(seqs["N1"]))
    GA_VS_NOT = pd.DataFrame([[0, 0], [0,0]], index=["GA", "other"], columns=["Apobec_context", "not_apobec_context"])


    res = pd.DataFrame()
    for son in ancestor_info:
        father = ancestor_info[son]
        if father == "root!":
            continue
        son_seq = seqs[son]
        father_seq = seqs[father]
        pattern = re.compile("^N\d*") #checks if the node is an internal node
        if pattern.findall(son) == []:
            node_type = "external"
        else:
            node_type = "internal"

        temp_res = pd.DataFrame()


        noGap_son_seq = []
        noGap_father_seq = []
        for i in range(len(son_seq)):
            if son_seq[i] != "-": #and father_seq[i] != "-":
                noGap_son_seq.append(son_seq[i])
                noGap_father_seq.append(father_seq[i])

        if len(noGap_son_seq) != len(noGap_father_seq):
            raise Error("not the same length")
        for i in range(len(noGap_father_seq)-1):
            if noGap_father_seq[i] != noGap_son_seq[i]:
                if noGap_father_seq[i-1] != noGap_son_seq[i-1] or noGap_father_seq[i+1] != noGap_son_seq[i+1]:
                    continue
                temp_res =  temp_res.append({"node_name":son, "mut":f"{noGap_father_seq[i]}{noGap_son_seq[i]}", "letter1":noGap_father_seq[i], "letter2":noGap_son_seq[i],
                                  "f0":noGap_father_seq[i-1], "f2":noGap_father_seq[i+1]}, ignore_index=True)


        temp_res["father_node_name"] = father
        temp_res["son_node_name"] = son
        temp_res["no_cutoff"] = node_type
        res = res.append(temp_res, ignore_index=True)


    res["base"] = basename
    family = basename.split("_")[0]
    res["family"] = family
    baltimore = get_baltimore_classifiaction(family)
    res["baltimore"] = baltimore
    res["GA"] = False
    res.loc[res.mut == "GA", "GA"] = True
    res["CU"] = False
    res.loc[res.mut == "CT", "CU"] = True
    res["GA_context"] = False
    res.loc[(res.f2 == "A"), "GA_context"] = True
    res.loc[(res.f2 == "G"), "GA_context"] = True
    res["CU_context"] = False
    res.loc[(res.f0 == "C"), "CU_context"] = True
    res.loc[(res.f0 == "T"), "CU_context"] = True
    res["count"] = 1

    branch = pd.read_csv(branch_file)
    res_merged = res.merge(branch, how='outer', on=['node_name', "base", "family", "baltimore"])

    for c in [0.05, 0.04, 0.03, 0.02, 0.01, 0.005, 0.001]:
        res_merged[str(c)] = "internal"
        res_merged.loc[res_merged.branch_len < c, str(c)] = "external"
    res_merged["branch"] = res_merged[str(cutoff)]
    print(res_merged)
    GA = res_merged.groupby(["branch", "GA", "GA_context",  "base", "family", "baltimore"], as_index=False).agg({"count": "count"})
    CU = res_merged.groupby(["branch", "CU", "CU_context", "base", "family", "baltimore"], as_index=False).agg({"count": "count"})

    mut_context = GA.loc[(GA.branch == "external") & (GA.GA == True) & (GA.GA_context == True), "count"].sum()
    mut_non_context = GA.loc[(GA.branch == "external") & (GA.GA == True) & (GA.GA_context == False), "count"].sum()
    non_mut_context = GA.loc[(GA.branch == "external") & (GA.GA == False) & (GA.GA_context == True), "count"].sum()
    non_mut_non_context = GA.loc[(GA.branch == "external") & (GA.GA == False) & (GA.GA_context == False), "count"].sum()

    con = np.array([[mut_context, mut_non_context], [non_mut_context, non_mut_non_context]])
    table = sm.stats.Table(con)
    rslt = table.test_nominal_association()
    pval = rslt.pvalue
    statistic = rslt.statistic
    resid = table.resid_pearson
    resid_mut_context = resid[0][0]
    resid_mut_non_context = resid[0][1]
    resid_non_mut_context = resid[1][0]
    resid_non_nut_non_context = resid[1][1]

    if mut_non_context * non_mut_context == 0:
        OR = 0
    else:
        OR = (mut_context * mut_non_context) / (mut_non_context * non_mut_context)
    cont_res = pd.DataFrame()
    cont_res = cont_res.append({"base": basename, "family": family, "baltimore": baltimore, "sub": "GA",
                                "mut_context": mut_context, "mut_non_context": mut_non_context,
                                "non_mut_context": non_mut_context,
                                "non_mut_non_context": non_mut_non_context,
                                "resid_mut_context": resid_mut_context,
                                "resid_mut_non_context": resid_mut_non_context,
                                "resid_non_mut_context": resid_non_mut_context,
                                "resid_non_nut_non_context": resid_non_nut_non_context,
                                "p_value_not_corrected": pval,
                                "statistic": statistic, "OR": OR},
                               ignore_index=True)

    mut_context = CU.loc[(CU.branch == "external") & (CU.CU == True) & (CU.CU_context == True), "count"].sum()
    mut_non_context = CU.loc[(CU.branch == "external") & (CU.CU == True) & (CU.CU_context == False), "count"].sum()
    non_mut_context = CU.loc[(CU.branch == "external") & (CU.CU == False) & (CU.CU_context == True), "count"].sum()
    non_mut_non_context = CU.loc[(CU.branch == "external") & (CU.CU == False) & (CU.CU_context == False), "count"].sum()

    con = np.array([[mut_context, mut_non_context], [non_mut_context, non_mut_non_context]])
    table = sm.stats.Table(con)
    rslt = table.test_nominal_association()
    pval = rslt.pvalue
    statistic = rslt.statistic
    resid = table.resid_pearson
    resid = table.resid_pearson
    resid_mut_context = resid[0][0]
    resid_mut_non_context = resid[0][1]
    resid_non_mut_context = resid[1][0]
    resid_non_nut_non_context = resid[1][1]

    if mut_non_context * non_mut_context == 0:
        OR = 0
    else:
        OR = (mut_context * mut_non_context) / (mut_non_context * non_mut_context)

    cont_res = cont_res.append({"base": basename, "family": family, "baltimore": baltimore, "sub": "CU",
                                "mut_context": mut_context, "mut_non_context": mut_non_context,
                                "non_mut_context": non_mut_context,
                                "non_mut_non_context": non_mut_non_context,
                                "resid_mut_context": resid_mut_context,
                                "resid_mut_non_context": resid_mut_non_context,
                                "resid_non_mut_context": resid_non_mut_context,
                                "resid_non_nut_non_context": resid_non_nut_non_context,
                                "p_value_not_corrected": pval,
                                "statistic": statistic, "OR": OR},
                               ignore_index=True)
    cont_res.to_csv(output)
    return cont_res