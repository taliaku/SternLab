#! /powerapps/share/python-anaconda-3.6/bin/python

import os
import glob
from file_utilities import check_dirname, check_filename
import pandas as pd
import re
from PAML_utilities import retrive_kappa_from_paml_output
from phyVirus.get_baltimore import get_baltimore_classifiaction
from scipy.stats.distributions import chi2
from statsmodels.stats.multitest import multipletests
import numpy as np
from io import StringIO
from Bio import AlignIO, SeqIO
from Bio.Seq import Seq
from Bio import Phylo
import tqdm



def write_params_file(param_file, input_aln, input_tree, output_res, output_log, output_tree, gtr_output=None,
                      fixed_beta=0, fixed_s = 0, fixed_ProbS = 1, init_probS = 0.01, fixed_tau=1, init_beta=0,  fixed_kappa=1, bblOpt=0, optimizeLineSearch=0):

    param_file = check_filename(param_file, Truefile=False)
    input_aln = check_filename(input_aln)
    input_tree = check_filename(input_tree)
    output_res = check_filename(output_res, Truefile=False)
    output_log = check_filename(output_log, Truefile=False)
    output_tree =  check_filename(output_tree, Truefile=False)

    kappa = 2.0
    if gtr_output != None:
        gtr_output = check_filename(gtr_output)
        kappa_from_gtr = retrive_kappa_from_paml_output(gtr_output)
        if kappa_from_gtr != None:
            kappa = kappa_from_gtr

    with open(input_aln, "r") as aln_handle:
        aln_data = aln_handle.readline()
        first_seq = aln_data.split(">")[1].strip()

    fixed_beta = int(fixed_beta)
    if not fixed_beta in [0, 1]:
        raise Exception("fixed beta must be 0 or 1")

    fixed_s = int(fixed_s)
    if not fixed_s in [0, 1]:
        raise Exception("fixed_s must be 0 or 1")


    params = open(param_file, "w")
    params.write("# input\n")
    params.write("_inSeqFile %s\n" % input_aln)
    params.write("_inTreeFile %s\n" % input_tree)
    params.write("_inQuerySeq %s\n" % first_seq)
    params.write("_rootAt %s\n" % first_seq)
    params.write("_useQueryFreqsAtRoot 0\n")
    params.write("\n")
    params.write("# output\n")
    params.write("_outResFile %s\n" % output_res)
    params.write("_logFile %s\n" % output_log)
    params.write("_outTreeFile %s\n" % output_tree)
    params.write("\n")
    params.write("# advanced: remove # to enable advanced parameter. Parameters are described in the ParaSel manual at https://www.sternadi.com/parasel\n")
    params.write("# _modelName hky\n")
    params.write("_fixedS %i\n" % fixed_s)
    params.write("_initS 1.0\n")
    params.write("_fixedProbS %i\n" % fixed_ProbS)
    params.write("_initProbS %f\n" % init_probS)
    params.write("_initKappa %f\n" % kappa)
    params.write("_fixedKappa %i\n" % fixed_kappa)
    params.write("_initAlpha 1\n")
    params.write("_fixedAlpha 0\n")
    params.write("_initBeta 0\n")
    params.write("_fixedBeta %i\n" % fixed_beta)
    params.write("_initTau 1\n")
    params.write("_fixedTau %i\n" % fixed_tau)
    params.write("_bblOpt %i\n" % bblOpt)
    params.write("_doMutationMapping 0\n")
    params.write("_optimizeLineSearch %i\n" % optimizeLineSearch)
    params.write("#_threshold 0.75\n")

    params.close()
    return params


def return_branch_type(row):
    """
    check for each node if it's internal or external
    :param row: dataframe row
    :return: internal or external
    """
    #if the node name is N## than it's an internal node
    if re.match("^N\d+", row["node_name"]) == None:
        return "external"
    return "internal"

def mut_or_not(row):
    if row["letter1"] == row["letter2"]:
        return "no_mut"
    return "mut"

def get_general_info(file):
    base = file.split("/")[-1].split(".dirSel.results.mutation.map")[0]
    family = base.split("_")[0]
    baltimore = get_baltimore_classifiaction(family)
    return(base, family, baltimore)

def analyze_dirSel_mutation_map(file, ratios_output=None, overwrite = False):
    file = check_filename(file)
    if ratios_output == None:
        ratios_output = file.split(".dirSel.results.mutation.map")[0] + ".dirSel.ratios"
    else:
        ratios_output = check_filename(ratios_output, Truefile=False)
    if os.path.isfile(ratios_output) and not overwrite:
        print("output file %s exists" % ratios_output)
        return
    base, family, baltimore = get_general_info(file)
    mapping = {0:"A", 1:"C", 2:"G", 3:"T"}
    df = pd.read_csv(file, sep="\t", index_col=False)
    if df.shape[0] == 0:
        print("%s - has not mutation mapping" % file)
        return
    df["branch"] = df.apply(return_branch_type, axis=1)
    df["count"] = 1
    summary = df.groupby(["letter1", "letter2", "branch"], as_index=False).agg({"count": "count"})
    for num,let in mapping.items():
        summary.loc[summary.letter1 == num, "letter1"] = let
        summary.loc[summary.letter2 == num, "letter2"] = let
    summary["basename"] = base
    summary["family"] = family
    summary["baltimore"] = baltimore

    ratios = pd.DataFrame()
    for branch in ["external", "internal"]:
        for nuc1 in mapping.values():
            count_from_nuc1 = sum(summary.loc[(summary.letter1 == nuc1) & (summary.branch == branch), "count"].values)
            for nuc2 in mapping.values():
                if nuc1 == nuc2:
                    continue
                mut = nuc1 + nuc2
                mut_count = list(summary.loc[(summary.letter1 == nuc1) &
                                    (summary.letter2 == nuc2) & (summary.branch == branch), "count"].values)
                if mut_count == []:
                    mut_count.append(1)
                ratio = float(mut_count[0]) /  count_from_nuc1


                ratios = ratios.append({"basename":base, "family":family, "baltimore":baltimore, "branch":branch,
                                        "substitution":mut, "ratio":ratio}, ignore_index=True)
    ratios.to_csv(ratios_output)


def extract_dirSel_parameters_single_file(f):
    with open(f, "r") as handle:
        data = handle.readlines()
        line = 3  # start from line 3
        params = {}
        while not "#Rate categories are" in data[line]:
            if ":" in data[line]:
                name = data[line].split(":")[0].split("#")[1].strip()
                value = float(data[line].split(":")[1].strip())
            else:
                name = data[line].split("=")[0].split("#")[1].strip()
                value = float(data[line].split("=")[1].strip())
            params[name] = value
            line += 1

    return params

def extract_dirSel_parameters(files, output=None, overwrite=False, model=None, to_return=False):
    sig_files = pd.read_csv(
        "/sternadi/home/volume3/taliakustin/phyVirus_analysis/significant_incomplete_purifying_datasets.csv")
    if output != None:
        output = check_filename(output, Truefile=False)
    df = pd.DataFrame()
    for f in files:
        if "Reo" in f or "check" in f:
            continue
        base = f.split("/")[-1].split(".dirSel")[0]
        if sig_files.loc[sig_files.base == base].empty:
            continue
        family = base.split("_")[0]
        baltimore = get_baltimore_classifiaction(family)
        with open(f, "r") as handle:
            data =  handle.read()
            data = data.split("Parameters are:\n")[1].split("\n")
            line = 0 #start from line 3
            params = {}
            while not "#Rate categories are" in data[line]:
                if ":" in data[line]:
                    name = data[line].split(":")[0].split("#")[1].strip()
                    value = float(data[line].split(":")[1].strip())
                else:
                    name = data[line].split("=")[0].split("#")[1].strip()
                    value = float(data[line].split("=")[1].strip())
                params[name] = value
                line += 1
        df = df.append({"baltimore":baltimore, "family":family, "basename":base, "model":model, **params}, ignore_index=True)
    if to_return:
        return(df)
    df.to_csv(output)

def merge_alternative_and_null_dfs(alter, null, output=None, degrees_of_freedom=1, files=True, to_return=False):
    if files:
        alter = check_filename(alter)
        null = check_filename(null)
        alter = pd.read_csv(alter)
        null = pd.read_csv(null)
    if output !=None:
        output = check_filename(output, Truefile=False)
    merged = pd.merge(alter, null, on=['basename', "family", "baltimore"], suffixes=["_alter", "_null"])
    merged["diff"] = merged["Log-likelihood_alter"] - merged["Log-likelihood_null"]
    merged["diff2"] = 2 * merged["diff"]
    merged["p_value_not_corrected"] = chi2.sf(merged["diff2"], degrees_of_freedom)
    merged["p_value"] = multipletests(merged["p_value_not_corrected"], method="fdr_bh")[1]
    merged["sig_by_chi2"] = np.where(merged["p_value"] <= 0.05, True, False)
    merged["sig_by_diff"] = np.where(merged["diff2"] >= 100, True, False)
    if to_return:
        return(merged)
    merged.to_csv(output)


def analyze_dirSel_mutation_map_to_nucs(file, ratios_output=None, summary_otuput = None, overwrite = False):
    file = check_filename(file)
    if ratios_output == None:
        ratios_output = file.split(".dirSel.results.mutation.map")[0] + ".dirSel.ratios_to_nucs"
        summary_output = file.split(".dirSel.results.mutation.map")[0] + ".dirSel.mutation_summary"
    else:
        ratios_output = check_filename(ratios_output, Truefile=False)
        summary_output = check_filename(summary_otuput)
    if os.path.isfile(ratios_output) and not overwrite:
        print("output file %s exists" % ratios_output)
        return
    base, family, baltimore = get_general_info(file)
    mapping = {0:"A", 1:"C", 2:"G", 3:"T"}
    df = pd.read_csv(file, sep="\t", index_col=False)
    if df.shape[0] == 0:
        print("%s - has not mutation mapping" % file)
        return
    df["branch"] = df.apply(return_branch_type, axis=1)
    df["count"] = 1
    df["mut"] = df.apply(mut_or_not, axis=1)

    summary = df.groupby(["letter1", "letter2", "branch", "mut"], as_index=False).agg({"count": "count"})
    for num, let in mapping.items():
        summary.loc[summary.letter1 == num, "letter1"] = let
        summary.loc[summary.letter2 == num, "letter2"] = let
    summary["basename"] = base
    summary["family"] = family
    summary["baltimore"] = baltimore

    ratios_to_nucs = pd.DataFrame()
    for branch in ["external", "internal"]:
        for nuc in mapping.values():
            to_nuc_rate = sum(summary.loc[(summary["letter2"]==nuc) & (summary["letter1"] != nuc) & (summary["branch"] == branch), "count"].values)  \
                          / sum(summary.loc[(summary["letter1"]!=nuc)   & (summary["branch"] == branch), "count"].values)

            ratios_to_nucs = ratios_to_nucs.append({"basename": base, "family": family, "baltimore": baltimore, "branch": branch,
                                    "nucleotide": nuc, "to_nucleotide_rate": to_nuc_rate}, ignore_index=True)
    print(ratios_output)
    ratios_to_nucs.to_csv(ratios_output)

def analyze_dirSel_mutation_map_to_nucs_updated(file, ratios_output=None, overwrite=False):
    file = check_filename(file)
    print(file)
    if ratios_output == None:
        ratios_output = file.split(".dirSel.results.mutation.map")[0] + ".dirSel.ratios_to_nucs_updated"
    else:
        ratios_output = check_filename(ratios_output, Truefile=False)
    if os.path.isfile(ratios_output) and not overwrite:
        print("output file %s exists" % ratios_output)
        return
    base, family, baltimore = get_general_info(file)

    mapping = {0: "A", 1: "C", 2: "G", 3: "T"}
    df = pd.read_csv(file, sep="\t", index_col=False)
    if df.shape[0] == 0:
        print("%s - has not mutation mapping" % file)
        return
    df["branch"] = df.apply(return_branch_type, axis=1)
    df["count"] = 1
    df["mut"] = df.apply(mut_or_not, axis=1)

    summary = df.groupby(["letter1", "letter2", "branch", "mut"], as_index=False).agg({"count": "count"})
    for num, let in mapping.items():
        summary.loc[summary.letter1 == num, "letter1"] = let
        summary.loc[summary.letter2 == num, "letter2"] = let
    summary["basename"] = base
    summary["family"] = family
    summary["baltimore"] = baltimore

    ratios_to_nucs = pd.DataFrame()
    for branch in ["external", "internal"]:
        for nuc in mapping.values():
            if (sum(summary.loc[(summary["mut"]=="mut") & (summary["branch"] == branch), "count"].values)) != 0:
                #rate of mutation to this nuc from all mutations to all nucs
                to_nuc_rate = sum(summary.loc[(summary["letter2"] == nuc) & (summary["letter1"] != nuc) & (
                summary["branch"] == branch), "count"].values) \
                              / sum(summary.loc[(summary["mut"]=="mut") & (summary["branch"] == branch), "count"].values)
            else:
                to_nuc_rate = 0
            #num of position that can go to nuc from all positions available to mutate
            potential = sum(summary.loc[(summary["letter1"] != nuc) & (summary["branch"] == branch), "count"].values)\
                        / sum(summary.loc[(summary["branch"] == branch), "count"].values)


            ratios_to_nucs = ratios_to_nucs.append(
                {"basename": base, "family": family, "baltimore": baltimore, "branch": branch,
                 "nuc": nuc, "rate": to_nuc_rate/potential}, ignore_index=True)
    ratios_to_nucs.to_csv(ratios_output)


def analyze_context_dirSel(file, letter1="G", letter2="A", before_or_after = "after", output=None, overwrite=False, save_or_return="save", aln_file=None):
    mutation = letter1 + letter2
    file = check_filename(file)
    if before_or_after not in ["before", "after"]:
        raise ("before_or_after has to be before / after")
    if save_or_return not in ["save", "return"]:
        raise ("save_or_return has to be save / return")

    if save_or_return == "save":
        if output == None:
            output = file.split(".dirSel.results.mutation.map")[0] + ".dirSel.context_%s%s_%s" % (letter1, letter2, before_or_after)
        else:
            output = check_filename(output, Truefile=False)
        if os.path.isfile(output) and not overwrite:
            print("output file %s exists" % output)
            return
    base, family, baltimore = get_general_info(file)
    dirSel_result = file.split(".mutation.map")[0]
    if aln_file == None:
        aln_file = glob.glob("/sternadi/home/volume3/taliakustin/phyVirus_analysis/aln/%s.*fas" % base)[0]
    else:
        aln_file = check_filename(aln_file)

    mapping = {0: "A", 1: "C", 2: "G", 3: "T"}
    df = pd.read_csv(file, sep="\t", index_col=False)
    df["branch"] = df.apply(return_branch_type, axis=1)
    df["count"] = 1
    df["mut"] = df.apply(mut_or_not, axis=1)
    for num, let in mapping.items():
        df.loc[df.letter1 == num, "letter1"] = let
        df.loc[df.letter2 == num, "letter2"] = let
    print(aln_file)

    #map between query seq and pos in aln
    with open(dirSel_result) as handle:
        data = handle.read()
        data = data.split("POSTERIOR PROBABILITY OF DIRECTIONAL SELECTION \n\t\t\tA\tC\tG\tT\t\n")[-1]
        data = "query\taln\tresidue\tA\tC\tG\tT\tOther\t\n" + data
        res = pd.read_csv(StringIO(data), sep="\t", encoding='utf8', index_col=False)


    muts_ext = df.loc[(df["letter1"] == letter1) & (df["letter2"] == letter2) & (df["branch"] == "external"), ["#Pos(query)", "node_name"]]
    if muts_ext.empty:
        print("no mutations in external branches from %s to %s" % (letter1, letter2))
        return
    node_names = set(muts_ext["node_name"])
    nps = {} #names and positions of mutations
    for n in node_names:
        pos_query = muts_ext.loc[muts_ext["node_name"] == n, "#Pos(query)"].values
        nps[n] = []
        for i in pos_query:
            nps[n].append(res.loc[res["query"] == i, "aln"].values[0] - 1)
    aln = AlignIO.read(aln_file, "fasta")
    if before_or_after == "after":
        context_nucs = {"%sA" % letter1:0, "%sC" % letter1:0, "%sG" % letter1:0, "%sT" % letter1:0}
        potential_contexts = {"%sA" % letter1:0, "%sC" % letter1:0, "%sG" % letter1:0, "%sT" % letter1:0}
    elif before_or_after == "before":
        context_nucs = {"A%s" % letter1:0, "C%s" % letter1:0, "G%s" % letter1:0, "T%s" % letter1:0}
        potential_contexts = {"A%s" % letter1:0, "C%s" % letter1:0, "G%s" % letter1:0, "T%s" % letter1:0}

    for s in aln:
        if s.id in nps:
            for i in nps[s.id]:
                if s.seq[i] != letter2:
                    #print (s.id, i, s.seq[i], letter2)
                    if s.seq[i] in ["W", "M"]:
                        continue
                    if "dirSel.aln" in aln_file:
                        pass
                    else:
                        raise ("Problem with mutation indexing")
                if before_or_after == "after":
                    if len(s.seq) == i+1:
                        print("last position %i" % i)
                        continue
                    elif str(s.seq[i+1:]).replace("-", "") == "":
                        print("last position %i" % i)
                        continue
                    context_nuc = str(s.seq[i + 1:]).replace("-", "")[0]
                    if context_nuc not in mapping.values():
                        continue
                    context_nucs["%s%s" % (letter1, context_nuc)] += 1
                elif before_or_after == "before":
                    if i == 0:
                        print("first position %i" % i)
                        continue
                    elif str(s.seq[:i-1]).replace("-", "") == "":
                        print("first position %i" % i)
                        continue
                    context_nuc = str(s.seq[:i-1]).replace("-", "")[-1]
                    if context_nuc not in mapping.values():
                        continue
                    context_nucs["%s%s" % (context_nuc, letter1)] += 1
        for context in potential_contexts.keys():
            potential_contexts[context] += str(s.seq).replace("-", "").count(context)


    context_rate = {}
    context_ratio = {}
    df = pd.DataFrame()
    print(context_nucs)
    print(potential_contexts)
    for context in context_nucs.keys():
        context_rate[context] = context_nucs[context] / sum(context_nucs.values())
        context_ratio[context] = context_rate[context] / potential_contexts[context]
        df = df.append({"baltimore": baltimore, "family": family, "basename": base, "mutation":mutation, "context":context,
                        "before_or_after":before_or_after, "ratio":context_ratio[context]},
                       ignore_index=True)

    if save_or_return == "save":
        df.to_csv(output)
    else:
        return(df)

def run_several_context_analysis_on_file(file, output=None, overwrite=False, mutations=["CT", "GA", "AG", "TC"],  before_or_after = ["before", "after"], aln_type="original", aln_file=None):
    file = check_filename(file)
    if aln_type not in ["original", "from_mut_map"]:
        raise ("aln_type has to be original / from_mut_map")
    if output == None:
        if aln_type == "original":
            output = file.split(".dirSel.results.mutation.map")[0] + ".dirSel.context"
            aln_file = None
        else:
            output = file.split(".dirSel.results.mutation.map")[0] + ".dirSel.context_aln_mut_map"
            if aln_file == None:
                aln_file = file.split(".dirSel.results.mutation.map")[0] + ".dirSel.aln"
            else:
                aln_file = check_filename(aln_file)
    else:
        output = check_filename(output, Truefile=False)
    if os.path.isfile(output) and not overwrite:
        print("output file %s exists" % output)
        return
    df = pd.DataFrame()
    for mut in mutations:
        letter1 = mut[0]
        letter2 = mut[1]
        for loc in before_or_after:
            print("doing mut %s with context %s" % (mut, loc))
            df = df.append(analyze_context_dirSel(file, letter1, letter2, loc, save_or_return="return", aln_file=aln_file))
    df.to_csv(output)


def get_aln_pos_from_query_pos(row, res):
    return res.loc[res["query"] == row["#Pos(query)"], "aln"].values[0] - 1

def change_aln_according_to_mutations(row, aln):
    seq = list(aln[row["node_name"]].seq)
    seq[row["aln_pos"]] = row["letter1"]
    aln[row["node_name"]].seq = "".join(seq)
    return aln

def alignment_from_mutation_map_info(file, output=None):
    file = check_filename(file)
    if output == None:
        output = file.split(".dirSel.results.mutation.map")[0] + ".dirSel.aln"
    else:
        output = check_filename(output, Truefile=False)
    if os.path.isfile(output) and not overwrite:
        print("output file %s exists" % output)
        return
    base, family, baltimore = get_general_info(file)

    dirSel_result = file.split(".mutation.map")[0]
    aln_file = glob.glob("/sternadi/home/volume3/taliakustin/phyVirus_analysis/aln/%s.*fas" % base)[0]
    mapping = {0: "A", 1: "C", 2: "G", 3: "T"}
    df = pd.read_csv(file, sep="\t", index_col=False)
    df["branch"] = df.apply(return_branch_type, axis=1)
    df["count"] = 1
    df["mut"] = df.apply(mut_or_not, axis=1)
    for num, let in mapping.items():
        df.loc[df.letter1 == num, "letter1"] = let
        df.loc[df.letter2 == num, "letter2"] = let

    with open(dirSel_result) as handle:
        data = handle.read()
        data = data.split("POSTERIOR PROBABILITY OF DIRECTIONAL SELECTION \n\t\t\tA\tC\tG\tT\t\n")[-1]
        data = "query\taln\tresidue\tA\tC\tG\tT\tOther\t\n" + data
        res = pd.read_csv(StringIO(data), sep="\t", encoding='utf8', index_col=False)

    muts_ext = df.loc[
        (df["branch"] == "external") & (df["letter1"] != df["letter2"]), ["#Pos(query)", "node_name", "letter1",
                                                                          "letter2"]]
    if muts_ext.empty:
        print("no mutations in external branches")
        return

    muts_ext["aln_pos"] = muts_ext.apply(get_aln_pos_from_query_pos, axis=1, res=res)

    aln = SeqIO.to_dict(SeqIO.parse(aln_file, "fasta"))
    for row in muts_ext.itertuples():
        seq_list = list(aln[row[2]].seq)
        if seq_list[row[5]] != row[4]:
            print(row)
            print(seq_list[row[5]])
        seq_list[row[5]] = row[3]
        aln[row[2]].seq = "".join(seq_list)
    for a in aln:
        if type(aln[a].seq) == str:
            aln[a].seq = Seq(aln[a].seq)
    with open(output, 'w') as handle:
        SeqIO.write(aln.values(), handle, 'fasta')

    print(output)

def analyze_dirSel_mutation_map_to_nucs_several_cutoffs(file, tree_file = None, cutoffs =["no_cutoff", 0.001, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4], ratios_output=None, overwrite=False):
    file = check_filename(file)
    print(file)
    if ratios_output == None:
        ratios_output = file.split(".dirSel.results.mutation.map")[0] + ".dirSel.ratios_to_nucs_cutoffs"
    else:
        ratios_output = check_filename(ratios_output, Truefile=False)
    if os.path.isfile(ratios_output) and not overwrite:
        print("output file %s exists" % ratios_output)
        return
    base, family, baltimore = get_general_info(file)


    if tree_file != None:
        tree_file = check_filename(tree_file)
    else:
        tree_file = glob.glob("/sternadi/home/volume3/taliakustin/phyVirus_analysis/tree/%s*tree*" % base)[0]

    mapping = {0: "A", 1: "C", 2: "G", 3: "T"}
    df = pd.read_csv(file, sep="\t", index_col=False)
    if df.shape[0] == 0:
        print("%s - has not mutation mapping" % file)
        return
    df["branch"] = df.apply(return_branch_type, axis=1)
    df["count"] = 1
    df["mut"] = df.apply(mut_or_not, axis=1)

    tree = Phylo.read(tree_file, "newick")
    clades = list(tree.find_clades())
    for c in clades:
        if c.name != None:
            df.loc[df.node_name == c.name, "branch_length"] = c.branch_length

    ratios_to_nucs = pd.DataFrame()
    for cutoff in cutoffs:
        print(cutoff)
        if cutoff != "no_cutoff":
            df.loc[df.branch_length > cutoff, "branch"] = "internal"
            df.loc[df.branch_length <= cutoff, "branch"] = "external"

        if df.loc[df.branch=="external"].shape[0] == 0:
            print("no external branches in cutoff %f" % cutoff)
            continue

        summary = df.groupby(["letter1", "letter2", "branch", "mut"], as_index=False).agg({"count": "count"})
        for num, let in mapping.items():
            summary.loc[summary.letter1 == num, "letter1"] = let
            summary.loc[summary.letter2 == num, "letter2"] = let
        summary["basename"] = base
        summary["family"] = family
        summary["baltimore"] = baltimore

        for branch in ["external", "internal"]:
            for nuc in mapping.values():
                if (sum(summary.loc[(summary["mut"] == "mut") & (summary["branch"] == branch), "count"].values)) != 0:
                    # rate of mutation to this nuc from all mutations to all nucs
                    to_nuc_rate = sum(summary.loc[(summary["letter2"] == nuc) & (summary["letter1"] != nuc) & (summary["branch"] == branch), "count"].values) \
                                  / sum(
                        summary.loc[(summary["mut"] == "mut") & (summary["branch"] == branch), "count"].values)
                else:
                    to_nuc_rate = 0
                # num of position that can go to nuc from all positions available to mutate
                potential = sum(
                    summary.loc[(summary["letter1"] != nuc) & (summary["branch"] == branch), "count"].values) \
                            / sum(summary.loc[(summary["branch"] == branch), "count"].values)


                ratios_to_nucs = ratios_to_nucs.append(
                    {"basename": base, "family": family, "baltimore": baltimore, "branch": branch,
                     "nuc": nuc, "rate": to_nuc_rate / potential, "cutoff_branch_length":cutoff}, ignore_index=True)
    ratios_to_nucs.to_csv(ratios_output)


def make_distributaion_of_rates(file, tree_file=None, cutoffs =["no_cutoff", 0.001, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4]):
    file = check_filename(file)
    print(file)
    """
    if ratios_output == None:
        ratios_output = file.split(".dirSel.results.mutation.map")[0] + ".dirSel.ratios_to_nucs_updated"
    else:
        ratios_output = check_filename(ratios_output, Truefile=False)
    
    if os.path.isfile(ratios_output) and not overwrite:
        print("output file %s exists" % ratios_output)
        return
    """
    base, family, baltimore = get_general_info(file)

    if tree_file != None:
        tree_file = check_filename(tree_file)
    else:
        tree_file = glob.glob("/sternadi/home/volume3/taliakustin/phyVirus_analysis/tree/%s*tree*" % base)[0]

    df = pd.read_csv(file, sep="\t", index_col=False)
    if df.shape[0] == 0:
        print("%s - has not mutation mapping" % file)
        return
    df["branch"] = df.apply(return_branch_type, axis=1)
    df["count"] = 1
    df["mut"] = df.apply(mut_or_not, axis=1)

    tree = Phylo.read(tree_file, "newick")
    clades = list(tree.find_clades())
    for c in clades:
        if c.name != None:
            df.loc[df.node_name == c.name, "branch_length"] = c.branch_length

    results = pd.DataFrame()
    node_names = set(df["node_name"])
    for n in tqdm.tqdm(node_names):
        df_temp = df.loc[df.node_name == n]
        for c in cutoffs:
            if c == "no_cutoff":
                branch = np.unique(df_temp.branch)[0]
            elif np.unique(df_temp.branch)[0] == "internal" or np.unique(df_temp.branch_length)[0] >= c:
                branch="internal"
            else:
                branch = "external"

            results_tmp = calculate_relative_ratio_for_nodes(df_temp)
            results_tmp["node_name"] = n
            results_tmp["branch"] = branch
            results_tmp["branch_length"] = np.unique(df_temp.branch_length)[0]
            results_tmp["cutoff"] = c

            results = results.append(results_tmp)
    results["basename"] = base
    results["family"] = family
    results["baltimore"] = baltimore

    return results


def calculate_relative_ratio_for_nodes(df):
    mapping = {0: "A", 1: "C", 2: "G", 3: "T"}
    #summary = df.groupby(["node_name", "letter1", "letter2", "branch", "mut"], as_index=False).agg({"count": "count"})
    summary = df.groupby(["letter1", "letter2", "branch", "mut"], as_index=False).agg({"count": "count"})
    if sum(summary.loc[summary.mut=="mut", "count"].values) == 0:
        return pd.DataFrame()

    for num, let in mapping.items():
        summary.loc[summary.letter1 == num, "letter1"] = let
        summary.loc[summary.letter2 == num, "letter2"] = let

    relative_ratios = pd.DataFrame()
    for branch in ["external", "internal"]:
        if sum((summary.loc[summary["branch"] == branch, "count"].values)) == 0:
            continue
        for nuc in mapping.values():
            if (sum(summary.loc[(summary["mut"] == "mut") & (summary["branch"] == branch), "count"].values)) != 0:
                # rate of mutation to this nuc from all mutations to all nucs
                to_nuc_rate = sum(summary.loc[(summary["letter2"] == nuc) & (summary["letter1"] != nuc) & (summary["branch"] == branch), "count"].values) \
                              / sum(
                    summary.loc[(summary["mut"] == "mut") & (summary["branch"] == branch), "count"].values)
            else:
                continue
            # num of position that can go to nuc from all positions available to mutate
            potential = sum(
                summary.loc[(summary["letter1"] != nuc) & (summary["branch"] == branch), "count"].values) \
                        / sum(summary.loc[(summary["branch"] == branch), "count"].values)


            relative_ratios = relative_ratios.append({"branch": branch,
                 "nuc": nuc, "rate": to_nuc_rate / potential, "to_nuc_rate":to_nuc_rate}, ignore_index=True)

    return(relative_ratios)


def write_params_file_PQR(param_file, input_aln, input_tree, output_res, output_log, output_tree, gtr_output=None, fixed_beta=0,
                          fixed_s = 0, fixed_ProbS = 1, init_probS = 0.01, fixed_tau=1, fixed_kappa=1, bblOpt=0):

    param_file = check_filename(param_file, Truefile=False)
    input_aln = check_filename(input_aln)
    input_tree = check_filename(input_tree)
    output_res = check_filename(output_res, Truefile=False)
    output_log = check_filename(output_log, Truefile=False)
    output_tree =  check_filename(output_tree, Truefile=False)

    kappa = 2.0
    if gtr_output != None:
        gtr_output = check_filename(gtr_output)
        kappa_from_gtr = retrive_kappa_from_paml_output(gtr_output)
        if kappa_from_gtr != None:
            kappa = kappa_from_gtr

    with open(input_aln, "r") as aln_handle:
        aln_data = aln_handle.readline()
        first_seq = aln_data.split(">")[1].strip()

    fixed_beta = int(fixed_beta)
    if not fixed_beta in [0, 1]:
        raise Exception("fixed beta must be 0 or 1")

    fixed_s = int(fixed_s)
    if not fixed_s in [0, 1]:
        raise Exception("fixed_s must be 0 or 1")


    params = open(param_file, "w")
    params.write("# input\n")
    params.write("_inSeqFile %s\n" % input_aln)
    params.write("_inTreeFile %s\n" % input_tree)
    params.write("_inQuerySeq %s\n" % first_seq)
    params.write("_rootAt %s\n" % first_seq)
    params.write("_useQueryFreqsAtRoot 0\n")
    params.write("\n")
    params.write("# output\n")
    params.write("_outResFile %s\n" % output_res)
    params.write("_logFile %s\n" % output_log)
    params.write("_outTreeFile %s\n" % output_tree)
    params.write("\n")
    params.write("# advanced: remove # to enable advanced parameter. Parameters are described in the ParaSel manual at https://www.sternadi.com/parasel\n")
    params.write("# _modelName hky\n")
    params.write("_fixedS %i\n" % fixed_s)
    params.write("_initS 1.0\n")
    params.write("_fixedProbS %i\n" % fixed_ProbS)
    params.write("_initProbS %f\n" % init_probS)
    params.write("_initKappa %f\n" % kappa)
    params.write("_fixedKappa %i\n" % fixed_kappa)
    params.write("_initAlpha 1\n")
    params.write("_fixedAlpha 0\n")
    params.write("_initBeta 0\n")
    params.write("_fixedBeta %i\n" % fixed_beta)
    params.write("_initTau 1\n")
    params.write("_fixedTau %i\n" % fixed_tau)
    params.write("_bblOpt %i\n" % bblOpt)
    params.write("_doMutationMapping 0\n")
    params.write("#_threshold 0.75\n")
    params.write("_optimizeLineSearch 0\n")

    params.close()
    return params