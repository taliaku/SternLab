#! /powerapps/share/python-anaconda-3.6/bin/python

import sys
sys.path.insert(0,'/sternadi/home/volume1/taliakustin/SternLab')
from optparse import OptionParser
import os
import glob
from file_utilities import check_filename, check_dirname
from dirSel_utilities import write_params_file
import pbs_runners
from netMHCpan_utilities import netMHCpan_to_csv
from pandas_utilities import merge_dfs
import pandas as pd
import numpy as np
from phylogenetic_utilities import label_trees_by_cutoff
from phyVirus.phyVirus_utilities import  get_codon_alignment_phy_from_basename, get_rooted_tree_from_basename
import PAML_utilities
import fastml_utilities
from phyVirus.get_baltimore import  get_baltimore_classifiaction

def main():
    no_cutoff = "/sternadi/home/volume3/taliakustin/phyVirus_analysis/codeml_results/codeml_results_no_cutoff.09.02.2020.csv"
    cutoffs = {"0.05":"/sternadi/home/volume3/taliakustin/phyVirus_analysis/codeml_results/codeml_results_cutoff_0.05.09.02.2020.csv",
               "0.04":"/sternadi/home/volume3/taliakustin/phyVirus_analysis/codeml_results/codeml_results_cutoff_0.04.09.02.2020.csv",
               "0.03":"/sternadi/home/volume3/taliakustin/phyVirus_analysis/codeml_results/codeml_results_cutoff_0.03.09.02.2020.csv",
               "0.02":"/sternadi/home/volume3/taliakustin/phyVirus_analysis/codeml_results/codeml_results_cutoff_0.02.09.02.2020.csv",
               "0.01":"/sternadi/home/volume3/taliakustin/phyVirus_analysis/codeml_results/codeml_results_cutoff_0.01.09.02.2020.csv",
               "0.005":"/sternadi/home/volume3/taliakustin/phyVirus_analysis/codeml_results/codeml_results_cutoff_0.005.09.02.2020.csv",
               "0.001":"/sternadi/home/volume3/taliakustin/phyVirus_analysis/codeml_results/codeml_results_cutoff_0.001.09.02.2020.csv"}


    all_cutoffs = pd.DataFrame()
    for c in cutoffs:
        temp = pd.read_csv(cutoffs[c])
        temp["cutoff"] = c
        all_cutoffs = all_cutoffs.append(temp)


    df = pd.read_csv(no_cutoff)
    incomp_sig_df = df.loc[(df.sig_by_chi2 == True) & (df.wewi_sig == "over_0")]

    incomp_sig_df["cutoff"] = "no_cutoff"

    incomp_bases = {b:"no_cutoff" for b in df.loc[(df.sig_by_chi2 == True) & (df.wewi_sig == "over_0"), "base"]}

    negative_bases = df.loc[(df.sig_by_chi2 == False) | ((df.sig_by_chi2 == True) & (df.wewi_sig == "under_0")) | (df.diff_wewi > 10), "base"]
    bases_to_run_again = {}
    cutoffs_values = [0.05, 0.04, 0.03, 0.02, 0.01, 0.005, 0.001]

    for b in negative_bases:
        if all_cutoffs.loc[all_cutoffs.base==b].empty == True:
            bases_to_run_again[b] = cutoffs_values
        elif all_cutoffs.loc[(all_cutoffs.base==b) & (all_cutoffs.sig_by_chi2 == True) & (all_cutoffs.wewi_sig == "over_0") & (all_cutoffs.diff_wewi < 100)].empty == False:
            c = max(all_cutoffs.loc[(all_cutoffs.base==b) & (all_cutoffs.sig_by_chi2 == True) & (all_cutoffs.wewi_sig == "over_0"), "cutoff"])
            temp = all_cutoffs.loc[(all_cutoffs.base==b) & (all_cutoffs.sig_by_chi2 == True) & (all_cutoffs.wewi_sig == "over_0") & (all_cutoffs.diff_wewi < 100)].iloc[0]
            incomp_sig_df = incomp_sig_df.append(temp, ignore_index=True)
            incomp_bases[b] = c
        else:
            cutoffs_run = list(all_cutoffs.loc[all_cutoffs.base == b, "cutoff"])
            cutoffs_to_run = []
            for c in cutoffs_values:
                if str(c) not in cutoffs_run:
                    cutoffs_to_run.append(c)
            if cutoffs_to_run==[]:
                continue
            bases_to_run_again[b] = cutoffs_to_run

    count = 0
    for b in bases_to_run_again:
        base=b +"."
        aln = get_codon_alignment_phy_from_basename(base)
        tree = get_rooted_tree_from_basename(base)
        for c in bases_to_run_again[b]:
            out_dir = f"/sternadi/home/volume3/taliakustin/phyVirus_analysis/rtl_{str(c)}/"
            labled_tree = label_trees_by_cutoff(tree, c, out_dir)
            if labled_tree == False:
                #print(f"no tree {b}")
                continue
            else:

                ctl = f"/sternadi/home/volume3/taliakustin/phyVirus_analysis/c{str(c)}/{b}.ctl"
                mlc = f"/sternadi/home/volume3/taliakustin/phyVirus_analysis/c{str(c)}/{b}.mlc"
                if os.path.isfile(ctl):
                    continue
                count += 1
                PAML_utilities.write_ctl_codeml_deleterious_load_file(ctl=ctl, seq=aln, tree=labled_tree, out=mlc, model=2)
                pbs_runners.codeml_runner(queue="duduhimem", ctl=ctl, alias=f"codeml_{c}", cmdname=b)
                print(f"run comdel {b} cutoff {c}")


    print(count)

    incomp_sig_df["mutation_count"] = 0
    context = pd.DataFrame()

    for b in incomp_bases:
        if "SIV" in b or "Reo" in b:
            continue
        print(b,incomp_bases[b] )
        family = b.split("_")[0]
        baltimore = get_baltimore_classifiaction(family)
        map_file = fastml_utilities.get_mutation_joint_df(b, infoDictionary={"base": b, "family": family,
                                                                                      "baltimore": baltimore},
                                                                overwrite=False)
        summary_output = fastml_utilities.summmerize_mutations(map_file, cutoff=incomp_bases[b], overwrite=False)
        print(summary_output)
        summary = pd.read_csv(summary_output)
        mutation_count = sum(summary.loc[(summary.branch == "external") & (summary.mut == "mut"), "count"])
        incomp_sig_df.loc[incomp_sig_df.base == b, "mutation_count"] = mutation_count
        context = context.append(fastml_utilities.get_mutation_context(b, cutoff=incomp_bases[b]), ignore_index=True)
        context.to_csv("/sternadi/home/volume3/taliakustin/phyVirus_analysis/context_from_fastml.csv")



    incomp_sig_df.to_csv("/sternadi/home/volume3/taliakustin/phyVirus_analysis/significant_incomplete_purifying_datasets.csv")
    context.to_csv("/sternadi/home/volume3/taliakustin/phyVirus_analysis/context_CT_from_fastml.csv")


if __name__ == "__main__":
    main()

