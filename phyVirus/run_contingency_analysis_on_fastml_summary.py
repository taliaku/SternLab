#! /powerapps/share/python-anaconda-3.6/bin/python

import sys
sys.path.insert(0,'/sternadi/home/volume1/taliakustin/SternLab')
from optparse import OptionParser
import os
import glob
from file_utilities import check_filename, check_dirname
from phyVirus.get_baltimore import get_baltimore_classifiaction
import pandas as pd
import numpy as np
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests
from statsmodels_utilities import parse_summary_MH_test
import tqdm



def main():
    sig_incomplete_file = "/sternadi/home/volume3/taliakustin/phyVirus_analysis/significant_incomplete_purifying_datasets.csv"
    summary_file = "/sternadi/home/volume3/taliakustin/phyVirus_analysis/fastml_midpoint_tree_joint_mutation_summary.csv"




    tables_all = {"A":{}, "C":{}, "G":{}, "U":{}}
    all_nucs = {"A":pd.DataFrame(), "C":pd.DataFrame(), "G":pd.DataFrame(), "T":pd.DataFrame()}


    sig_incomplete = pd.read_csv(sig_incomplete_file)
    basenames = list(sig_incomplete.loc[:, "base"])
    print(len(basenames))
    count = 0

    for base in basenames:
        family = base.split("_")[0]
        baltimore = get_baltimore_classifiaction(family)
        cutoff = list(sig_incomplete.loc[sig_incomplete.base==base, "cutoff"])[0]
        summary_file = glob.glob(f"/sternadi/home/volume3/taliakustin/phyVirus_analysis/fastml_midpoint_tree/{base}.joint.mutation_summary")
        if summary_file == []:
            count += 1
            continue
        else:
            summary_file = summary_file[0]
        summary = pd.read_csv(summary_file)
        mutation_count = sum(summary.loc[(summary.branch=="external") & (summary.mut == "mut"), "count"])
        if  mutation_count < 50:
           continue


        for n in all_nucs.keys():
            internal_to = summary.loc[(summary.branch =="internal") & (summary.mut == "mut") & (summary.letter2 == n), "count"].sum()
            internal_non_to = summary.loc[(summary.branch =="internal") & (summary.mut == "mut") & (summary.letter2 != n), "count"].sum()
            external_to = summary.loc[(summary.branch =="external") & (summary.mut == "mut") & (summary.letter2 == n), "count"].sum()
            external_non_to = summary.loc[(summary.branch =="external") & (summary.mut == "mut") & (summary.letter2 != n), "count"].sum()
            con = np.array([[internal_to,internal_non_to], [external_to, external_non_to]])
            table = sm.stats.Table(con)
            rslt = table.test_nominal_association()
            pval = rslt.pvalue
            statistic = rslt.statistic
            resid = table.resid_pearson
            resid_internal_to =resid[0][0]
            resid_internal_non_to = resid[0][1]
            resid_external_to = resid[1][0]
            resid_external_non_to = resid[1][1]
            if internal_non_to*external_to == 0:
                OR = 0
            else:
                OR = (internal_to*external_non_to) / (internal_non_to*external_to)
            new_nuc = n
            if new_nuc == "T":
                new_nuc = "U"
            all_nucs[n] = all_nucs[n].append({"base":base, "family":family, "baltimore":baltimore, "nuc":new_nuc,
                              "internal_to":internal_to, "internal_non_to":internal_non_to,
                              "external_to":external_to, "external_non_to":external_non_to,
                              "resid_internal_to":resid_internal_to, "resid_internal_non_to":resid_internal_non_to,
                              "resid_external_to":resid_external_to, "resid_external_non_to":resid_external_non_to,
                              "p_value_not_corrected":pval,
                                      "statistic":statistic, "OR":OR, "nuc":new_nuc, "cutoff":cutoff, "mutation_count":mutation_count},
                                     ignore_index=True)

            if family in tables_all[new_nuc]:
                tables_all[new_nuc][family].append(con)
            else:
                tables_all[new_nuc][family] = [con]


    all_nucs_df = pd.DataFrame()

    for n in all_nucs:
        new_nuc = n
        if new_nuc == "T":
            new_nuc = "U"
        all_nucs[n]["p_value"] = multipletests(all_nucs[n]["p_value_not_corrected"], method="fdr_bh")[1]
        all_nucs_df = all_nucs_df.append(all_nucs[n], ignore_index=True)


    all_output = "/sternadi/home/volume3/taliakustin/phyVirus_analysis/contigency_table_significant_incomplete_selection_results.11.02.2020.csv"

    all_nucs_df.to_csv(all_output)
    HM_all = pd.DataFrame()

    for n in tables_all:
        res_families = pd.DataFrame()
        for family in tables_all[n].keys():
            st = sm.stats.StratifiedTable(tables_all[n][family])
            res = parse_summary_MH_test(st.summary())
            res["family"] = family
            res["nuc"] = n
            res["baltimore"] = get_baltimore_classifiaction(family)

            res_families = res_families.append(res, ignore_index=True)
        res_families["p_value"] = multipletests(res_families["OR_1_pvalue"], method="fdr_bh")[1]
        HM_all = HM_all.append(res_families, ignore_index=True)



    all_HM_output = "/sternadi/home/volume3/taliakustin/phyVirus_analysis/HM_significant_incomplete_selection_results.11.02.2020.csv"
    HM_all.to_csv(all_HM_output)

    print(f"no summaries for {count} basenamse")

if __name__ == "__main__":
    main()

