#! /powerapps/share/python-anaconda-3.6/bin/python

import sys
sys.path.insert(0,'/sternadi/home/volume1/taliakustin/SternLab')
from optparse import OptionParser
import os
import glob
from file_utilities import check_filename, check_dirname
from phyVirus_utilities import get_basename, get_baltimore_classifiaction
import fastml_utilities
from dirSel_utilities import extract_dirSel_parameters_single_file, extract_dirSel_parameters, merge_alternative_and_null_dfs
from scipy.stats.distributions import chi2
from statsmodels.stats.multitest import multipletests
import pandas as pd



def main():
    files = glob.glob("/sternadi/home/volume3/taliakustin/phyVirus_analysis/fasta/*fasta")
    cutoffs = ["no_cutoff", "0.1", "0.05", "0.01", "0.005", "0.001"]

    null = glob.glob("/sternadi/home/volume3/taliakustin/phyVirus_analysis/beta_null/*results")
    null_results = extract_dirSel_parameters(null, to_return=True, model="null")
    sig_all = pd.DataFrame()
    sig_basenames = []
    for c in cutoffs:
        if c == "no_cutoff":
            cutoffs = glob.glob("/sternadi/home/volume3/taliakustin/phyVirus_analysis/beta_alt/*results")
        else:
            if c == "0.1":
                continue
            cutoffs = glob.glob("/sternadi/home/volume3/taliakustin/phyVirus_analysis/beta_alt_%s/*results" % (c))
        print(c, len(cutoffs))
        cutoff_results = extract_dirSel_parameters(cutoffs, to_return=True, model=c)
        merged = merge_alternative_and_null_dfs(cutoff_results, null_results, to_return=True, files=False, degrees_of_freedom=1)
        if c =="no_cutoff":
            no_cutoff_merged = merged
        #merged.to_csv("/sternadi/home/volume3/taliakustin/phyVirus_analysis/tmp.csv")
        sig =  merged.loc[merged.sig_by_chi2 == True]
        sig = sig[~sig.basename.isin(sig_basenames)]
        sig["cutoff"] = c
        sig_basenames = sig_basenames + list(sig.basename)
        sig_all = sig_all.append(sig, ignore_index=True)
    #get_non_sig
    no_sig = no_cutoff_merged[~no_cutoff_merged.basename.isin(sig_basenames)]
    no_sig["cutoff"] = "no_cutoff"
    sig_all = sig_all.append(no_sig, ignore_index=True)

    sig_all.to_csv("/sternadi/home/volume3/taliakustin/phyVirus_analysis/beta_cutoffs_free_tau_results.csv")





if __name__ == "__main__":
    main()

