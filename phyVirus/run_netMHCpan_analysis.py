#! /powerapps/share/python-anaconda-3.6/bin/python

import sys
sys.path.insert(0,'/sternadi/home/volume1/taliakustin/SternLab')
from optparse import OptionParser
import os
import glob
from file_utilities import check_filename, check_dirname
from dirSel_utilities import write_params_file
import pbs_runners
from netMHCpan_utilities import   go_over_netMHCpan_results, make_peptide_csv_for_each_seq, \
    merge_peptide_file_with_netMHCpan4_results, merge_peptide_file_with_netMHCpan4_results_get_null_results, go_over_MHC_results_by_netMHC_file, go_over_MHC_results_by_fasta_file, go_over_MHC_results_by_fasta_with_set

def main():
    parser = OptionParser("usage: %prog [options]")
    parser.add_option("-f", "--file", dest="file", help="fasta file")
    (options, args) = parser.parse_args()
    file = options.file

    go_over_MHC_results_by_fasta_with_set(file)

if __name__ == "__main__":
    main()

