#! /powerapps/share/python-anaconda-3.6/bin/python

import sys
sys.path.insert(0,'/sternadi/home/volume1/taliakustin/SternLab')
from optparse import OptionParser
import os
import glob
from file_utilities import check_filename, check_dirname
from dirSel_utilities import write_params_file
import pbs_runners
import pandas as pd

def main():
    parser = OptionParser("usage: %prog [options]")
    parser.add_option("-f", "--file", dest="file", help="file")



    (options, args) = parser.parse_args()
    file = options.file
    output = file.split("_results.csv")[0] + "_results_no_duplicates.csv"
    if os.path.isfile(output):
        return
    net = pd.read_csv(file)
    net = net.loc[:,["allele", "base", "basename", "bind", "core", "icore", "peptide", "pos","rank"]]
    net = net.drop_duplicates()
    net.to_csv(output)

if __name__ == "__main__":
    main()

