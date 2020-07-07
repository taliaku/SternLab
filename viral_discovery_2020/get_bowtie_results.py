#! /powerapps/share/python-anaconda-3.6/bin/python

import sys
sys.path.insert(0,'/sternadi/home/volume1/taliakustin/SternLab')
from optparse import OptionParser
import os
import glob
from file_utilities import check_filename, check_dirname
from bowtie2_utilities import read_bowtie2_summary
import pandas as pd
def main():
    parser = OptionParser("usage: %prog [options]")
    parser.add_option("-d", "--d", dest="directory", help="directory with bowtie result files")
    parser.add_option("-n", "--name", dest="name", help="analysis name")
    parser.add_option("-o", "--output", dest="output", help="output file")



    (options, args) = parser.parse_args()
    directory = options.directory
    name = options.name
    output = options.output
    output = check_filename(output, Truefile=False)
    df = pd.DataFrame()

    directory = check_dirname(directory)
    files = glob.glob(f"{directory}/*power8.tau.ac.il.OU")
    for f in files:
        print(f)
        results = read_bowtie2_summary(f, viral_discovery=True)
        results["name"] = name
        df = df.append(results, ignore_index=True)
    df.to_csv(output)
    print(f"{output}")
if __name__ == "__main__":
    main()

