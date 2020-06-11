#! /powerapps/share/python-anaconda-3.6/bin/python

import sys
sys.path.insert(0,'/sternadi/home/volume2/noam/SternLab')
from optparse import OptionParser
import os
import glob
import argparse
from Bio import SeqIO
import re
from file_utilities import check_filename, check_dirname
import seqFileAnalyzer
import pandas as pd

def main(args):
    genomic_dir = check_dirname(args.genomic)
    files = glob.glob(f"{genomic_dir}/*.aln.best.fas")
    freqs_df = pd.DataFrame()
    for f in files:
        segment = int(f.split("segment_")[1].split(".")[0])
        freqs = seqFileAnalyzer.base_frequencies(f)
        freqs["segment"] = segment
        freqs_df = freqs_df.append(freqs, ignore_index=True)
    print(freqs_df)



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-g", "--genomic", type=str,
                        help="directory with genomic alignments ", required=True)
    args = parser.parse_args()
    main(args)

