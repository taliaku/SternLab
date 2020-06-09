#! /powerapps/share/python-anaconda-3.6/bin/python

import sys
sys.path.insert(0,'/sternadi/home/volume1/taliakustin/SternLab')
from optparse import OptionParser
import os
import glob
import seqFileAnalyzer
from file_utilities import check_filename, check_dirname

DIRSEL_PATH="/sternadi/home/volume1/taliakustin/software/phylogenyCode/programs/directionalSelection/directionalSelection"

def main():
    parser = OptionParser("usage: %prog [options]")
    parser.add_option("-f", "--fasta", dest="fasta_file", help="fasta file")
    parser.add_option("-o", "--output", dest="output_dir", help="output directory")

    (options, args) = parser.parse_args()
    fasta = options.fasta_file
    output_dir = options.output_dir
    fasta = check_filename(fasta)
    output_dir = check_dirname(output_dir)
    seqFileAnalyzer.analyze_nuc_frequencies_from_fasta(fasta, output_dir=output_dir)



if __name__ == "__main__":
    main()

