#! /powerapps/share/python-anaconda-3.6/bin/python

import sys
sys.path.insert(0,'/sternadi/home/volume1/taliakustin/SternLab')
from optparse import OptionParser
import os
import glob
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
import re
from file_utilities import check_filename, check_dirname
import seqFileAnalyzer
import pandas as pd

def main(args):
    genomic_dir = "/sternadi/home/volume3/taliakustin/TiLV2020/genome/"
    cds_dir = "/sternadi/home/volume3/taliakustin/TiLV2020/cds/"
    non_coding_dir = "/sternadi/home/volume3/taliakustin/TiLV2020/non_coding/"
    genomis_files = glob.glob(f"{genomic_dir}/segment*.fasta")
    cds_files = glob.glob(f"{cds_dir}/segment*.fasta")
    for g in genomis_files:
        base = g.split("/")[-1].split(".fasta")[0]
        c = glob.glob(f"{cds_dir}{base}_CDS.fasta")[0]
        genome = list(SeqIO.parse(g, "fasta"))
        cds = list(SeqIO.parse(c, "fasta"))
        print(c)
        uncoding = []
        for i in range(len(genome)):
            if cds[i].id not in genome[i].id:
                raise(TypeError("Not the same order of genome and CDS"))
            temp_seq = genome[i]
            temp_seq.seq = Seq(str(genome[i].seq).replace(str(cds[i].seq), ""))
            if temp_seq.seq != "":

                uncoding.append(temp_seq)
        if uncoding != []:
            output = f"{non_coding_dir}{base}_NON_CODING.fasta"
            SeqIO.write(uncoding, output, "fasta")






if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    #parser.add_argument("-g", "--genomic", type=str, help="directory with genomic alignments ", required=True)
    args = parser.parse_args()
    main(args)

