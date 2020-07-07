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
    gb_file = check_filename(args.gb)
    fasta_file = check_filename(args.fasta)
    gb = list(SeqIO.parse(gb_file, "gb"))
    fasta = list(SeqIO.parse(fasta_file, "fasta"))
    print(len(gb), len(fasta))
    new_seqs = []
    j = 0
    for i in range(len(fasta)):
        while not gb[j].id in fasta[i].id:
            j += 1

        s = gb[i]
        for feature in s.features:
            if feature.type == "source":
                mol_type = feature.qualifiers["mol_type"][0]
                if mol_type == "genomic RNA" or mol_type == "genomic DNA":
                    temp_seq = fasta[i]
                    temp_seq.seq = temp_seq.seq.reverse_complement()
                    new_seqs.append(temp_seq)
                elif mol_type == "mRNA" or mol_type == "viral cRNA":
                    new_seqs.append(fasta[i])
                else:
                    print(mol_type)
                    new_seqs.append(fasta[i])
    output = fasta_file.split(".fasta")[0] + "_correctly_oriented.fasta"
    SeqIO.write(new_seqs, output, "fasta")








if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-g", "--gb", type=str, help="gb file", required=True)
    parser.add_argument("-f", "--fasta", type=str, help="fasta file", required=True)
    args = parser.parse_args()
    main(args)

