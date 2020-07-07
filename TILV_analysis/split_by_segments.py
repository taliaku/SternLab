#! /powerapps/share/python-anaconda-3.6/bin/python

import sys
sys.path.insert(0,'/sternadi/home/volume2/noam/SternLab')
from optparse import OptionParser
import os
import glob
import argparse
from Bio import SeqIO
import re
from file_utilities import check_filename


def main(args):
    file = check_filename(args.file)
    cds = args.cds
    outdir = os.path.dirname(file)
    fasta = SeqIO.parse(file, "fasta")
    r = re.compile("segment ([A-Za-z0-9]+)")
    segments = {}

    for f in fasta:
        s = []
        #print(r.findall(f.description))
        try:
            s = r.findall(f.description)[0]
            if s not in segments.keys():
                segments[s] = []
            segments[s].append(f)
        except:
            print(f.description)
            print(r.findall(f.description))


    for i in segments.keys():
        if cds:
            output_file = f"{outdir}/segment_{i}_CDS.fasta"
        else:
            output_file = f"{outdir}/segment_{i}.fasta"
        SeqIO.write(segments[i], output_file, "fasta")



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--file", type=str,
                        help="input file TiLV sequences", required=True)
    parser.add_argument("-c", "--cds", action="store_true")
    args = parser.parse_args()
    main(args)

