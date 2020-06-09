#! /powerapps/share/python-anaconda-3.6/bin/python

import sys
sys.path.insert(0,'/sternadi/home/volume1/taliakustin/SternLab')
from optparse import OptionParser
import os
import glob
from file_utilities import check_filename, check_dirname
from phyVirus_utilities import get_basename, get_baltimore_classifiaction
import fastml_utilities
from Bio import SeqIO
import pandas as pd

def main():
    parser = OptionParser("usage: %prog [options]")
    parser.add_option("-f", "--file", dest="fasta", help="fasta file")

    (options, args) = parser.parse_args()

    fasta = options.fasta
    fasta = check_filename(fasta)

    seqs = list(SeqIO.parse(fasta, "fasta"))
    family = fasta.split("/")[-1].split(".fasta")[0]
    baltimore = get_baltimore_classifiaction(family)
    family_files = glob.glob("/sternadi/home/volume3/taliakustin/phyVirus_analysis/fasta/{}*".format(family))
    df = pd.DataFrame()
    for s in seqs:
        num_of_matches = 0
        for f in family_files:
            family_seqs = list(SeqIO.parse(f, "fasta"))
            base = get_basename(f)
            for fs in family_seqs:
                if s.seq == fs.seq:
                    df = df.append({"family":family, "base":base, "baltimore":baltimore, "seq":fs.description, "description":s.description}, ignore_index=True)
        if num_of_matches > 1:
            print(s.description)
    df.to_csv("/sternadi/home/volume3/taliakustin/phyVirus_analysis/{}_gi_info.csv".format(family))







if __name__ == "__main__":
    main()

