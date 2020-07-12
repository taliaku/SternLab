#! /powerapps/share/python-anaconda-3.6/bin/python

import sys
sys.path.insert(0,'/sternadi/home/volume1/taliakustin/SternLab')
import argparse
import os
import glob
from Bio import SeqIO
import pandas as pd
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--fasta", type=str,
                        help="fasta", required=True)

    args = parser.parse_args()
    fasta_file = args.fasta

    out_dir = "/sternadi/home/volume3/taliakustin/phyVirus_analysis/human_coding_seqs/"
    base = fasta_file.split("/")[-1].split(".fasta")[0]
    seqs = list(SeqIO.parse(fasta_file, "fasta"))


    gi_all = "/sternadi/home/volume3/taliakustin/phyVirus_analysis/gi_info_edited.csv"
    gi_flu = "/sternadi/home/volume3/taliakustin/phyVirus_analysis/Orthomyxo_gi_info_final.csv"

    orthomyxo = False
    if "Orthomyxo" in fasta_file:
        gi = gi_flu
        orthomyxo = True
    else:
        gi = gi_all
    df = pd.read_csv(gi)

    human_seqs = []
    count = 0
    for s in seqs:
        if len(s.seq) % 3 != 0:
            continue
        id = s.id
        seq_row = df[df["seq"] == id]
        if seq_row.empty:
            print(f"NO ID {id}")
        if orthomyxo:
            host = seq_row["class"].values[0]
        else:
            host = seq_row["host_class"].values[0]
        if host == "Human" or host == "Homo sapiens":
            human_seqs.append(s)
            count += 1
    print(f"in {base}: {count} out of {len(seqs)} are human")
    if count > 0:
        outfile = out_dir + base + "_human_seqs.fasta"
        SeqIO.write(human_seqs, outfile, "fasta")
        print(f"saved file {outfile}")
    else:
        print("no human seqs")


if __name__ == "__main__":
    main()

