#! /powerapps/share/python-anaconda-3.6/bin/python

import sys
sys.path.insert(0,'/sternadi/home/volume1/taliakustin/SternLab')
import tqdm
import argparse
import os
import glob
from Bio import SeqIO
import pandas as pd
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--family", type=str,
                        help="family", required=True)

    args = parser.parse_args()
    family = args.family

    replace_dic = {"Arena":"Ar", "Filo":"Fi", "Flavi":"Fl", "Hanta":"Ha", "Orthomyxo":"Or", "Paramyxo":"Pa",
                   "Peribunya":"Pe", "Phenui":"Ph", "Rhabdo":"Rh", "Calci":"Ca", "Corona":"Co", "Picorna":"Pi", "Toga":"To"}

    files = glob.glob(f"/sternadi/home/volume3/taliakustin/phyVirus_analysis/fasta/{family}*fasta")


    coding_seqs_by_host = {}

    gi_all_file = "/sternadi/home/volume3/taliakustin/phyVirus_analysis/gi_info_edited.csv"
    gi_flu_file = "/sternadi/home/volume3/taliakustin/phyVirus_analysis/Orthomyxo_gi_info_final.csv"

    gi_all = pd.read_csv(gi_all_file)
    gi_flu = pd.read_csv(gi_flu_file)

    count_not_triplets = 0

    for f in (files):
        count = 0
        print(f)
        if "Pox" in f or "Herpes" in f or "Reo" in f or "HIV" in f or "SIV" in f or "Nairo" in f:
            continue
        seqs = list(SeqIO.parse(f, "fasta"))
        for s in seqs:
            if len(s.seq) % 3 != 0:
                count_not_triplets += 1
                continue
            if s.seq[:3] != "ATG":
                count_not_triplets += 1
                continue
            id = s.id
            if family == "Orthomyxo":
                seq_row = gi_flu[gi_flu["seq"] == id]

                if seq_row.empty:
                    print(f"NO ID {id}")
                    count += 1
                    continue
                host = seq_row["class"].values[0]
            else:
                seq_row = gi_all[gi_all["seq"] == id]
                if seq_row.empty:
                    print(f"NO ID {id}")
                host = seq_row["host_class"].values[0]
            if host != "Human":
                continue
            if not host in coding_seqs_by_host.keys():
                coding_seqs_by_host[host] = []

            coding_seqs_by_host[host].append(s)
        if count > 0:
            print(f"{count} uknown seqs for {f}")
    output_dir = "/sternadi/home/volume3/taliakustin/phyVirus_analysis/coding_seqs_by_host/"
    for host in coding_seqs_by_host.keys():
        output = f"{output_dir}/{family}_{host}.fasta"
        SeqIO.write(coding_seqs_by_host[host], output, "fasta")
    print(f"number of seqs not in triplets {count_not_triplets}")
    print("Done")




if __name__ == "__main__":
    main()

