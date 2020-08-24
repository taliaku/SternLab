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


    files = glob.glob(f"/sternadi/home/volume3/taliakustin/phyVirus_analysis/final_dataset/fasta/*fasta")
    df = pd.DataFrame()
    gi_all_file = "/sternadi/home/volume3/taliakustin/phyVirus_analysis/gi_info_edited.csv"
    gi_flu_file = "/sternadi/home/volume3/taliakustin/phyVirus_analysis/Orthomyxo_gi_info_final.csv"

    gi_all = pd.read_csv(gi_all_file)
    gi_flu = pd.read_csv(gi_flu_file)

    res = []
    count_not_triplets = 0
    index = 1
    for f in (files):
        count = 0
        print(f)
        if "Pox" in f or "Herpes" in f or "Reo" in f or "HIV" in f or "SIV" in f or "Nairo" in f or "Hepe" in f:
            continue
        seqs = list(SeqIO.parse(f, "fasta"))
        for s in seqs:
            id = s.id
            if len(s.seq) % 3 != 0:
                count_not_triplets += 1
                continue
            #if s.seq[:3] != "ATG":
            #    count_not_triplets += 1
            #    continue
            if  "Orthomyxo" in f:
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
            s.description = ""
            df = df.append({"seq": s.id, "index":str(index)}, ignore_index=True)
            s.id = str(index)
            index += 1
            res.append(s)

    SeqIO.write(res, "/sternadi/home/volume3/taliakustin/phyVirus_analysis/final_dataset/coding_seqs_human/all_human_seqs.fasta", "fasta")
    df.to_csv("/sternadi/home/volume3/taliakustin/phyVirus_analysis/final_dataset/coding_seqs_human/index_Human.csv")
    print(f"number of seqs not in triplets {count_not_triplets}")
    print("Done")




if __name__ == "__main__":
    main()

