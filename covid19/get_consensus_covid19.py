#! /powerapps/share/python-anaconda-3.6/bin/python

import sys
sys.path.insert(0,'/sternadi/home/volume1/taliakustin/SternLab')
from optparse import OptionParser
import os
import glob
import tqdm
import pandas as pd
def main():
    parser = OptionParser("usage: %prog [options]")
    parser.add_option("-d", "--d", dest="directory", help="directory with freqs files")


    (options, args) = parser.parse_args()
    directory = options.directory
    freqs_files = glob.glob(f"{directory}/*freqs.csv")
    con_all = ""
    mutations = pd.DataFrame()
    missing_poss = pd.DataFrame()
    end_pos = 29892

    for f in tqdm.tqdm(freqs_files):
        if "all.freqs.csv" in f:
            continue
        df = pd.read_csv(f)
        df_rank0 = df[(df["rank"] == 0) & (df["ref_base"] != "-") & (df["base"] != "-") & (df["ref_position"] >= 55) & (
                    df["ref_position"] <= 29836)]
        loc = int(df_rank0.head(1)["ref_position"]) - 1
        con = (count) * "N"
        sample = freq.split("/")[-1].split("_")[0]
        for index, row in df_rank0.iterrows():
            new_loc = int(row["ref_position"])
            coverage = int(row["coverage"])
            ref_base = str(row["ref_base"])
            base = str(row["base"])
            if new_loc - loc != 1:
                con += "N" * (new_loc - loc - 1)
                missing_poss = missing_poss.append({"start": loc, "end": new_loc, "sample": sample}, ignore_index=True)
            if coverage < 5:
                if ref_base == base:
                    con += base
                else:
                    con += "N"
            else:
                con += base
                if ref_base != base:
                    mutations = mutations.append(
                        {"ref_position": loc, "ref_base": ref_base, "base": base, "coverage": coverage,
                         "sample": sample}, ignore_index=True)
            loc = new_loc
        con += (end_pos - loc - 1) * "N"
        con_all += f">Israel/{sample}/2020\n{con}\n"

    with open(f"{directory}/consensus_all.fasta", "w") as handle:
            handle.write(con_all)
    missing_poss.to_csv(f"{directory}/missing_positions_all.csv")
    mutations.to_csv(f"{directory}/mutations_all.csv")



if __name__ == "__main__":
    main()

