#! /powerapps/share/python-anaconda-3.6/bin/python

import sys
sys.path.insert(0,'/sternadi/home/volume2/noam/SternLab')
from optparse import OptionParser
import os
import glob
import tqdm
import pandas as pd
import re
from freqs_utilities import estimate_insertion_freq


def main():
    parser = OptionParser("usage: %prog [options]")
    parser.add_option("-d", "--d", dest="directory", help="directory with freqs files")
    parser.add_option('-o', '--o', dest='output_basename')
    parser.add_option('-m', '--mask', action='store_true', default=False, dest='mask')

    (options, args) = parser.parse_args()
    directory = options.directory
    output_basename = options.output_basename
    freqs_files = glob.glob(f"{directory}/*freqs.csv")
    con_all = ""
    mutations = pd.DataFrame()
    mutations_turned_ns = pd.DataFrame()
    missing_poss = pd.DataFrame()
    end_pos = 29903

    for f in tqdm.tqdm(freqs_files):
        #if "all.freqs.csv" in f or 'all.csv' in f or '2089839' not in f:
        if "all.freqs.csv" in f or 'all.csv' in f or 'all_freqs.csv' in f:
            continue
        sample = f.split("/")[-1].split("_")[0]
        df = pd.read_csv(f)
        df = df[(df.coverage > 10)]
        df_snps_and_dels = df[(df["rank"] == 0) & (df["ref_base"] != "-") & (df["ref_position"] >= 55) & (df["ref_position"] <= 29836)]
        df = estimate_insertion_freq(df)
        df_insertions = df[(df["estimated_freq"] >= 0.8) & (df["ref_base"] == "-") & (df["base"] != "-") & (df["ref_position"] >= 55) & (df["ref_position"] <= 29836) & (df['coverage'] > 10)]
        df_changes = pd.concat([df_snps_and_dels, df_insertions], sort=False).sort_values('ref_position')
        count = int(df_changes.head(1)["ref_position"])
        loc = int(df_changes.head(1)["ref_position"]) - 1
        con = (count) * "N"    
        for index, row in df_changes.iterrows():
            new_loc = int(row["ref_position"])
            coverage = int(row["coverage"])
            ref_base = str(row["ref_base"])
            base = str(row["base"])
            frequency = float(row['frequency'])
            
            if ref_base != '-': # not insertion
                if new_loc - loc != 1:
                    con += "N" * (new_loc - loc - 1)
                    missing_poss = missing_poss.append({"start": loc, "end": new_loc, "sample": sample}, ignore_index=True)
                
                if base == ref_base:
                    con += base
                else: # rank0 ref_base != base
                    if coverage < 5 or frequency < 0.8:
                        con += "N"
                        mutations_turned_ns = mutations_turned_ns.append({"ref_position":new_loc, "ref_base":ref_base, "base":base, "coverage":coverage, 'frequency':frequency, "sample":sample}, ignore_index=True)
                    else: # freq >= 0.8 and coverage >= 5
                        con += base
                        mutations = mutations.append({"ref_position":new_loc, "ref_base":ref_base, "base":base, "coverage":coverage, 'frequency':frequency, "sample":sample}, ignore_index=True)
                loc = new_loc
                
            if ref_base == '-' and new_loc - loc < 1: # insertion 
                con += base
                mutations = mutations.append({"ref_position":new_loc, "ref_base":ref_base, "base":base, "coverage":coverage, 'frequency':frequency, "sample":sample}, ignore_index=True)
                
        con += (end_pos - loc - 1) * "N"
        if options.mask:
            con = post_processing(con)
        con_all += f">Israel/{sample}/2021\n{con}\n"

    with open(f"{directory}/{output_basename}_consensus_all.fasta", "w") as handle:
            handle.write(con_all)
    missing_poss.to_csv(f"{directory}/{output_basename}_missing_positions_all.csv")
    mutations.to_csv(f"{directory}/{output_basename}_mutations_all.csv")
    mutations_turned_ns.to_csv(f"{directory}/{output_basename}_mutations_turned_ns_all.csv")

def post_processing(con_seq):
    def repl(m):
        return 'N' * len(m.group())

    return re.sub('N[AGCT-]{1,20}N', repl, (re.sub('N[AGCT-]{1,20}N', repl, con_seq)))

if __name__ == "__main__":
    main()

