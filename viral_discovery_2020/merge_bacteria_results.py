#! /powerapps/share/python-anaconda-3.6/bin/python

import sys
sys.path.insert(0,'/sternadi/home/volume1/taliakustin/SternLab')
from optparse import OptionParser
import os
import glob
from file_utilities import check_filename, check_dirname
from bowtie2_utilities import read_bowtie2_summary
import pandas as pd
import tqdm
def main():
    parser = OptionParser("usage: %prog [options]")
    parser.add_option("-d", "--d", dest="directory", help="directory with mapping of bacteria file")

    tax_names_file = "/sternadi/datasets/volume1/bacteriaDB/hmpdacc/bacteria_DB_no_viruses_correct.unique.csv"
    tax_names = pd.read_csv(tax_names_file)


    (options, args) = parser.parse_args()
    directory = options.directory
    files = glob.glob(f"{directory}/*all_taxon_count.txt")
    all = pd.DataFrame()
    for f in tqdm.tqdm(files):

        base = f.split("/")[-1].split(".mapped_bacteria_all_taxon_count.txt")[0]
        patient = base.split("_")[0]
        gen_mat = base.split("_")[1]
        lane = int(base.split("_")[3].split("L00")[1].split(".")[0])

        with open(f, "r") as handle:
            counts = handle.readlines()
        counts = [{"count":int(i.strip().split(" ")[0]), "tax":i.strip().split(" ")[-1].split("\n")[0]} for i in counts]
        counts_df = pd.DataFrame(counts)
        counts_df["patient"] = patient
        counts_df["gen_mat"] = gen_mat
        counts_df["lane"] = lane

        merged = counts_df.merge(tax_names, on="tax")
        merged["count"] = merged["count"] / 2



        all = all.append(merged, ignore_index=True)


    all.to_csv("/sternadi/nobackup/volume1/talia_temp/viral_discovery2/mapped_to_bacteria_filtered_20/all_mapped_to_bacteria_with_taxname.csv")

if __name__ == "__main__":
    main()

