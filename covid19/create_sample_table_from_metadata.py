#! /powerapps/share/python-anaconda-3.6/bin/python

import sys
sys.path.insert(0,'/sternadi/home/volume1/taliakustin/SternLab')
from optparse import OptionParser
import os
import glob
import tqdm
import pandas as pd
import numpy as np

def main():
    parser = OptionParser("usage: %prog [options]")
    parser.add_option("-m", "--m", dest="metadata", help="metadata file")
    parser.add_option("-i", "--i", dest="include", help="file containing sequences that are included in the analysis")
    parser.add_option("-c", "--c", dest="country", help="country to filter (default Israel)", default="Israel")
    parser.add_option("-a", "--a", dest="authors", help="authors to filter (default Stern Lab)", default="Stern Lab")


    (options, args) = parser.parse_args()
    metadata_file = options.metadata
    include_file = options.include
    country = options.country
    authors = options.authors
    metadata = pd.read_csv(metadata_file, sep="\t")

    with open(include_file, "r")as handle:
        include = handle.read().split("\n")[:-1]

    #metadata_filtered = metadata.loc[(metadata["country"] == country) & (metadata["authors"] == authors) & (metadata["strain"].isin(include))]
    metadata_filtered = metadata.loc[(metadata["country"] == country) & (metadata["authors"] == authors)]

    print(metadata_filtered.groupby(['originating_lab'])["originating_lab"].count())




    #AGES
    ages = pd.DataFrame()
    no_age = metadata_filtered.loc[metadata_filtered["age"]=="?"].shape[0]
    ages = ages.append({"age":"unknown", "sample_num":no_age}, ignore_index=True)

    known_age = metadata_filtered.loc[metadata_filtered["age"] != "?"]
    known_age["age"] = pd.to_numeric(known_age["age"])

    for i in range(0, 90, 10):
        ages = ages.append({"age": f"{i}-{i+9}",
                            "sample_num": known_age.loc[(known_age["age"] >= i) & (known_age["age"] <= i+9)].shape[0]},
                           ignore_index=True)
    ages = ages.append({"age": f"{i} and up",
                        "sample_num": known_age.loc[(known_age["age"] >= i+10)].shape[0]},
                       ignore_index=True)

    print(ages)
    print(metadata_filtered.groupby(['sex'])["sex"].count())



if __name__ == "__main__":
    main()

