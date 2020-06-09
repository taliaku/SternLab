#! /powerapps/share/python-anaconda-3.6/bin/python

import sys
sys.path.insert(0,'/sternadi/home/volume1/taliakustin/SternLab')
from optparse import OptionParser
import os
import glob
import pandas as pd
import tqdm
def main():
    parser = OptionParser("usage: %prog [options]")
    parser.add_option("-f", "--file", dest="file")

    (options, args) = parser.parse_args()
    file = options.file
    if "SIV" in file or "check" in file:
        print("SIV")
        return
    run(file)

def run(file):
    base = file.split("/")[-1].split("-translated")[0]
    directory = f"/sternadi/nobackup/volume1/talia_temp/netMHC_by_seq_250/{base}"
    print(directory)
    if not os.path.isdir(directory):
        print("directory doen't exist")
        return
    results = glob.glob(f"{directory}/*")
    print(len(results))
    output = f"/sternadi/nobackup/volume1/talia_temp/netMHC_by_seq_250/{base}.all_alleles.csv"
    if os.path.isfile(output):
        print("exists")
        return
    if len(results) < 248:
        print("not finished")
        return
    res_all = pd.DataFrame()
    for res_file in (results):
        if os.stat(res_file).st_size == 0:
            continue
        print(res_file)
        allele = res_file.split("/")[-1].split(".")[1]
        try:
            res = pd.read_csv(res_file, skiprows=1, sep="\t")
        except pd.io.common.EmptyDataError:
            continue
        res = res.loc[res.NB == 1]
        res["allele"] = allele
        res['property'] = res['Ave'].apply(lambda x: 'strong' if x > 0.5 else 'weak')
        res_all = res_all.append(res, ignore_index=True)
    res_all.to_csv(output)
    







if __name__ == "__main__":
    main()

