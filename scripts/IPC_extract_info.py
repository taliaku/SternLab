#! /powerapps/share/python-anaconda-3.6/bin/python

from optparse import OptionParser
import glob
import os
import numpy as np
import pandas as pd
def main():
    parser = OptionParser("usage: %prog [options]")
    parser.add_option("-d", "--dir", dest="dir", help="directory")
    parser.add_option("-o", "--output", dest="output_file", help="tree output file")
    (options, args) = parser.parse_args()
    dir = options.dir
    output = options.output_file
    files = glob.glob(f"{dir}/*pI.txt")
    res = {"base":[], "seq":[], "averageIP":[]}
    for f in files:
       base = os.path.basename(f).split(".aln")[0].split(".fasta")[0].split(".")[0]
       with open(f, "r") as handle:
        data = handle.readlines()
        data_len = len(data)
        l = 0
        while l < data_len:
            if ">" not in data[l]:
                l += 1
            else:
                seq = data[l].split(">")[1].split(" weight")[0]
                mean = np.mean([float(i) for i in data[l+1].strip().split(",")])
                res["base"].append(base)
                res["seq"].append(seq)
                res["averageIP"].append(mean)
                l += 2

    res=pd.DataFrame(res)
    res.to_csv(output)





if __name__ == "__main__":
    main()
