#! /powerapps/share/python-anaconda-3.6/bin/python

import sys
sys.path.insert(0,'/sternadi/home/volume1/taliakustin/SternLab')
from optparse import OptionParser
import os
import glob
import pandas as pd
def main():
    parser = OptionParser("usage: %prog [options]")
    parser.add_option("-f", "--file", dest="file")


    (options, args) = parser.parse_args()
    file = options.file

    output = file.split(".csv")[0] + "_with_nucs.csv"
    if os.path.isfile(output):
        print(f"output exists: {output}")
        #return

    pep = pd.read_csv(file)
    pep["A_nuc"] = pep.apply(return_counts, args="A", axis=1)
    pep["C_nuc"] = pep.apply(return_counts, args="C", axis=1)
    pep["G_nuc"] = pep.apply(return_counts, args="G", axis=1)
    pep["U_nuc"] = pep.apply(return_counts, args="U", axis=1)



    pep.to_csv(output)


def return_counts(x, nuc):
    aas_nucs = {"F":{"A":0, "C":0, "G":0, "U":2}, "L":{"A":0, "C":0, "G":0, "U":1}, "I":{"A":1, "C":0, "G":0, "U":1},
            "M":{"A":1, "C":0, "G":1, "U":1}, "V":{"A":0, "C":0, "G":1, "U":1}, "S":{"A":0, "C":1, "G":0, "U":1},
            "P":{"A":0, "C":2, "G":0, "U":0}, "T":{"A":1, "C":1, "G":0, "U":0}, "A":{"A":0, "C":1, "G":1, "U":0},
            "Y":{"A":1, "C":0, "G":0, "U":1}, "H":{"A":1, "C":1, "G":0, "U":0}, "Q":{"A":1, "C":1, "G":0, "U":0},
            "N":{"A":2, "C":0, "G":0, "U":0}, "K":{"A":2, "C":0, "G":0, "U":0}, "D":{"A":1, "C":0, "G":1, "U":0},
            "E":{"A":1, "C":0, "G":1, "U":0}, "C":{"A":0, "C":0, "G":1, "U":1}, "W":{"A":0, "C":0, "G":2, "U":1},
            "R":{"A":0, "C":0, "G":1, "U":0}, "G":{"A":0, "C":0, "G":2, "U":0}}
    count = 0
    for i in aas_nucs.keys():
        if x[i] == None:
            continue
        count += x[i] * aas_nucs[i][nuc]
    return count









if __name__ == "__main__":
    main()

