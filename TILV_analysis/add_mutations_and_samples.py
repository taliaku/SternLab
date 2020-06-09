#! /usr/local/python_anaconda/bin/python3.4
import sys
from optparse import OptionParser
sys.path.insert(0, "/sternadi/home/volume1/taliakustin/SternLab/")

from file_utilities import check_filename, check_dirname
import pandas as pd
import os
import glob
from Bio.Seq import Seq


def main():
    parser = OptionParser("usage: %prog [options]")
    parser.add_option("-f", "--freqs", dest="freqs_file", help="freqs_file")
    parser.add_option("-t", "--type", dest="type", help="type of filenamen\n1 - sample_X_segment_X\n"
                                                        "2 - PX-SX")

    (options, args) = parser.parse_args()

    freqs_file = options.freqs_file
    output = freqs_file.split(".freqs")[0] + "_with_mutations_and_segment.csv"
    type = options.type

    freqs = pd.read_csv(freqs_file, sep="\t")
    if type == "1":
        sample = int(freqs_file.split("sample_")[1].split("_")[0].split(".")[0])
        segment = int(freqs_file.split("segment_")[1].split("_")[0].split(".")[0])
    elif type == "2":
        sample = int(freqs_file.split("/")[-1].split("P")[1].split("-")[0])
        segment = int(freqs_file.split("/")[-1].split(".")[0].split("S")[1])

    freqs["sample"] = sample
    freqs["segment"] = segment
    print(output)
    add_mutation_to_freq_file(output, freqs=freqs)



def add_mutation_to_freq_file(output, freqs_file = None, freqs = None):
    # assumes that position 1 is the beginning of the CDS
    # removes positions that at the beginning or at the end that are not part of a full codon
    if freqs_file == None and type(freqs) == "NoneType":
        raise Exception("Need to specify or freqs file path or a freqs pandas object")
    elif freqs_file != None and freqs != None:
        print(freqs_file, freqs)
        print(type(freqs))
        raise Exception("Need to specify or freqs file path OR a freqs pandas object - only one!")
    elif freqs_file != None:
        freqs = pd.read_csv(freqs_file, sep="\t")
    freqs = freqs[freqs.Pos % 1 == 0] #removes insertions
    freqs = freqs[freqs.Base != "-"] #removes deletions
    freqs = freqs[freqs.Ref != "-"] #removes deletions

    freqs.reset_index(drop=True, inplace=True)

    first_pos = int(freqs.loc[1].Pos) #gets the first position in the right frameshift
    if first_pos == 1:
        start_from = first_pos
    elif first_pos % 3 == 1:
        start_from = first_pos
    elif first_pos % 3 == 2:
        start_from = first_pos + 2
    elif first_pos % 3 == 0:
        start_from = first_pos + 1

    freqs["Mutation_type"] = None
    freqs["wt_aa"] = None
    freqs["mut_aa"] = None
    freqs["wt_codon"] = None
    freqs["mut_codon"] = None
    freqs["Mutation"] = None

    for pos in range(start_from, int(max(freqs.Pos)), 3): # add mutation information
        temp = freqs.loc[freqs['Pos'].isin([pos, pos+1, pos+2])]
        if len(temp) != 12: #12 - is 4 * 3 [A, C, G, T]

            continue
        first = temp.iloc[0].Ref
        second = temp.iloc[4].Ref
        third = temp.iloc[8].Ref
        wt_codon = "".join([first, second, third])
        wt_aa = str(Seq(wt_codon).translate())

        pos = temp.iloc[0].Pos
        for n in range(0, 12):
            ref_base = temp.iloc[n].Ref
            mut_base = temp.iloc[n].Base
            if n <= 3:
                mut_codon = "".join([mut_base, second, third])
                current_pos = pos
            elif n > 3 and n <= 7:
                mut_codon = "".join([first, mut_base, third])
                current_pos = pos + 1
            elif n > 7  :
                mut_codon = "".join([first, second, mut_base])
                current_pos = pos + 2

            mut_aa = str(Seq(mut_codon).translate())

            if wt_codon == mut_codon:
                mutation_type = "consensus"
            elif wt_aa == mut_aa:
                mutation_type = "synonymous"
            elif wt_aa != "*" and mut_aa == "*":
                mutation_type = "stop"
            else:

                mutation_type = "missense"

            freqs.loc[(freqs["Pos"] == current_pos) & (freqs["Base"] == mut_base), "Mutation_type"] = mutation_type
            freqs.loc[(freqs["Pos"] == current_pos) & (freqs["Base"] == mut_base), "wt_aa"] = wt_aa
            freqs.loc[(freqs["Pos"] == current_pos) & (freqs["Base"] == mut_base), "mut_aa"] = mut_aa
            freqs.loc[(freqs["Pos"] == current_pos) & (freqs["Base"] == mut_base), "wt_codon"] = wt_codon
            freqs.loc[(freqs["Pos"] == current_pos) & (freqs["Base"] == mut_base), "mut_codon"] = mut_codon
            freqs.loc[(freqs["Pos"] == current_pos) & (freqs["Base"] == mut_base), "Mutation"] = ref_base + mut_base

    freqs = freqs[freqs.Mutation_type.notnull()] #removes Nones - rows at the beginning and the end
    freqs.to_csv(output, index=False)
    return freqs


if __name__ == "__main__":
    main()