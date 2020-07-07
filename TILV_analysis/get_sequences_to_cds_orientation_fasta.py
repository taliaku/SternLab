#! /powerapps/share/python-anaconda-3.6/bin/python

import sys
sys.path.insert(0,'/sternadi/home/volume1/taliakustin/SternLab')
from optparse import OptionParser
import os
import glob
import argparse
from Bio import SeqIO, Seq, SeqRecord
import re
from file_utilities import check_filename, check_dirname
import seqFileAnalyzer
import pandas as pd

def main(args):
    fasta_file = check_filename(args.fasta)
    gb_file = check_filename(args.gb)
    fasta = list(SeqIO.parse(fasta_file, "fasta"))
    output = check_filename(args.output, Truefile=False)
    gb =  list(SeqIO.parse(gb_file, "fasta"))
    base = fasta_file.split("/")[-1].split(".fasta")[0]
    new_fasta_file = []
    index_text = ""
    ids = []
    count = 0
    didnt_work = 0
    ids_fasta = [s.id for s in fasta]
    ids_gb = [s.id for s in gb]

    for f in fasta:
        #print(f.description)
        try:
            index = ids_gb.index(f.description)
        except:
            part_id = "|".join(f.description.split("|")[:5])
            for index in range(len(ids_gb)):
                if part_id in ids_gb[index]:
                    break


        if gb[index].seq in f.seq:
            pass
        elif gb[index].seq.reverse_complement() in f.seq:
            f.seq = f.seq.reverse_complement()
            count += 1
        else:
            found = False
            if (ids_gb.count(f.description)) > 1: #more than 1 cds record for this id
                cds_ids = [i for i in range(len(ids_gb)) if ids_gb[i] == f.id]
                for i in cds_ids:
                    if gb[i].seq in f.seq:
                        found = True
                        break
                    elif gb[i].seq.reverse_complement() in f.seq:
                        f.seq = f.seq.reverse_complement()
                        found = True
                        break
            if not found:
                print(f.id)
                print(gb[index].id)
                print(f.seq)
                print(gb[index].seq)
                didnt_work += 1

        new_id = f.id.split("|")[1] + f"_{base}"
        if new_id in ids:
            print("PROBLEM!!!")
        ids.append(new_id)
        index_text += f"{new_id} {f.description}\n"
        f.id = new_id
        f.description = ""
        new_fasta_file.append(f)
    print(f"reverse complemented {count} sequences (out of {len(fasta)} seqneces)")
    print(f"{didnt_work} didn't work")
    SeqIO.write(new_fasta_file, output, "fasta")
    output_index = output.split(".fasta")[0] + ".index.txt"
    with open(output_index, "w") as handle:
        handle.write(index_text)






if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--fasta", type=str,
                        help="fasta file", required=True)
    parser.add_argument("-g", "--gb", type=str,
                        help="gene bank file", required=True)
    parser.add_argument("-o", "--output", type=str,
                        help="output fasta file", required=True)
    args = parser.parse_args()
    main(args)

