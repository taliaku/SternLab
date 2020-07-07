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
    new_fasta_file = []
    index_text = ""
    ids = []
    count = 0
    no_cds = 0
    ambisense_problems = 0
    for g in SeqIO.parse(gb_file, "gb"):
        for f in fasta:
            if g.id in f.id and g.description in f.description:
                orientation = []
                cdss = []
                for feature in g.features:
                    if feature.type == "source":
                        mol_type = feature.qualifiers["mol_type"][0]
                    if feature.type == "CDS":
                        is_CDS = True
                        orientation.append(feature.location.strand)
                        cdss.append(feature.location.extract(g.seq))
                if len(cdss) == 0:
                    no_cds += 1
                    continue
                if len(set(orientation)) != 1:
                    #print(f"{f.id} is Ambisense")
                    if mol_type == "viral cRNA" or mol_type=="mRNA":
                        new_id = f.id.split(" ")[0]
                        if new_id in ids:
                            print("PROBLEM!!!")
                            continue
                        ids.append(new_id)
                        index_text += f"{new_id} {f.description}\n"
                        f.id = new_id
                        f.description = ""
                        new_fasta_file.append(f)
                    elif mol_type == "genomic RNA" or mol_type=="genomic DNA":
                        f.seq = f.seq.reverse_complement()
                        count += 1
                        new_id = f.id.split(" ")[0]
                        if new_id in ids:
                            print("PROBLEM!!!")
                            continue
                        ids.append(new_id)
                        index_text += f"{new_id} {f.description}\n"
                        f.id = new_id
                        f.description = ""
                        new_fasta_file.append(f)
                    else:
                        print(mol_type)
                        raise TypeError("problem with ambisense")
                        ambisense_problems+=1
                        continue
                else:
                    for cds in cdss:
                        if cds in f.seq:
                            new_id = f.id.split(" ")[0]
                            if new_id in ids:
                                print("PROBLEM!!!")
                                continue
                            ids.append(new_id)
                            index_text += f"{new_id} {f.description}\n"
                            f.id = new_id
                            f.description = ""
                            new_fasta_file.append(f)
                            break
                        elif cds.reverse_complement() in f.seq:
                            f.seq = f.seq.reverse_complement()
                            count += 1

                            new_id = f.id.split(" ")[0]
                            if new_id in ids:
                                print("PROBLEM comp!!!")
                            ids.append(new_id)
                            index_text += f"{new_id} {f.description}\n"
                            f.id = new_id
                            f.description = ""
                            new_fasta_file.append(f)
                            break






    print(f"reverse complemented {count} sequences (out of {len(fasta)} seqneces)")
    print(f"number of sequences in new file {len(new_fasta_file)}")
    print(f"sequences with no CDS {no_cds}")
    print(f"sequences with ambisense problem {ambisense_problems}")
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

