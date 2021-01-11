#! /powerapps/share/python-anaconda-3.6/bin/python

import os,sys,inspect
parentdir = os.path.dirname(os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe()))))
sys.path.insert(0,parentdir)
import argparse
from Bio import SeqIO


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--fasta", type=str,
                        help="aligned fasta path", required=True)
    parser.add_argument("-p", "--percent_Ns", type=int, default=5,
                        help="N percent threshold, sequences with more are removed", required=False)

    args = parser.parse_args()
    fasta_file = args.fasta
    percent_Ns = args.percent_Ns

    fasta = list(SeqIO.parse(fasta_file, "fasta"))
    new_fasta = []
    aln_len = len(fasta[0].seq) - fasta[0].seq.count('-')

    for f in fasta:
        gaps = f.seq.count("-")
        if gaps == aln_len:
            continue
        Ns = f.seq.count("N") + f.seq.count("n")
        if gaps + Ns == aln_len: #if all seq is gaps and Ns
            continue
        len_without_gaps = aln_len - gaps
        percent_Ns_seq = (Ns / len_without_gaps) * 100

        if percent_Ns_seq >= percent_Ns:
            continue
        else:
            new_fasta.append(f)
    output_file = fasta_file.replace(".fasta", f".under_{percent_Ns}_percent_Ns.fasta").replace(".aln", f".under_{percent_Ns}_percent_Ns.aln")
    SeqIO.write(new_fasta, output_file, "fasta")

if __name__ == "__main__":
    main()