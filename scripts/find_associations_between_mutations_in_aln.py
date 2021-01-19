#! /usr/local/python_anaconda/bin/python3.4
import pandas as pd
import matplotlib
matplotlib.use('agg')
import seaborn as sns
import matplotlib.pyplot as plt
import os,sys,inspect
parentdir = os.path.dirname(os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe()))))
sys.path.insert(0,parentdir)
import argparse
from file_utilities import check_filename
from Bio import SeqIO, AlignIO
from scipy.stats import fisher_exact, chi2_contingency


def find_association_between_two_mutations(aln, p1, p2, v1, v2):
    aln_p12 =  aln[:, p1-1:p1] + aln[:, p2-1:p2]
    aln_p12_list = [str(s.seq) for s in aln_p12]
    both = aln_p12_list.count(v1+v2)
    first = [str(s.seq) for s in aln[:,p1-1:p1]].count(v1)
    second = [str(s.seq) for s in aln[:,p2-1:p2]].count(v2)
    none = len(aln) - both - first - second

    if both == 0:
        return None
    g, p, dof, expctd = chi2_contingency([[none, first], [second, both]])
    return(p, none, first, second, both)


def list_int(values):
    return [int(i) for i in values.split(',')]

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", dest="input", help="input alignment path")
    parser.add_argument("-o", "--output", dest="output", help="output csv path")
    parser.add_argument("-p", "--positions_to_check", dest="positions", help="positions to check", type=list_int, default="")
    parser.add_argument("-w", "--WT", dest="WT", help="WT sequence")

    args = parser.parse_args()


    alignment_file = args.input
    output = args.output
    alignment_file = check_filename(alignment_file)
    output = check_filename(output, Truefile=False)
    positions = args.positions
    positions.sort()
    WT = args.WT

    for s in SeqIO.parse(alignment_file, "fasta"):
        if s.id == WT:
            WT_seq = s.seq
    print(WT_seq)

    aln = AlignIO.read(alignment_file, "fasta")
    res = {"p1":[], "p2":[], "v1":[], "v2":[], "none":[], "first":[], "second":[], "both":[],  "pvalue":[]}
    for p1 in positions:
        for p2 in positions:
            if  p2 <= p1:
                continue
            print(p1, p2)
            p1_unique = list(set(aln[:,p1-1]))
            p2_unique =list(set(aln[:,p2-1]))
            p1_WT = WT_seq[p1-1]
            p2_WT = WT_seq[p2 - 1]
            for v1 in p1_unique:
                if v1 == p1_WT:
                    continue
                for v2 in p2_unique:
                    if v2 == p2_WT:
                        continue


                    res_temp = find_association_between_two_mutations(aln, p1, p2, v1, v2)
                    if res_temp == None:
                        continue
                    p, none, first, second, both = res_temp
                    res["p1"].append(p1)
                    res["p2"].append(p2)
                    res["v1"].append(v1)
                    res["v2"].append(v2)
                    res["none"].append(none)
                    res["first"].append(first)
                    res["second"].append(second)
                    res["both"].append(both)
                    #res["OR"].append(OR)
                    res["pvalue"].append(p)


    res = pd.DataFrame(res)
    res.to_csv(output)









if __name__ == '__main__':
    main()
