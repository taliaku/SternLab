#! /powerapps/share/python-anaconda-3.6/bin/python

import sys
sys.path.insert(0,'/sternadi/home/volume1/taliakustin/SternLab')
from optparse import OptionParser
import os
import glob
from file_utilities import check_filename, check_dirname
from dirSel_utilities import write_params_file
import pbs_runners
from netMHCpan_utilities import netMHCpan_to_csv
from pandas_utilities import merge_dfs
import pandas as pd
import numpy as np
from phylogenetic_utilities import label_trees_by_cutoff
from phyVirus.phyVirus_utilities import  get_codon_alignment_phy_from_basename, get_rooted_tree_from_basename
import PAML_utilities
import fastml_utilities
from phyVirus.get_baltimore import  get_baltimore_classifiaction
from Bio import Entrez

HOSTS = {"Homo sapiens":[], "Mammalia":[], "Aves":[], "Insecta":[], "Actinopterygii":[], "Arachnida":[], "Reptilia":[], "Lepidosauria":[]}

codons = ['CGA', 'CGC', 'CGG', 'CGU', 'AGA', 'AGG', 'CUA', 'CUC', 'CUG', 'CUU', 'UUA', 'UUG', 'UCA', 'UCC', 'UCG',
          'UCU', 'AGC', 'AGU', 'ACA', 'ACC', 'ACG', 'ACU', 'CCA', 'CCC', 'CCG', 'CCU', 'GCA', 'GCC', 'GCG', 'GCU',
          'GGA', 'GGC', 'GGG', 'GGU', 'GUA', 'GUC', 'GUG', 'GUU', 'AAA', 'AAG', 'AAC', 'AAU', 'CAA', 'CAG', 'CAC',
          'CAU', 'GAA', 'GAG', 'GAC', 'GAU', 'UAC', 'UAU', 'UGC', 'UGU', 'UUC', 'UUU', 'AUA', 'AUC', 'AUU', 'AUG',
          'UGG', 'UAA', 'UAG', 'UGA']


def return_nuc_content_from_codon_data(l):
    l = l.split(" ")[:-1]
    nucs = {"A":0, "C":0, "G":0, "U":0}
    for i in range(len(codons)):
        for n in nucs:
            nucs[n] += codons[i].count(n) * int(l[i])
    return(nucs)



def main():
    Entrez.email = "taliakustin@gmail.com"
    files = glob.glob("/sternadi/home/volume3/taliakustin/phyVirus_analysis/hosts/*spsum")
    res_all = pd.DataFrame()
    for f in files:
        res = pd.DataFrame()
        out_file = f.split(".spsum")[0] + ".out.csv"
        if os.path.isfile(out_file):
            res_temp = pd.read_csv(out_file)
            res_all = res_all.append(res_temp, ignore_index=True)
            continue
        print(f"working on {f}")
        data = open(f, "r").read()
        data = data.split("\n")
        for i in range(0,len(data), 2):
            if "mitochondrion" in data[i]:
                continue
            if "chloroplast" in data[i]:
                continue
            if data[i] == "":
                continue
            num_of_seqs = int(data[i].split(" ")[-1])
            if num_of_seqs < 50:
                continue
            tax_id = data[i].split(":")[0]
            name = data[i].split(":")[1]
            print(tax_id)
            try:
                handle = Entrez.efetch(db="taxonomy", id=tax_id, retmode="xml")
                records = Entrez.parse(handle)
                records = list(records)
                if len(records) != 1:
                    print("Error!! more than 1 tax result")
                r = records[0]
                lineage = r["Lineage"]

                for h in HOSTS:
                    if h in lineage:
                        nucs = return_nuc_content_from_codon_data(data[i+1])
                        nucs["host"] = h
                        nucs["tax_id"] = tax_id
                        nucs["cds_num"] = num_of_seqs
                        nucs["name"] = name
                        res = res.append(nucs, ignore_index=True)
                        break
            except:
                nucs = return_nuc_content_from_codon_data(data[i + 1])
                nucs["host"] = "Unknown"
                nucs["tax_id"] = tax_id
                nucs["cds_num"] = num_of_seqs
                nucs["name"] = name
                res = res.append(nucs, ignore_index=True)

        res.to_csv(out_file)
        res_all = res_all.append(res, ignore_index=True)
    res_all.to_csv("/sternadi/home/volume3/taliakustin/phyVirus_analysis/hosts/host_count_nuc_coding_seq_info.csv")







if __name__ == "__main__":
    main()

