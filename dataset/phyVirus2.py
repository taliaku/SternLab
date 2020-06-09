#! /usr/local/python_anaconda/bin/python3.4

from optparse import OptionParser
from file_utilities import check_dirname, check_filename, change_filename
import glob
from seqFileAnalyzer import get_consensus_percentage
import pandas as pd
from PAML_utilities import write_ctl_file


def make_dirs_for_each_groups_and_split_fasta():
    files = glob.glob("/sternadi/home/volume3/taliakustin/phyVirus2/*fasta")
    for f in files:
        base = f.split("/")[-1].split(".fasta")[0]
        os.system("mkdir /sternadi/home/volume3/taliakustin/phyVirus2/split_seqs/%s/" % base)
        os.system("cp %s %s" % (f, "/sternadi/home/volume3/taliakustin/phyVirus2/split_seqs/%s/%s.fasta" % (base, base)))
        seqFileTools.split_fasta_file_per_seq("/sternadi/home/volume3/taliakustin/phyVirus2/split_seqs/%s/%s.fasta" % (base, base))


def run_blast_on_all_files():
    files = glob.glob("/sternadi/home/volume3/taliakustin/phyVirus2/split_seqs/*/*fasta")
    for f in files:
        base = f.split("/")[-1].split("_")[0]
        db = "/sternadi/home/volume3/taliakustin/phyVirus2/nuc_db/%s.fasta" % base
        out = f.split("fasta")[0] + "dc-megablast_out.txt"
        pbs_runners.script_runner(
            "/sternadi/home/volume1/shared/tools/ncbi-blast-2.2.30+/bin/blastn -query %s -out %s -db %s -outfmt 6 -task dc-megablast" % (
            f, out, db), alias="dc-megablast-all")



groups = {}

group_counter = 1

for f in tqdm.tqdm(files):
    with open(f, "r") as file:
    data = file.readlines()
    res = []
    for l in data:
        if l == data[0]:
            res.append(l.split("\t")[0])
        res.append(l.split("\t")[1])
    exist = False
    selected = 0
    for r in res:
        for g in groups:
            if r in groups[g]:
            exist = True
            selected = g
    if exist:
        for r in res:
            if r in groups[selected]:
                continue
            else:
                groups[selected].append(r)
    else:
        groups[group_counter] = []
        for r in res:
            groups[group_counter].append(r)
        group_counter += 1
        print("group+1 = %i" % group_counter)




    seqs = list(SeqIO.parse("/sternadi/home/volume3/taliakustin/phyVirus2/Arena.fasta", format="fasta"))
    groups_proteins = []

    for s in seqs:
             for g in groups:
                     if s.description in groups[g]:
                             if g not in groups_proteins.keys():
                                 groups_proteins[g] = []
                     groups_proteins[g].append(s)

    for g in groups_proteins:
             SeqIO.write(groups_proteins[g], outfiles + str(g) + ".fasta", format="fasta")


        
