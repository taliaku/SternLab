#! /usr/local/python_anaconda/bin/python3.4

import os
import sys
sys.path.append('/sternadi/home/volume1/shared/SternLab')
from blast_utilities import blast_to_mutations_list, blast_to_df
import argparse
import pandas as pd

def analyze_sanger(input_fastas_file, percent_identity, reference_path):    
    os.system('/sternadi/home/volume1/shared/tools/ncbi-blast-2.2.30+/bin/makeblastdb -in ' + input_fastas_file + ' -dbtype nucl')
    os.system('/sternadi/home/volume1/shared/tools/ncbi-blast-2.2.30+/bin/blastn -query ' + args.reference + ' -task blastn -db ' + input_fastas_file + ' -outfmt "6 sseqid qstart qend qstrand sstart send sstrand length btop" -num_alignments 1000000 -dust no -soft_masking F -perc_identity ' + str(percent_identity) + ' -evalue 1e-07 -out ' + input_fastas_file + '.blast')
    
    mutations = blast_to_mutations_list(input_fastas_file + '.blast')
    blast_df = blast_to_df(input_fastas_file + '.blast')
    df = pd.merge(mutations, blast_df, on='read', how='right')#.fillna(0)
    df = df[['read', 'position', 'ref', 'base']]
    df = df.rename(columns={'read':'file'})
    df = df.drop_duplicates()
    df.to_csv(input_fastas_file + '.mutations_list.csv', index=False)
    return  'Finished'

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", type=str, help="fasta file with sequences to compare to referece",
                        required=True)
    parser.add_argument("-d", "--percent_identity", type=str, help='percent identity for blast, default 85', default="85", required=False)
    parser.add_argument("-r", "--reference", type=str, help='path to reference', required=False)
    args = parser.parse_args()
    if not vars(args):
        parser.print_help()
        parser.exit(1)
    print(analyze_sanger(args.input, args.percent_identity, args.reference))