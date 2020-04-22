#! /usr/local/python_anaconda/bin/python3.4


import os
import sys
sys.path.append('/sternadi/home/volume1/shared/SternLab')
from blast_utilities import blast_to_mutations_list, blast_to_df
import argparse
import pandas as pd


def analyze_sanger(input_folder, output_folder, percent_identity, reference_path):
#    if not os.path.isdir(input_folder):
#        return 'Input folder does not exist'
#    if output_folder == False:
#        if input_folder.endswith('/'):
#            input_folder = input_folder[:-1]
#        output_folder = input_folder + '_out'
#    if not os.path.exists(output_folder):
#        os.mkdir(output_folder)
#    
#    files = [f for f in os.listdir(input_folder) if f.endswith('.fasta')]
#    sangers_string = ''
#    for f in files:
#        sangers_string += '>' + input_folder + '/' + f
#        with open(input_folder + '/' + f) as f:
#            s = f.read()
#            sangers_string += s[s.find('\n'):]  + '\n'
#    with open(output_folder + '/sangers.fasta', 'w') as f:
#        f.write(sangers_string)
    
    os.system('/sternadi/home/volume1/shared/tools/ncbi-blast-2.2.30+/bin/makeblastdb -in ' + output_folder + '/sangers.fasta' + ' -dbtype nucl')
    os.system('/sternadi/home/volume1/shared/tools/ncbi-blast-2.2.30+/bin/blastn -query ' + args.reference + ' -task blastn -db ' + output_folder + '/sangers.fasta' + ' -outfmt "6 sseqid qstart qend qstrand sstart send sstrand length btop" -num_alignments 1000000 -dust no -soft_masking F -perc_identity ' + str(percent_identity) + ' -evalue 1e-07 -out ' + output_folder + '/sangers.blast')
    
    mutations = blast_to_mutations_list(output_folder + '/sangers.blast')
    blast_df = blast_to_df(output_folder + '/sangers.blast')
    df = pd.merge(mutations, blast_df, on='read', how='right')#.fillna(0)
    df = df[['read', 'position', 'ref', 'base']]
    df = df.rename(columns={'read':'file'})
    df = df.drop_duplicates()
    df.to_csv(output_folder + '/sangers.blast.mutations_list.csv', index=False)
    return  'Finished'
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input_dir", type=str, help="input directory with .fasta files",
                        required=True)
    parser.add_argument("-o", "--output", type=str, help="a path to an output directory, optional", default=False, required=False)
    parser.add_argument("-d", "--percent_identity", type=str, help='percent identity for blast, default 85', default="85", required=False)
    parser.add_argument("-r", "--reference", type=str, help='path to reference', default="/sternadi/home/volume2/noam/ms2_reference/noam_concensus_ref.txt", required=False)
    args = parser.parse_args()
    if not vars(args):
        parser.print_help()
        parser.exit(1)

    print(analyze_sanger(args.input_dir, args.output, args.percent_identity, args.reference))
