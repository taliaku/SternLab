import pandas as pd
import math
from tqdm import tqdm
import argparse
import re

def compare_fastas_to_ref(fastas, ref_seq_name, output_csv, remove_edges=True):
    # fasta needs to be aligned, ref_seq_name is one of the sequences in the aligned fasta.
    with open(fastas) as f:
        f = f.read()
    f = f.split('>')
    f = {i.split('\n')[0]:''.join(i.split('\n')[1:]) for i in f if i != ''}

    # remove read that somehow is shorter after msa...
    good_length = len(f[ref_seq_name])
    bad_seqs = [i for i in f if len(f[i]) != good_length]
    for i in bad_seqs:
        print('dropping seq ', i)
        f.pop(i)

    diffs = []
    ref_pos = 0
    for i in tqdm(range(len(f[ref_seq_name]))):
        # handle mismatches and deletions
        if f[ref_seq_name][i] != '-':
            ref_pos = math.floor(ref_pos) + 1
            for sample in f:
                if f[sample][i] != f[ref_seq_name][i]:
                    diffs.append((ref_pos, sample, f[ref_seq_name][i], f[sample][i], i))

        # handle insertions
        else:
            ref_pos += 0.001
            for sample in f:
                if f[sample][i] != f[ref_seq_name][i]:
                    diffs.append((ref_pos, sample, f[ref_seq_name][i], f[sample][i], i))

    df = pd.DataFrame(diffs, columns=['position', 'sample', 'ref_base', 'base', 'fasta_position'])
    print('diff count: {}'.format(len(diffs)))

    if remove_edges == 'y':
        # remove leading and trailing Ns and -s per sequence
        edges = []
        for sample in tqdm(f):
            start_pos = re.search(r'[^Nn-]', f[sample]).start()
            end_pos = good_length - re.search(r'[^Nn-]', f[sample][::-1]).start()
            edges.append((sample, start_pos, end_pos))

        edges = pd.DataFrame(edges, columns=['sample', 'start_pos', 'end_pos'])
        print(edges)
        df = pd.merge(df, edges, on='sample')
        df = df[(df.fasta_position >= df.start_pos) & (df.fasta_position <= df.end_pos)]

    df.drop(columns=['fasta_position'], inplace=True)
    df.sort_values(['sample', 'position']).to_csv(output_csv, index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", type=str, help="aligned fasta file",
                        required=True)
    parser.add_argument("-r", "--reference", type=str, help='reference sequence name, one of the sequences in the aligned fasta', required=True)
    parser.add_argument("-o", "--output", type=str, help='output csv', required=True)
    parser.add_argument("-e", '--remove_edges', type=str, help='remove the ends with deletions and Ns? y or n', default='n', required=False)
    args = parser.parse_args()
    if not vars(args):
        parser.print_help()
        parser.exit(1)
    compare_fastas_to_ref(args.input, args.reference, args.output, args.remove_edges)
