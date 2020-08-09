import pandas as pd
import math
from tqdm import tqdm
import argparse

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
        if f[ref_seq_name][i] != '-':
            ref_pos = math.floor(ref_pos) + 1
            for sample in f:
                if f[sample][i] != f[ref_seq_name][i]:
                    diffs.append((ref_pos, sample, f[ref_seq_name][i], f[sample][i]))

        # handle insertions
        else:
            ref_pos += 0.001
            for sample in f:
                if f[sample][i] != f[ref_seq_name][i]:
                    diffs.append((ref_pos, sample, f[ref_seq_name][i], f[sample][i]))

    df = pd.DataFrame(diffs, columns=['position', 'sample', 'ref_base', 'base'])
    print('diff count: {}'.format(len(diffs)))

    if remove_edges == 'y':
        edges = []
        for sample in tqdm(df['sample'].drop_duplicates()):

            # TODO- drops first\last diff. fix
            a = df[df['sample'] == sample]
            a['position_diff'] = -a.position.diff(periods=-1)
            start_pos = a[(a.position_diff > 1)].position.min()
            a['position_diff'] = a.position.diff()
            end_pos = a[a.position_diff > 1].position.max()

            start_pos = start_pos if not math.isnan(start_pos) else df.position.min()
            end_pos = end_pos if not math.isnan(end_pos) else df.position.max()
            # TODO- should this happen? (start_pos\end_pos == NaN)

            edges.append((sample, start_pos, end_pos))

        edges = pd.DataFrame(edges, columns=['sample', 'start_pos', 'end_pos'])
        print(edges)
        df = pd.merge(df, edges, on='sample')
        df = df[(df.position > df.start_pos) & (df.position < df.end_pos)]

    return df.to_csv(output_csv, index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", type=str, help="aligned fasta file",
                        required=True)
    parser.add_argument("-r", "--reference", type=str, help='reference sequence name, one of the sequences in the aligned fasta', required=False)
    parser.add_argument("-o", "--output", type=str, help='output csv', required=False)
    parser.add_argument("-e", '--remove_edges', type=str, help='remove the ends with deletions and Ns? y or n', default='n', required=False)
    args = parser.parse_args()
    if not vars(args):
        parser.print_help()
        parser.exit(1)
    print(compare_fastas_to_ref(args.input, args.reference, args.output, args.remove_edges))
