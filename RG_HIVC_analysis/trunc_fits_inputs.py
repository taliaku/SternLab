import argparse

import pandas as pd


def trunc_file(args):
    gen_threshold = 350
    if args.gen != None:
        gen_threshold = int(args.gen)

    df = pd.read_csv(args.in_file, sep='\t')
    # filter Gen > 400
    df = df[df['Gen'] < gen_threshold]
    out_file = args.in_file + '.trunc'
    df.to_csv(out_file, sep='\t', encoding='utf-8', index=False)



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--in_file", type=str, required=True)
    parser.add_argument("-g", "--gen", type=str, required=False)
    args = parser.parse_args()
    trunc_file(args)