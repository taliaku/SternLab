import pandas as pd
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('files', help='files to concat', type=argparse.FileType('r'), nargs='+')
    parser.add_argument('-o', '--out_file', help='out file', type=str)

    args = parser.parse_args()
    df = pd.concat([pd.read_csv(i, low_memory=False) for i in args.files])
    df.to_csv(args.out_file, index=False)
