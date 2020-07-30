import os
import sys
import argparse
import numpy as np
import pandas as pd
import multiprocessing as mp
from scipy.stats import fisher_exact


def main(args):
    freqs_file = args.freqs_file
    blast_output = args.blast_output
    mutations_all = args.mutations_all
    output_folder = args.output_folder
    os.makedirs(output_folder, exist_ok=True)
    input_x = str(args.position)
    if '-' in input_x:
        start_pos, end_pos = input_x.split('-')
        pool = mp.Pool(processes=100)
        results = {pos: pool.apply(get_variant, args=(pos, freqs_file, blast_output, mutations_all, output_folder))
                   for pos in range(int(start_pos), int(end_pos))}
    else:
        results = {input_x: get_variant(input_x=int(input_x), freqs_file=freqs_file, blast_output=blast_output,
                                        mutations_all=mutations_all, output_folder=output_folder)}
    for pos, output_strings in results.items():
        if len(output_strings) != 0:
            with open(os.path.join(output_folder, f"{pos}.txt"), 'w') as text_file:
                for line in output_strings.values():
                    print(line, file=text_file)


def get_variant(input_x, freqs_file, blast_output, mutations_all, output_folder):
    with open(os.path.join(output_folder, f"{input_x}.test"), 'w') as text_file:
        print(f"started variant {input_x}", file=text_file)
    freqs = pd.read_csv(freqs_file, sep="\t")
    freqs = freqs[freqs['Pos'] == np.round(freqs['Pos'])]  #remove insertions
    if (input_x < freqs["Pos"].min()) or (input_x > freqs["Pos"].max()):
        return {}

    all_mappings = pd.read_csv(blast_output, names=["read_id", "start", "end", 'read_start', 'read_end',
                                                         'plus_or_minus', 'length', 'mutations'], sep="\t")
    all_mutations = pd.read_csv(mutations_all, names=["pos", "read_id", "mutant", "read_positions"],
                                skiprows=[0], dtype=str, sep="\t")
    all_mutations = all_mutations[all_mutations.pos != 'ref_pos']
    all_mutations['pos'] = all_mutations['pos'].astype(int)
    cons = freqs[(freqs["Rank"] == 0)
                 & (freqs["Base"] != "-")]
    cons.insert(0, "pos", pd.to_numeric(cons.loc[:, "Pos"]))

    all_mutations = pd.merge(all_mutations, cons[["pos","Ref"]], on="pos")
    #Identify co-occurring variants up to max overlap length - 250 bases
    variants_combinations = range(input_x+1, input_x+250)
    output_dict = {}
    for y in variants_combinations:
        x = input_x
        maps_for_two_pos = all_mappings[(all_mappings["start"] <= x) & (all_mappings["end"] >= y)]
        grouped = maps_for_two_pos.groupby('read_id')
        grouped = grouped.filter(lambda x: len(x) == 2)
        merged = pd.merge(pd.DataFrame({"read_id":grouped["read_id"].unique()}), all_mutations[all_mutations["pos"]==x][["pos","read_id"]], on="read_id", how="left")
        merged = pd.merge(merged, all_mutations[all_mutations["pos"]==y][["pos","read_id"]], on="read_id", how="left")
        x_label = "pos_" + str(x)
        y_label = "pos_" + str(y)
        merged[x_label] = np.where(merged["pos_x"] == x, 1, 0)
        merged[y_label] = np.where(merged["pos_y"] == y, 1, 0)
        ct = pd.crosstab(merged[x_label], merged[y_label])
        if ct.shape == (2,2):
            fisher_test = fisher_exact(ct, alternative='greater')
            output_dict[y] = '\t'.join([str(x) for x in [x, y, fisher_test[0], fisher_test[1], ct[1][1]*1.0/(ct[0][0]+ct[0][1]+ct[1][0]+ct[1][1])]])
        else:
            output_dict[y] = '\t'.join([str(x) for x in [x, y, 0.0, 1.0, 0.0]])
    return output_dict


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-b", "--blast_output", type=str, help="all BLAST output for this sample")
    parser.add_argument("-m", "--mutations_all", type=str, help="mutations_all.txt file (filtered from text)")
    parser.add_argument("-p", "--position", help="The position to consider pairs, entering a range like '30-100' would "
                                                 "iterate over all positions in range")
    parser.add_argument("-f", "--freqs_file", type=str, help="freqs file")
    parser.add_argument("-o", "--output_folder", type=str, help="where to output files")
    args = parser.parse_args(sys.argv[1:])
    main(args)
