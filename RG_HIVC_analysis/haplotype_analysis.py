import pandas as pd
import numpy as np
import math
import matplotlib.pyplot as plt
import itertools
import sys
import argparse
from scipy.stats import chi2_contingency, fisher_exact
import seaborn as sns;

from RG_HIVC_analysis.constants import orig_first_samples, orig_high_quality_samples, RT_short_ET86_interval, \
    control_excluded_patients, control_bl_samples
from freqs_utilities import change_ref_to_consensus

sns.set_context("poster")

# TODO- finish reviewing all different pipeline output types.
# TODO- learn blast files & algorithm (2nd iteration- understand pipeline's full flow)

# TODO- go over inclarities in flow (with maoz)

def mutations_association(args):
    # position to handle (will process interval of 250 positions ahead)
    # TODO- for AA muts version- change code to run over all x in [1939,2629] interval, each time comparing 2\1\0 places ahead according to position in codon (no need to run only 3rd position)
    input_x=args.position

    freqs=pd.read_csv(args.freqs_file, sep="\t")
    freqs = freqs[freqs['Pos'] == np.round(freqs['Pos'])] #remove insertions #TODO- maoz- why not simply ref != '-'
    if (input_x < freqs["Pos"].min()) or (input_x > freqs["Pos"].max()):
        sys.exit()

    # blast files (all .fasta.blast files joined together)
    all_mappings=pd.read_csv(args.blast_output, names=["read_id","start","end"], sep="\t")
    # summary of all observed mutations from ref, including mappings to origin reads
    all_mutations=pd.read_csv(args.mutations_all, names=["pos","read_id","mutant","read_positions"], sep="\t")  #all delta's between single reads to the reference. displays the mutated base (origin is known by ref), with its position

    cons = freqs[(freqs["Rank"] == 0)
                 & (freqs["Base"] != "-")]
    cons.insert(0, "pos", pd.to_numeric(cons.loc[:,"Pos"]))

    all_mutations = pd.merge(all_mutations, cons[["pos","Ref"]], on="pos") # adding Ref\Cons to all_mutations
    #remove C>A and G>T #TODO- common seq errors- can consider filtering
    #all_mutations = all_mutations[~(((all_mutations["Ref"]=="C")&(all_mutations["mutant"]=="A")) | ((all_mutations["Ref"]=="G")&(all_mutations["mutant"]=="T")))]

    #variants=all_mutations["pos"].unique()

    variants_combinations=range(input_x+1,input_x+2) # x-> (x+1,x+2) instead of (x+1,x+250)

    for y in variants_combinations:
        x=input_x
        maps_for_two_pos = all_mappings[(all_mappings["start"] <= x) & (all_mappings["end"] >= y)] # reads surrounding the [x,y] interval
        merge_read_id = pd.DataFrame({"read_id": maps_for_two_pos["read_id"].unique()})
        merge_x = all_mutations[all_mutations["pos"] == x][["pos", "read_id"]]
        merged = pd.merge(merge_read_id, merge_x, on="read_id", how="left")
        merge_y = all_mutations[all_mutations["pos"] == y][["pos", "read_id"]]
        merged = pd.merge(merged, merge_y, on="read_id", how="left")

        x_label = "pos_" + str(x)
        y_label = "pos_" + str(y)
        merged[x_label] = np.where(merged["pos_x"] == x, 1, 0)
        merged[y_label] = np.where(merged["pos_y"] == y, 1, 0)
        ct = pd.crosstab(merged[x_label], merged[y_label])
        if ct.shape == (2,2):
            fisher_test = fisher_exact(ct, alternative='greater') # greater- aims to the bottom-right cell? # TODO- review fisher's test
            print('\t'.join([str(x) for x in [x, y, fisher_test[0], fisher_test[1], ct[1][1]*1.0/(ct[0][0]+ct[0][1]+ct[1][0]+ct[1][1])]]))
        else:
            print('\t'.join([str(x) for x in [x, y, 0.0, 1.0, 0.0]])) # statistic: odds ratio- ratio of diagonals, p-value, *shared_freq*
            # list of all connections of pair mutations-
            # 1. should be grouped into cliques with one joint freq?- and mapped to freq plot with color per clique?
            # 2. adjacent positions should be registered with joint freq, and given as input to aa_mut script?


def freq_plot(freq_df, fname, plot_header=''):

    print('plotting freq_plot')
    g = sns.relplot(x= "Pos",
                    y= "Freq",
                    col='sample_id',
                    # col_order='years_since_infection', # chronological presentation
                    # # hue='sample_id',
                    hue='mutation_type',
                    col_wrap=7,
                    # join=True,
                    data=freq_df)


    # plot adjustments
    g.set(yscale="log")
    plt.ylim(5e-4, 1)
    g.fig.suptitle(plot_header, y= 0.1)
    # g.set_ylabels("mutation_rate")
    # g.set_xlabels("ET first sample (Years)")
    # g.set_xticklabels(rotation=45, fontsize=11)

    # extract plot
    # plt.show()
    plt.savefig(fname=fname)
    # g.savefig('')

def dist_plot(freq_df):
    # g = sns.jointplot(x= "Pos", y= "Freq", data=freq_df, kind='kde')

    # plot 1
    # f, ax = plt.subplots(figsize=(10, 6))
    # sns.kdeplot(freq_df.Pos, freq_df.Freq, ax=ax)
    # sns.rugplot(freq_df.Freq, vertical=True, ax=ax)

    # plot 2
    # samples = ['X83354_S92','504201_S45','X100748_S73','504223_S67']
    samples = ['87415_S75','504194_S38','504202_S46','TASPX119494_S74']
    # f, axes = plt.subplots(2, math.ceil(len(samples)/2), figsize=(6, 6))
    # for i in range(len(samples)):
    #     sns.jointplot(x='Pos', y='Freq',
    #                  data=freq_df[freq_df['sample_id'].isin([samples[i]])],
    #                  kind= 'kde',
    #                  legend=False,
    #                  ax=axes[round(i/len(samples))][i % math.ceil(len(samples)/2)])


    # plot 3
    plt.xscale('log')
    f, axes = plt.subplots(1, len(samples), figsize=(len(samples)*6, 6))
    for i in range(len(samples)):
        sns.distplot(freq_df[freq_df['sample_id'].isin([samples[i]])]['Freq'],
                     ax=axes[i])
        # sns.distplot(freq_df.Freq);

    # plt.ylim(10**-4, 1)
    # ax.set_xticks(range(0, 9000, 1000))
    plt.show()

def dist_plot_v2(freq_df, fig_path, run_name, file_header=''):
    print('plotting dist_plot_v2')
    for sample in freq_df.sample_id.unique():
    # for sample in ['504234_S78','504230_S74','504233_S77']:
        freqs_per_sample = freq_df[(freq_df.sample_id == sample)]['Freq']
        print(len(freqs_per_sample))
        ax = sns.distplot(freqs_per_sample)
        plot_header = file_header + '_' + str(sample)
        ax.set_title(plot_header)
        # ax.set(xscale="log")

        # plt.show()
        dist_fig_path = fig_path + '/mut_freq_distribution_figs/' + run_name
        plt.savefig(dist_fig_path + '/' + plot_header + '_std_scale.png')
        plt.cla()

def main_plots():
    # get freq files
    fig_path = '/Users/omer/PycharmProjects/SternLab/RG_HIVC_analysis/freq_trajectories/'
    run_name = 'orig_high'
    freq_df = pd.read_csv('/Users/omer/PycharmProjects/SternLab/RG_HIVC_analysis/runs/%s/unified_freqs_filtered_verbose.csv' % run_name)

    # force ref=consensus in each sample (rank 0 == ref)
    # TODO- redundant for new unified_freqs_filtered_verbose.csv
    dfs = []
    for sample in freq_df.sample_id.unique():
        df_sample = freq_df[freq_df['sample_id'] == sample]
        # TODO- remove indels before change_ref_to_consensus()
        # df_sample_con_as_ref = change_ref_to_consensus(df_sample)
        dfs.append(df_sample_con_as_ref)

    freq_df = pd.concat(dfs)


    # filters
    freq_df = freq_df[freq_df['Rank'] == 1]
    freq_df = freq_df[freq_df['Read_count'] > 100]
    freq_df = freq_df[freq_df['Freq'] != 0]
    # freq_df = freq_df[freq_df['Freq'] > 0.01]
    freq_df = freq_df[freq_df['Prob'] >= 0.95]
    # freq_df = freq_df[(freq_df["Pos"] >= RT_short_ET86_interval[0]) & (freq_df["Pos"] <= RT_short_ET86_interval[1])]

    # transitions only
    # freq_df["mutation_type"] = freq_df['Ref'] + freq_df['Base']
    # freq_df = freq_df[freq_df["mutation_type"].isin(['GA', 'AG', 'GG', 'AA', 'CT', 'TC', 'CC', 'TT'])]
    freq_df["mutation_type"] = freq_df['Ref'] + '->' + freq_df['Base']
    mismathces = freq_df[freq_df["mutation_type"].isin(['G->G', 'A->A', 'C->C', 'T->T'])]
    freq_df = freq_df[freq_df["mutation_type"].isin(['G->A', 'A->G', 'G->G', 'A->A', 'C->T', 'T->C', 'C->C', 'T->T'])]

    # choose patients\samples
    # patients = ['29447']
    # patients = control_excluded_patients
    # freq_df = freq_df[freq_df['ind_id'].isin(patients)]
    # freq_df = freq_df[~freq_df['ind_id'].isin(patients)]

    samples = orig_first_samples
    # samples = control_bl_samples

    # samples = orig_high_quality_samples
    # samples = ['87415_S75', '504194_S38', '504202_S46', 'TASPX119494_S74'] # patient 15664
    # samples = set(orig_first_samples).intersection(set(orig_high_quality_samples))

    freq_df = freq_df[freq_df['sample_id'].isin(samples)]

    file_header= '%s_rank1_cov_100_prob_95_transitions_founders_fixed_consensus' % run_name
    # plot_header= '%s/Rank1, transitions only, founder samples, cov>100 per pos, prob >=0.95' % run_name
    plot_header= file_header

    fname = fig_path + '/muts_per_position_' + file_header
    # report underlying statistics
    mut_per_sample_stats = freq_df.groupby('sample_id').count()['Rank']
    print('#mut per sample stats: \n' + str(mut_per_sample_stats.describe()))
    mut_per_type_stats = freq_df.groupby('mutation_type').count()['Rank']
    print('#mut per type stats: \n' + str(mut_per_type_stats))

    # exporting to file
    f = open(fname + '.txt', 'w')
    f.write('#mut per sample stats: \n' + str(mut_per_sample_stats.describe()))
    f.write('#mut per type stats: \n' + str(mut_per_type_stats))
    f.close()

    # plot1
    plot_fname = fname + '.png'
    freq_plot(freq_df, plot_fname, plot_header)

    # plot2
    # dist_plot_v2(freq_df, fig_path, run_name, file_header)


def main_mutations_association():
    global args
    parser = argparse.ArgumentParser()
    parser.add_argument("blast_output", type=str, help="all BLAST output for this sample")
    parser.add_argument("mutations_all", type=str, help="mutations_all.txt file (filtered from text)")
    parser.add_argument("position", type=int, help="The position to consider pairs")
    parser.add_argument("freqs_file", type=str, help="freqs file")
    args = parser.parse_args(sys.argv[1:])

    # inferring co-occurence for a single position, with all the rest
    mutations_association(args)


if __name__ == "__main__":
    main_plots()
    # main_mutations_association()