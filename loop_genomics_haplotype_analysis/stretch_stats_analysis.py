import pandas as pd

# TODO- re-verify entire method
def stretches_stats_per_sample(stretch_file, freq_file):
    pairs_by_stretch = pd.read_csv(stretch_file, sep='\t')
    # print(pairs_by_stretch.head())
    print(pairs_by_stretch.shape)

    # filter pair duplication- keep only 'forward' side
    pairs_by_stretch = pairs_by_stretch[pairs_by_stretch['Pos1'] < pairs_by_stretch['Pos2']]

    # flatten pairs- each mutation to a separated line
    pos1 = pairs_by_stretch[['Stretch', 'Pos1', 'Freq', 'meandist']]
    pos2 = pairs_by_stretch[['Stretch', 'Pos2', 'Freq', 'meandist']]
    mutations_by_stretch = pd.concat([pos1.rename(columns={'Pos1':'Pos'}),
                                     pos2.rename(columns={'Pos2':'Pos'})])

    # TODO- understand how to handle position duplicates
    # check_duplicates = mutations_by_stretch[mutations_by_stretch.duplicated(['Stretch', 'Pos'])]
    mutations_by_stretch = mutations_by_stretch.sort_values(by=['Pos']).drop_duplicates(['Stretch', 'Pos'])
    # print(mutations_by_stretch.head())
    print(mutations_by_stretch.shape)

    # Add mutation identity- Rank0 -> Rank1 (heuristic)
    freqs = pd.read_csv(freq_file, sep='\t')
    ref = freqs[freqs['Rank'] == 0][['Pos','Base']]
    base = freqs[freqs['Rank'] == 1][['Pos','Base']]

    mutations_by_stretch = pd.merge(mutations_by_stretch, ref, how = 'left', on= 'Pos').rename(columns={'Base':'Rank0'})
    mutations_by_stretch = pd.merge(mutations_by_stretch, base, how = 'left', on= 'Pos').rename(columns={'Base':'Rank1'})
    mutations_by_stretch['mutation'] = mutations_by_stretch['Rank0'] + mutations_by_stretch['Rank1']
    print(mutations_by_stretch.shape)
    # print(mutations_by_stretch.head())

    # mutation count
    count = mutations_by_stretch.groupby('Stretch').count()['Pos'].sort_values(ascending=False)
    print(count.head())

    # strech length
    min = mutations_by_stretch.groupby('Stretch').min()['Pos']
    max = mutations_by_stretch.groupby('Stretch').max()['Pos']

    # filter stretches
    top_stretches_count = count.head().index
    print(top_stretches_count.tolist())
    mutations_by_stretch = mutations_by_stretch[mutations_by_stretch['Stretch'].isin(top_stretches_count)]

    # TODO- filter by meandist? currently not
    freq = pairs_by_stretch.groupby('Stretch').first()['meandist'].sort_values(ascending=False)
    # print(freq.head())
    # top_stretches_freq = freq.head().index
    # print(top_stretches_freq.tolist())
    # mutations_by_stretch = mutations_by_stretch[mutations_by_stretch['Stretch'].isin(top_stretches_freq)]

    # mutation type summary
    mutation_types_summary = mutations_by_stretch.groupby(['Stretch', 'mutation']).count()['Pos'].sort_values(ascending=False).sort_index(level='Stretch', sort_remaining=False)
    mutation_types_summary = pd.merge(mutation_types_summary.reset_index(), count.reset_index(), on='Stretch', how='left')
    mutation_types_summary = mutation_types_summary.rename(columns={'Pos_x':'count', 'Pos_y':'total'})
    mutation_types_summary['ratio'] = mutation_types_summary['count'] / mutation_types_summary['total']

    # min, max, length columns
    mutation_types_summary = pd.merge(mutation_types_summary, min.reset_index(), on='Stretch',how='left')
    mutation_types_summary = mutation_types_summary.rename(columns={'Pos':'first_pos'})
    mutation_types_summary = pd.merge(mutation_types_summary, max.reset_index(), on='Stretch',how='left')
    mutation_types_summary = mutation_types_summary.rename(columns={'Pos':'last_pos'})
    mutation_types_summary['length'] = mutation_types_summary['last_pos'] - mutation_types_summary['first_pos']
    mutation_types_summary = pd.merge(mutation_types_summary, freq.reset_index(), on='Stretch', how='left')

    return mutation_types_summary


def combine_stretches_stats_all_sampels():
    global mutation_types_summary
    all_summaries = []
    for sample in env_samples:
        print('Running: {}'.format(sample))
        stretch_file = '{}/{}s_output/{}/analysis/stretches.csv'.format(root_dir, sample[:3], sample)
        freq_file = '{}/{}s_output/{}/joined.freqs'.format(root_dir, sample[:3], sample)
        mutation_types_summary = stretches_stats_per_sample(stretch_file, freq_file)

        # extract sample summary
        # mutation_types_summary.to_csv('{}/{}_mutation_summary.csv'.format(output_dir, sample), index=False)

        # global summary
        mutation_types_summary['Stretch'] = sample + '_' + mutation_types_summary['Stretch'].astype(str)
        mutation_types_summary_p = mutation_types_summary.pivot(index='Stretch',
                                                                columns='mutation',
                                                                values='ratio')
        mutation_types_summary_p = pd.merge(mutation_types_summary_p.reset_index(), mutation_types_summary[
            ['Stretch', 'total', 'length', 'first_pos', 'last_pos', 'meandist']].drop_duplicates(), on='Stretch',
                                            how='left')

        all_summaries.append(mutation_types_summary_p)
    # merge
    global_summary = pd.concat(all_summaries, sort=False)
    # export
    # TODO- rename columns & set order
    # print(global_summary)
    global_summary.to_csv('{}/global_summary_env_v3.csv'.format(output_dir), index=False)
    # mutation types in rows instead of columns-
    global_summary = global_summary.transpose()
    new_header = global_summary.iloc[0]
    global_summary = global_summary[1:]
    global_summary.columns = new_header
    global_summary.to_csv('{}/global_summary_env_v3_for_pres.csv'.format(output_dir))


def hiv_shafer_analysis():
    global mutation_types_summary, root_dir, output_dir, env_samples
    test = False
    if test:
        mutation_types_summary = stretches_stats_per_sample(
            # '/sternadi/home/volume2/ita/haplotypes/envs_output/env1/analysis/stretches.csv',
            # '/sternadi/home/volume2/ita/haplotypes/envs_output/env1/joined.freqs'
            '/Volumes/STERNADILABHOME$/volume2/ita/haplotypes/envs_output/env11/analysis/stretches.csv',
            '/Volumes/STERNADILABHOME$/volume2/ita/haplotypes/envs_output/env11/joined.freqs'
        )
        print(mutation_types_summary)

        # extract
        # mutation_types_summary.to_csv('/Volumes/STERNADILABHOME$/volume2/ita/haplotypes/envs_output/env11/analysis/mutation_type_summary.csv')
        mutation_types_summary.to_csv(
            '/Users/omer/Documents/Lab/HIV_shafer/mutation_type_summaries/mutation_type_summary.csv', index=False)

    else:
        root_dir = '/Volumes/STERNADILABHOME$/volume2/ita/haplotypes'
        output_dir = '/Volumes/STERNADILABHOME$/volume1/shared/analysis/HIV_shafer/accungs_haplotype_inference/mutation_type_summaries'

        samples = ['env1', 'env10', 'env11', 'env12', 'env13', 'env14', 'env15', 'env2', 'env3', 'env4', 'env5', 'env6',
                   'env7', 'env8', 'env9',
                   'gag1', 'gag10', 'gag11', 'gag12', 'gag13', 'gag14', 'gag15', 'gag2', 'gag3', 'gag4', 'gag5', 'gag6',
                   'gag7', 'gag8', 'gag9']

        env_samples = ['env1', 'env10', 'env11', 'env12', 'env13', 'env14', 'env15', 'env2', 'env3', 'env4', 'env5',
                       'env6', 'env7', 'env8', 'env9']
        gag_samples = ['gag10', 'gag11', 'gag12', 'gag13', 'gag14', 'gag15', 'gag2', 'gag3', 'gag5', 'gag6', 'gag8',
                       'gag9']

        # TODO- use full samples (stretch file of env+gag combined)

        combine_stretches_stats_all_sampels()


if __name__ == '__main__':

    hiv_shafer_analysis()

