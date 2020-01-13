#!/usr/bin/env python

import sys
#sys.path.append('/sternadi/home/volume1/shared/SternLab')
from blast_utilities import blast_to_mutations_list, blast_to_df
import pandas as pd
import os

def main(blasts_input_folder, output_folder):
    analysis = {'reads_from_blast': 0, 'reads_whole_amplicon': 0}
    primer_ids = pd.DataFrame()
    mutations = pd.DataFrame()
    blast_files = [os.path.join(blasts_input_folder, f) for f in os.listdir(blasts_input_folder) if f.endswith('.blast')]
    for f in blast_files:
        print(f)
        ##mutations_list = blast_to_mutations_list(f, f + '.mutations_list')
        blast_df = blast_to_df(f)
        analysis['reads_from_blast'] = analysis['reads_from_blast'] + len(blast_df)
        # get only reads where whole amplicon was sequenced
        blast_df = blast_df[((blast_df.start_ref == 1) & (blast_df.end_ref == 262))]
        ##blast_df = blast_df[((blast_df.end_ref == 262) & (blast_df.start_ref < 50))]

        blast_df = blast_df[['read']]
        analysis['reads_whole_amplicon'] = analysis['reads_whole_amplicon'] + len(blast_df)
        
        # get reads with matching primer ids
        mutations_list = pd.read_csv(f + '.mutations_list')
        n_mutations = mutations_list[(mutations_list.ref == 'N')].copy()
        # filter only reads that sequenced whole amplicon
        n_mutations = pd.merge(n_mutations, blast_df, on='read', how='inner')
        partial_primer_ids = get_primer_ids(n_mutations)
        # keep only primer ids without deletions or insertions
        partial_primer_ids = partial_primer_ids[~(partial_primer_ids['primer_id'].str.contains('-')) & (partial_primer_ids.primer_id.str.len() == 15)]        
        primer_ids = pd.concat([primer_ids, partial_primer_ids])

        # Check only 1764 mutations (position 164 in this reference)
        partial_mutations = mutations_list[(mutations_list.position == 64)].copy()
        mutations = pd.concat([mutations, partial_mutations])
        
    analysis['reads_with_primer_id'] = len(primer_ids)
    analysis['primer_ids'] = len(primer_ids[['primer_id']].drop_duplicates())
    # left inner merge so reads with no 1764 mutation are kept as well
    df = pd.merge(primer_ids, mutations, on='read', how='left').fillna(0)
    df['primer_id_counted'] = df.groupby('primer_id')['read'].transform('count')
    analysis['primer_ids_counted_more_than_once'] = len(df[df.primer_id_counted > 1]['primer_id'].drop_duplicates())
    analysis['reads_per_unique_primer_ids'] = get_reads_per_unique_primer_id(df)
    df = get_position_concensus(df)
    create_graph2(df, output_folder)
    analysis['concensus_stats'] = concensus_stats(df)
    # concensus when 90 percent of reads per primer_id have same base
    df = df[(df.base_percent_per_primer_id >= 0.9)]
    ##df = df[(df.base_percent_per_primer_id > 0.5)]
    df = df[['primer_id', 'primer_id_counted', 'base', 'base_per_primer_id_counted', 'base_percent_per_primer_id']].drop_duplicates()
    analysis['primer_ids_with_concensus'] = len(df)
    analysis['primer_ids__count_more_than_one_concensus'] = len(df[(df.primer_id_counted > 1)])
    analysis['base_frequencies'] = df.groupby('base')['primer_id'].count() / len(df)
    analysis['primer_id_count_groupby_concensus'] = df.groupby('base')['primer_id_counted'].describe()
    create_graphs(analysis, output_folder)
    open(output_folder + '1664primer_id_results.txt', 'w').write(pretty_print_analysis(analysis))
    return df
   


def get_primer_ids(df):
    df = df.sort_values(['read', 'position'])
    df['primer_id'] = df.groupby('read')['base'].transform(''.join)
    df = df[['read', 'primer_id']].drop_duplicates()
    return df

def get_position_concensus(df):
    df['base_per_primer_id_counted'] = df.groupby(['primer_id', 'base'])['read'].transform('count')
    df = df[['primer_id', 'base', 'primer_id_counted', 'base_per_primer_id_counted']].drop_duplicates()
    df['base_percent_per_primer_id'] = df['base_per_primer_id_counted'] / df['primer_id_counted']
    df['base_percent_per_pid_rounded'] = df.base_percent_per_primer_id.round(2)
    return df

def concensus_stats(df):
    df['max_percent_base'] = df.groupby('primer_id')['base_percent_per_primer_id'].transform('max')
    df = df[['primer_id', 'max_percent_base']].drop_duplicates()
    return df.max_percent_base.describe()

def get_reads_per_unique_primer_id(df):
    df = df[['primer_id', 'primer_id_counted']].drop_duplicates()
    return df.groupby('primer_id_counted').count()

def pretty_print_analysis(analysis):
    s = 'Reads from blast: %s\nReads with whole amplicon: %s\nReads with primer id: %s\nUnique primer ids: %s\nReads per unique primer id: \n%s\n\nPrimer ids with 1764 concensus: %s, %s percent,\nConcensus percent out of primer ids seen more than once: %s percent\n\nBase frequencies: %s\n\n Primer id count stats, grouped by concensus:\n%s\n\nConcensus_stats:\n%s]'
    data = (analysis['reads_from_blast'], analysis['reads_whole_amplicon'], analysis['reads_with_primer_id'], analysis['primer_ids'], str(analysis['reads_per_unique_primer_ids']), analysis['primer_ids_with_concensus'],  analysis['primer_ids_with_concensus'] / analysis['primer_ids'] * 100, analysis['primer_ids__count_more_than_one_concensus'] / analysis['primer_ids_counted_more_than_once'] * 100, str(analysis['base_frequencies']), str(analysis['primer_id_count_groupby_concensus']), analysis['concensus_stats'])
    return (s % data)

def create_graphs(analysis, output_folder):
    f = analysis['reads_per_unique_primer_ids'].plot(kind='bar', logy=True)
    fig = f.get_figure()
    fig.savefig(output_folder + '1664reads_per_unique_primer_id.png')
    f.clear()
    return
    
def create_graph2(df, output_folder):
    df['max_percent_rounded'] = df.groupby('primer_id')['base_percent_per_pid_rounded'].transform('max')
    df = df[['primer_id', 'primer_id_counted', 'max_percent_rounded']].drop_duplicates()
    f2 = df.groupby('max_percent_rounded')['primer_id'].count().plot(kind='bar')
    fig2 = f2.get_figure()
    fig2.savefig(output_folder + '1664primer_ids.png')
    f2.clear()
    return

'''
main('X:/volume2/noam/primer_id_250418/data/T0-37A_pipeline_output/tmp/', 'X:/volume2/noam/primer_id_250418/data/T0-37A_primer_id_output/')
main('X:/volume2/noam/primer_id_250418/data/T15-37A_pipeline_output/tmp/', 'X:/volume2/noam/primer_id_250418/data/T15-37A_primer_id_output/')
main('X:/volume2/noam/primer_id_250418/data/T30-37A_pipeline_output/tmp/', 'X:/volume2/noam/primer_id_250418/data/T30-37A_primer_id_output/')

main('X:/volume2/noam/primer_id_250418/data/T0-37B_pipeline_output/tmp/', 'X:/volume2/noam/primer_id_250418/data/T0-37B_primer_id_output/')
main('X:/volume2/noam/primer_id_250418/data/T15-37B_pipeline_output/tmp/', 'X:/volume2/noam/primer_id_250418/data/T15-37B_primer_id_output/')
main('X:/volume2/noam/primer_id_250418/data/T30-37B_pipeline_output/tmp/', 'X:/volume2/noam/primer_id_250418/data/T30-37B_primer_id_output/')

main('X:/volume2/noam/primer_id_250418/data/T0-WT_pipeline_output/tmp/', 'X:/volume2/noam/primer_id_250418/data/T0-WT_primer_id_output/')
main('X:/volume2/noam/primer_id_250418/data/T15-WT_pipeline_output/tmp/', 'X:/volume2/noam/primer_id_250418/data/T15-WT_primer_id_output/')
main('X:/volume2/noam/primer_id_250418/data/T30-WT_pipeline_output/tmp/', 'X:/volume2/noam/primer_id_250418/data/T30-WT_primer_id_output/')
main('X:/volume2/noam/primer_id_250418/data/T30-WT2_pipeline_output/tmp/', 'X:/volume2/noam/primer_id_250418/data/T30-WT2_primer_id_output/')
'''