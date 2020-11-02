#! /usr/local/python_anaconda/bin/python3.4


import pandas as pd
import argparse
import numpy as np


# TODO- all fixes -> code review
def count_haplotypes(input_chosen_mutations,
                     input_blast_df,
                     input_mutation_df,
                     minimal_mutation_frequency,
                     substitution_error_cutoff,
                     deletion_error_cutoff,
                     output_file):
    '''
    This function separates the reads into strains. A strain is defined as any combination 
    of the mutations or lack of mutations that are provided to the function.
    
    Format for the csv with the mutations to check strains for: every mutation should 
    have its own row, the csv should have a header row titled "variant", and the mutations 
    should be written in the following format: "A1664.0G". For example:
    "variant
    A1664.0G
    A535.0G
    T1764.0-"
    The output file variants_chosen.csv from association_test_variant.py can be used here.
    '''

    recognized_mutations = pd.read_csv(input_chosen_mutations)['variant'].tolist()
    recognized_positions = [float(p[1:-1]) for p in recognized_mutations]
    blast_df = pd.read_csv(input_blast_df)
    # choose only reads that were mapped only once in blast
    blast_df['read_count'] = blast_df.groupby('read')['start_ref'].transform('count')
    blast_df = blast_df[(blast_df.read_count == 1)]
    print('len(blast_df): {}'.format(len(blast_df)))

    # choose only reads that are mapped from at least start_pos_read to end_pos_read
    # TODO- different logic possible?
    #  -> currently simply split run to 3 sections
    # print(blast_df.describe())
    blast_df = blast_df[(blast_df.start_ref < min(recognized_positions)) & (blast_df.end_ref > max(recognized_positions))]
    blast_df = blast_df[['read', 'start_ref', 'end_ref']]
    print('len(blast_df): {}'.format(len(blast_df)))

    mutations_df = pd.read_csv(input_mutation_df)
    mutations_df = mutations_df[(mutations_df.ref != '-')] # remove insertions
    print('len(mutations_df): {}'.format(len(mutations_df)))

    # TODO- keep useful prints

    # drop reads containing a variation that is not the recognized mutation or the WT in the positions we are checking combinations for.
    # TODO-
    #  dropping too many reads. Issue is:
    #       1. any recognized_position with incorrect reference (happens mostly with deletions?)- rank0 is marked as mutation -> dropping most reads
    #       2. positions with un-declared rank2/rank3 too?
    #  Solution1: quick fix- remove problematic positions
    #  Solution1: real fix- run consensus until consolidation (fixing freq file wont work!)
    #  Solution2: add all ranks. but verify rest-of-code compatibility
    reads_with_chosen_mutations = mutations_df[(mutations_df.position.isin(recognized_positions)) & (mutations_df.full_mutation.isin(recognized_mutations))].read.unique()
    print('len(reads_with_chosen_mutations): {}'.format(len(reads_with_chosen_mutations)))

    reads_to_drop = mutations_df[(mutations_df.position.isin(recognized_positions)) & ~(mutations_df.full_mutation.isin(recognized_mutations))].read.tolist()
    # TODO- reads_to_drop = list(set(reads_to_drop))
    #  or simply add .read.unique().tolist() in above row
    print('len(set(reads_to_drop)): {}'.format(len(set(reads_to_drop))))

    # debug problematic positions
    reads_to_drop_with_details = mutations_df[(mutations_df.position.isin(recognized_positions)) & ~(mutations_df.full_mutation.isin(recognized_mutations))]
    problematic_positions = reads_to_drop_with_details[['position', 'read']].drop_duplicates()
    problematic_positions = problematic_positions.groupby('position').read.count().reset_index().rename(
        columns={'read': 'reads_to_drop_per_position'}).sort_values('reads_to_drop_per_position', ascending=False)
    print(problematic_positions.head())
    print('len(problematic_positions): {}'.format(len(problematic_positions)))

    blast_df = blast_df[~(blast_df.read.isin(reads_to_drop))]
    print('len(blast_df): {}'.format(len(blast_df)))

    # filter mutations_df to chosen mutations+reliable reads
    mutations_df = mutations_df[mutations_df.full_mutation.isin(recognized_mutations)]
    print('len(mutations_df): {}'.format(len(mutations_df)))
    # TODO- why merge? is this more efficient than [mutations_df.read.isin(blast_df.read)]?
    df = pd.merge(mutations_df, blast_df[['read']], how='right', on='read')
    print('len(df): {}'.format(len(df)))

    # filter mutations below freq threshold
    # simple frequencies for every mutation, keep only mutations over minimal mutation frequency (default 0)
    # TODO- convert for loop to pandas ops
    if minimal_mutation_frequency != 0:
        base_variants = []
        for m in df.full_mutation.drop_duplicates().dropna().tolist():
            total_read_count = len(df.read.drop_duplicates())
            reads_with_current_mutation = len(df[df.full_mutation == m].read.drop_duplicates())
            m_freq = (reads_with_current_mutation / total_read_count)
            if m_freq >= minimal_mutation_frequency:
                base_variants.append(m)

        mutations_df = mutations_df[mutations_df.full_mutation.isin(base_variants)]
        print('len(mutations_df): {}'.format(len(mutations_df)))
        df = pd.merge(mutations_df, blast_df[['read']], how='right', on='read')
        print('len(df): {}'.format(len(df)))

    else:
        base_variants = df.full_mutation.drop_duplicates().dropna().tolist()

    # produce mutations_per_read
    df = df.sort_values('position')
    df['full_mutation'] = df.full_mutation.astype(str)
    df['mutations_on_read'] = df.groupby('read')['full_mutation'].transform(', '.join)
    df_reads_with_mutations = df[['mutations_on_read', 'read']].drop_duplicates()

    print('duplicated reads?: {}'.format(df_reads_with_mutations['read'].duplicated().any()))

    strains_df = df_reads_with_mutations.groupby('mutations_on_read').read.count().reset_index().rename(columns={'read':'read_count'}).sort_values('read_count', ascending=False)
    strains_df['strain_frequency'] = strains_df.read_count / strains_df.read_count.sum()
    strains_df['mutations_on_strain'] = strains_df.mutations_on_read.str.replace('nan', 'WT')


    ####### now decide which strains are believable
    reliable_strains = []
    # classify strains within each base variant (each of the chosen mutation)
    for i in base_variants:
        print(i)
        strains_containing_variant = strains_df[strains_df.mutations_on_strain.str.contains(i)].copy().reset_index(drop=True)
        strains_containing_variant['base_variant'] = i
        strains_containing_variant['relative_frequency'] = strains_containing_variant.strain_frequency / strains_containing_variant.strain_frequency.sum()
        strains_containing_variant = strains_containing_variant.sort_values('relative_frequency', ascending=False)
        strains_containing_variant = choose_believable_strains(strains_containing_variant, substitution_error_cutoff, deletion_error_cutoff)
        reliable_strains.append(strains_containing_variant)

    # WT as base variant
    strains_containing_variant = strains_df.copy().reset_index(drop=True) #TODO- wt is sum of all base_variants... understand better later
    strains_containing_variant['base_variant'] = 'WT'
    strains_containing_variant['relative_frequency'] = strains_containing_variant.strain_frequency / strains_containing_variant.strain_frequency.sum()
    strains_containing_variant = strains_containing_variant.sort_values('relative_frequency', ascending=False)
    strains_containing_variant = choose_believable_strains(strains_containing_variant, substitution_error_cutoff, deletion_error_cutoff)
    reliable_strains.append(strains_containing_variant)

    reliable_strains = pd.concat(reliable_strains)
    reliable_strains.to_csv(output_file, index=False)

    # create results with only believable strains, and recalculate their relative frequency.
    df_reliable = reliable_strains[reliable_strains.believable == True][['mutations_on_strain', 'read_count']].drop_duplicates().copy()
    df_reliable['frequency'] = df_reliable.read_count / df_reliable.read_count.sum()
    df_reliable = df_reliable.rename(columns={'mutations_on_strain':'strain'})
    df_reliable.sort_values('frequency', ascending=False).to_csv(output_file + '.reliable.csv', index=False)
    return

# TODO- verify correctness after changes
# TODO- irrelevant? filtering is based distance is incorrect here?
def choose_believable_strains(df, substitution_error_cutoff, deletion_error_cutoff):
    df['believable'] = False
    df['closest_strain'] = None
    df['smallest_diff'] = None
    df.at[0, 'believable'] = True
    for i in range(1, len(df)):
        strain1 = df.at[i, 'mutations_on_strain']
        df1 = df[:i].copy()
        strains1 = df1[df1.believable == True].mutations_on_strain.tolist()
        closest_strain = None
        smallest_diff = None

        # find closest to strain1 from all strains so far
        for s in strains1:
            # list all diffs of: strain1 and current strain
            diffs = list(set.symmetric_difference(set(s.split(', ')), set(strain1.split(', '))))
            if 'WT' in diffs:
                diffs.remove('WT')
            if closest_strain == None or len(smallest_diff) > len(diffs):
                closest_strain = s
                smallest_diff = diffs
            elif len(smallest_diff) == len(diffs):
                deletions_in_smallest_diff = [i[-1] for i in smallest_diff].count('-')
                deletions_in_diff = [i[-1] for i in diffs].count('-')
                if deletions_in_smallest_diff >= deletions_in_diff:
                    smallest_diff = diffs

        df.at[i, 'closest_strain'] = closest_strain
        df.at[i, 'smallest_diff'] = smallest_diff

        # determine believability
        mutations_in_smallest_diff = len(smallest_diff)
        deletions_in_smallest_diff = [i[-1] for i in smallest_diff].count('-')
        if (( substitution_error_cutoff**(mutations_in_smallest_diff - deletions_in_smallest_diff)) * ( deletion_error_cutoff**(deletions_in_smallest_diff)) ) \
                <= df.at[i, 'relative_frequency']:
            df.at[i, 'believable'] = True

    return df


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-b", "--input_blast_df", type=str, help="path to blasts df csv", required=True)
    parser.add_argument('-m', '--input_mutation_df', type=str, help='path to mutations df csv', required=True)
    parser.add_argument('-p', '--input_chosen_mutations', type=str,
                        help='path to csv file with mutations to separate into strains. Every mutation should have its own row, the header row titled "variant", and the mutations should be written in the following format: "A1664.0G". The output file variants_chosen.csv from association_test_variant.py can be used here.',
                        required=True)
    parser.add_argument("-o", "--output_file", type=str, help="a path to an output file", required=True)
    parser.add_argument('-d', '--deletion_error_cutoff', type=float, required=True)
    parser.add_argument('-s', '--substitution_error_cutoff', type=float, required=True)
    parser.add_argument('-f', '--minimal_mutation_frequency', type=float, required=False, default=0,
                        help='frequency cutoff for a single mutation. Only mutations that are on the list and appear at least in this frequency in the population will be included in the strain analysis.')
    args = parser.parse_args()
    if not vars(args):
        parser.print_help()
        parser.exit(1)
    input_chosen_mutations = args.input_chosen_mutations
    input_blast_df = args.input_blast_df
    input_mutation_df = args.input_mutation_df
    minimal_mutation_frequency = args.minimal_mutation_frequency
    substitution_error_cutoff = args.substitution_error_cutoff
    deletion_error_cutoff = args.deletion_error_cutoff
    output_file = args.output_file
    count_haplotypes(input_chosen_mutations, input_blast_df, input_mutation_df, minimal_mutation_frequency,
                     substitution_error_cutoff, deletion_error_cutoff, output_file)


def test_main():
    samples = ['env3']
    # samples = ['env10', 'env11', 'env12', 'env13', 'env14', 'env15', 'env4', 'env5', 'env6','env7', 'env8', 'env9']
    for sample in samples:
        print('Sample: {}'.format(sample))

        sample_dir = '/Volumes/STERNADILABHOME$/volume1/shared/analysis/HIV_shafer/associvar/{}'.format(sample)
        run_name = 'contig1_freq1_all_ranks_valid_ref_only'
        # run_name = 'contig1'

        count_haplotypes('{}/chosen_mutations_{}.csv'.format(sample_dir, run_name),
                         '{}/blasts.csv'.format(sample_dir),
                         '{}/mutations.csv'.format(sample_dir),
                         0,
                         0,
                         0,
                         '{}/associvar_strain_{}_{}.csv'.format(sample_dir, sample, run_name))


if __name__ == "__main__":
    # TODO- uncomment after test runs
    main()

    # TODO- tmp for test runs
    # test_main()
