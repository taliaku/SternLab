import sys

import pandas as pd
import glob
import os

from RG_HIVC_analysis.constants import RT_short_ET86_interval, gag_ET86_interval, pol_ET86_interval, env_ET86_interval


def create_coverage_stats_table(run_direcory='/sternadi/home/volume1/shared/analysis/HIV_ravi_gupta/run3'):
    freq_files = glob.glob('%s/*.freqs' % run_direcory)
    # freq_files = glob.glob('%s/**/*.freqs' % run_direcory)
    # coverage_stats = pd.DataFrame(columns=['sample_id', 'coverage_median', 'long_coverage', 'gag_cov', 'pol_cov', 'env_cov', 'file_size'])
    coverage_stats = pd.DataFrame(columns=['sample_id', 'RT_cov_median'])

    i=0
    for file in freq_files:
        sample_id = os.path.splitext(os.path.basename(file))[0]
        # cov_median = extract_coverage_median(file)
        RT_cov_median = extract_coverage_median(file, interval=RT_short_ET86_interval)

        # long_cov, gag_cov, pol_cov, env_cov = extract_long_cov(file)
        # file_size = os.path.getsize("/sternadi/datasets/volume2/HIV_ravi_gupta_processed/"+sample_id+"/"+sample_id+"_R1.fastq.gz")
        # file_size += os.path.getsize(
        #     "/sternadi/datasets/volume2/HIV_ravi_gupta_processed/" + sample_id + "/" + sample_id + "_R2.fastq.gz")
        # file_size /= 1024

        # row = [sample_id] + [cov_median] + [long_cov] + [gag_cov] + [pol_cov] + [env_cov] + [file_size]
        row = [sample_id] + [RT_cov_median]
        # print(row)
        coverage_stats.loc[i] = row
        i=i+1

    # coverage_stats = coverage_stats.sort_values(by='coverage_median', ascending=False)
    return coverage_stats

def extract_coverage_median(freq_file, interval = (0, sys.maxsize)):
    df = pd.read_csv(freq_file, sep='\t')
    df_filtered = df.loc[df['Read_count'] > 100]
    df_filtered = df_filtered[(df_filtered["Pos"] >= interval[0]) & (df_filtered["Pos"] <= interval[1])]
    cov_median = df_filtered['Read_count'].median()
    return cov_median

def extract_long_cov(freq_file):
    df = pd.read_csv(freq_file, sep='\t')
    df = df.drop_duplicates("Pos")
    df = df.loc[df['Read_count'] > 1000]
    long_cov = df['Pos'].count()

    gag_cov = df['Pos'].loc[(df['Pos'] > gag_ET86_interval[0]) & (df['Pos'] < gag_ET86_interval[1])].count()
    pol_cov = df['Pos'].loc[(df['Pos'] > pol_ET86_interval[0]) & (df['Pos'] < pol_ET86_interval[1])].count()
    env_cov = df['Pos'].loc[(df['Pos'] > env_ET86_interval[0]) & (df['Pos'] < env_ET86_interval[1])].count()
    return long_cov, gag_cov, pol_cov, env_cov


def compare_coverage_mad_increase():
    freq_files_run2 = glob.glob('/sternadi/home/volume1/shared/analysis/HIV_ravi_gupta/run2/**/*.freqs')
    medians = pd.DataFrame(columns=['sample_id', 'coverage_run2', 'coverage_run3', 'delta_percent'])

    i=0
    for file_run2 in freq_files_run2:
        sample_id = str(file_run2).split("/")[8]
        cov_mad2 = extract_cov_mad(file_run2)
        file_run3 = glob.glob(file_run2.replace('run2','run3'))[0]
        cov_mad3 = extract_cov_mad(file_run3)
        delta_percent = (cov_mad3/cov_mad2 - 1) * 100

        row = [sample_id] + [cov_mad2] + [cov_mad3] + [delta_percent]
        print(row)
        # medians.loc[i] = row
        i=i+1

    medians.sort_values(by='coverage_run3', ascending=False)

def extract_cov_mad(freq_file):
    df = pd.read_csv(freq_file, sep='\t')
    df_filtered = df.loc[df['Read_count'] > 100]
    cov_mad = df_filtered['Read_count'].mad()
    return cov_mad



def create_general_stats_table():
    summary_files = glob.glob('/sternadi/home/volume1/shared/analysis/HIV_ravi_gupta/run3/**/pipeline_summary.txt')
    general_stats = pd.DataFrame(columns=['sample_id', 'reads_mapped_to_ref', 'reads_in_freq_count', 'reads_mapped_once', 'bases_called'])

    i = 0
    for summary_file in summary_files:
        sample_id = str(summary_file).split("/")[8]
        with open(summary_file, 'r') as file:
            content = file.read().split()
            reads_mapped_to_ref = int(content[content.index('reference:')+1])
            reads_in_freq_count = int(content[content.index('count:')+1])
            reads_mapped_once = int(content[content.index('once:')+1])
            bases_called = int(content[content.index('called:')+1])


        row = [sample_id] + [reads_mapped_to_ref] + [reads_in_freq_count] + [reads_mapped_once] + [bases_called]
        # print(row)
        general_stats.loc[i] = row
        i = i + 1

    general_stats = general_stats.sort_values(by='reads_mapped_to_ref', ascending=False)
    return general_stats

def mergre_summary_tables():
    dates_vl = pd.read_csv('/Users/omer/PycharmProjects/SternLab/RG_data_analysis/dates_vl_stats.csv', sep=',')
    samples_format_conversion = pd.read_csv('/Users/omer/PycharmProjects/SternLab/RG_data_analysis/samples_format_conversion.csv', sep=',')
    pi_rates = pd.read_csv('/Users/omer/PycharmProjects/SternLab/RG_data_analysis/pi_rates_ZA04_2.csv', sep=',')
    coverage_stats = pd.read_csv('/Users/omer/PycharmProjects/SternLab/RG_data_analysis/coverage_stats_ET86_2.csv', sep=',')

    dates = dates_vl.set_index('FASTQ_name')
    conv = samples_format_conversion.set_index('table_fastq')
    join1 = dates.join(conv).set_index('sample_id')
    pi_rates = pi_rates.set_index('sample_id')
    cov = coverage_stats.set_index('sample_id')
    join2 = pi_rates.join(cov)

    final = join1.join(join2)
    final['sample_date'] = pd.to_datetime(final['sample_date'], format='%d/%m/%Y')
    final = final.sort_values(by=['ind_id', 'sample_date'])

    print(final)
    final.to_csv(path_or_buf='/Users/omer/PycharmProjects/SternLab/RG_data_analysis/final_ET86_pol_cov.csv')

    print(final.loc['130945_S2'])


def get_coverage_distribution_from_unified_freq_df():
    unified_freq_df_2s = pd.read_csv('/Users/omer/PycharmProjects/SternLab/RG_HIVC_analysis/ET86_2s/unified_freqs_filtered_verbose.csv')
    unified_freq_df_2s_muts = pd.read_csv('/Users/omer/PycharmProjects/SternLab/RG_HIVC_analysis/ET86_2s/with_muts_coordinated/unified_freqs_filtered_verbose.csv')
    unified_freq_df_4s = pd.read_csv('/Users/omer/PycharmProjects/SternLab/RG_HIVC_analysis/ET86_4s/unified_freqs_filtered_verbose.csv')
    unified_freq_df_4s_muts = pd.read_csv('/Users/omer/PycharmProjects/SternLab/RG_HIVC_analysis/ET86_4s/with_muts_and_con_as_ref/unified_freqs_filtered_verbose.csv')

    freq_series = [unified_freq_df_2s, unified_freq_df_2s_muts, unified_freq_df_4s, unified_freq_df_4s_muts]

    for f in freq_series:
        print(f.shape)
        f = f[(f["Pos"] >= RT_short_ET86_interval[0]) & (f["Pos"] <= RT_short_ET86_interval[1])]
        coverage_median = f.groupby('ind_id').mean()['Read_count']
        print(coverage_median)

        # g = sns.distplot(f['Read_count'])
        # g = sns.boxplot(x='ind_id', y='Read_count', data=f)
        #
        # plt.show()

    # df_final = reduce(lambda left, right: pd.merge(left, right, on='ind_id'), freq_series)
    # print(df_final)

    # b = unified_freq_df[unified_freq_df["Read_count"] > 1000]
    # long_cov = b.groupby('ind_id').count()['Read_count']
    # print(long_cov)


def get_general_df_stats():
    freq_df = pd.read_csv('/Users/omer/PycharmProjects/SternLab/RG_HIVC_analysis/ET86_4s/unified_freqs_filtered_verbose.csv')

    # pd.set_option('display.width', 600)
    pd.set_option('display.max_columns', 16)

    print(freq_df.describe(include='all'))


if __name__ == "__main__":
    # parser = OptionParser("usage: %prog [options]\nTry running %prog --help for more information")
    # parser.add_option("-d", "--run_dir", dest="run_dir", help="rundir")
    # parser.add_option("-n", "--run_name", dest="run_name", help="run name")
    # (options, args) = parser.parse_args()
    # run_dir = options.run_dir
    # run_name = options.run_name

    run_name = 'orig_low'
    run_dir = '/Users/omer/PycharmProjects/SternLab/RG_HIVC_analysis/runs/' + run_name

    coverage_stats = create_coverage_stats_table(run_direcory=run_dir)

    cov_stats_dir = '/Users/omer/PycharmProjects/SternLab/RG_HIVC_analysis/output_tables/cov_stats'
    coverage_stats.to_csv(path_or_buf=('%s/%s.csv' % (cov_stats_dir, run_name)), index=False)

    # get_coverage_distribution_from_unified_freq_df()
    # get_general_df_stats()


def verify_chosen_stats(unified_freq_df):
    pd.set_option('display.max_rows', 250)
    chosen_stat = unified_freq_df[(unified_freq_df['ind_id'] == 22097) & (unified_freq_df['Mutation_type'] == 'synonymous')]
    # The mut measure is defined as rate per position, that's why avreage on pos is ok.
    # aggregation on all combinations of positions+muts able to generate specific kind of mutation
    # what if there's only one pair pos+muts that can generate stop mutation? its value is the mut value? and if there's 2?
    #          -> final value we generate is a repesentitive of what?
    # total coverage collected on those positions- serves as sanity measure
    chosen_stat = chosen_stat[chosen_stat['sample_id'] == '504188_S32']
    # chosen_stat = chosen_stat[chosen_stat['sample_id'] == '504211_S55']
    chosen_stat = chosen_stat[['Freq', 'Read_count']]
    chosen_stat['count_for_position'] = chosen_stat['Freq'] * chosen_stat['Read_count']

    print(chosen_stat)

    print('mean is {}'.format(chosen_stat['Freq'].mean()))
    print('weighted_mean is {}'.format(chosen_stat['count_for_position'].sum()/ chosen_stat['Read_count'].sum()))


def potential_analysis(unified_freq_df):
    pd.set_option('display.width', 600)
    pd.set_option('display.max_columns', 16)
    pd.set_option('display.max_rows', 250)
    potential_analysis = unified_freq_df.groupby(['sample_id', 'Mutation_type'])['Pos'].agg(['count'])
    print(potential_analysis)


