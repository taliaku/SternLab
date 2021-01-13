import os

import pandas as pd
import glob
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

#from RG_HIVC_analysis import constants
#from RG_HIVC_analysis.constants import gag_ET86_interval, pol_ET86_interval, env_ET86_interval, orig_excluded_samples, control_excluded_patients
from scripts.diversity_analysis import pi_diversity_calc, apply_pi_related_filters

sns.set_context("poster")


def get_simple_diversity_stats():
    freq_files = glob.glob('/Users/omer/PycharmProjects/SternLab/RG_data_analysis/freq_files/*/*.freqs')
    basic_muts_count_stats = pd.DataFrame(
        columns=['sample_id', 'major_subs_count'])

    i = 0
    for file in freq_files:
        sample_id = str(file).split("/")[7]
        major_subs_count = count_major_subs(file)

        row = [sample_id] + [major_subs_count]
        # print(row)
        basic_muts_count_stats.loc[i] = row
        i = i + 1

    basic_muts_count_stats = basic_muts_count_stats.sort_values(by='major_subs_count', ascending=False)
    return basic_muts_count_stats


def count_major_subs(freq_file):
    df = pd.read_csv(freq_file, sep='\t')
    subs_count = df['Pos'].loc[(df['Read_count'] > 100) & (df['Base'] != df['Ref']) & (df['Rank'] == 0) & (df['Base'] != '-') & (df['Ref'] != '-') ].count()
    return subs_count


def generate_pi_rates_summary():
    """Deprecated"""

    # freq_files = glob.glob('/Users/omer/PycharmProjects/SternLab/RG_data_analysis/freq_files_ZA04_2/*')
    freq_files = glob.glob('/Users/omer/PycharmProjects/SternLab/RG_HIVC_analysis/ET86_2s/*.freqs')
    pi_diversity_rates = pd.DataFrame(
        columns=['sample_id', 'global', 'gag', 'pol', 'env'])

    # TODO- use append instead of loc[i]
    i = 0
    for file in freq_files:
        sample_id = os.path.splitext(os.path.basename(file))[0]
        if sample_id in orig_excluded_samples:
            print('Excluded sample: {} - Skipping'.format(sample_id))
            continue
        print('Handling sample: ' + sample_id)

        freq_df = pd.read_csv(file, sep='\t')
        filtered_freq_df = apply_pi_related_filters(freq_df, frequency_threshold= 0.001, coverage_threshold= 100)
        global_pi_rate = pi_diversity_calc(data=filtered_freq_df)

        # hxb2_file = glob.glob('/Users/omer/PycharmProjects/SternLab/RG_data_analysis/freq_files_HXB2_2/{}/*.freqs'.format(sample_id))[0]
        # freq_df_hxb2 = pd.read_csv(hxb2_file, sep='\t')
        # et86_file = glob.glob('/Users/omer/PycharmProjects/SternLab/RG_data_analysis/ET86_2s/{}.freqs'.format(sample_id))[0]
        # freq_df_et86 = pd.read_csv(et86_file, sep='\t')
        # global_pi_rate2 = pi_diversity_calc(data=freq_df_et86, min_read_count= 1000, freq_threshold= 0)

        gag_pi_rate = pi_diversity_calc(data=filtered_freq_df, interval= gag_ET86_interval)
        pol_pi_rate = pi_diversity_calc(data=filtered_freq_df, interval= pol_ET86_interval)
        env_pi_rate = pi_diversity_calc(data=filtered_freq_df, interval= env_ET86_interval)

        row = [sample_id] + [global_pi_rate] + [gag_pi_rate] + [pol_pi_rate] + [env_pi_rate]
        # print(row)
        pi_diversity_rates.loc[i] = row
        i = i + 1

    pi_diversity_rates = pi_diversity_rates.sort_values(by='global', ascending=False)
    print(pi_diversity_rates)
    pi_diversity_rates.to_csv('/Users/omer/PycharmProjects/SternLab/RG_HIVC_analysis/output_tables/pi_rates_ET86_2.csv', index=False)
    return pi_diversity_rates


def plot_diversity_by_time():
    """Deprecated"""

    pd.set_option('display.width', 600)  # TODO- remove
    pd.set_option('display.max_columns', 16)  # TODO- remove

    # joining diversity values with patient & date info
    samples_to_patient_and_dates = pd.read_csv('/Users/omer/PycharmProjects/SternLab/RG_HIVC_analysis/output_tables/final_ZA04.csv',
                                sep=',')[['sample_id', 'ind_id', 'sample_date']]
    samples_to_patient_and_dates = samples_to_patient_and_dates.set_index('sample_id')
    pi_rates = pd.read_csv('/Users/omer/PycharmProjects/SternLab/RG_HIVC_analysis/output_tables/pi_rates_ET86_2_median.csv', sep=',').set_index('sample_id')

    pis_by_ind = samples_to_patient_and_dates.join(pi_rates)

    # sorting by patient + sampling date
    pis_by_ind['sample_date'] = pd.to_datetime(pis_by_ind['sample_date'], format='%d/%m/%Y')
    ordered_pis_by_ind = pis_by_ind.sort_values(by=['ind_id', 'sample_date'])

    # coverting "sample_date" to "years_since_infection"
    first_samples_dates = ordered_pis_by_ind.groupby('ind_id').first().reset_index()
    first_samples_dates = first_samples_dates[['ind_id', 'sample_date']]
    # print(first_samples_dates)

    ordered_pis_by_ind = ordered_pis_by_ind.merge(first_samples_dates,
                                        on='ind_id',
                                        how='left',
                                        sort=False,
                                        suffixes= ('','_r'))
    ordered_pis_by_ind['years_since_infection'] = (ordered_pis_by_ind['sample_date'] - ordered_pis_by_ind['sample_date_r']) / np.timedelta64(1, 'Y')
    print(ordered_pis_by_ind)

    # generating plot
    ordered_pis_by_ind = ordered_pis_by_ind.melt(id_vars= ('ind_id', 'years_since_infection'),
                        value_vars= ('global', 'gag', 'pol', 'env'),
                        var_name='regions',  value_name='pi_diversity'
                        )
    ordered_pis_by_ind = ordered_pis_by_ind.sort_values(by=['ind_id', 'years_since_infection'])
    # print(ordered_pis_by_ind[ordered_pis_by_ind['ind_id'] == 16207])

    g = sns.relplot(
        x='years_since_infection',
        y='pi_diversity',
        col='ind_id',
        hue='regions',
        data=ordered_pis_by_ind,
        col_wrap=5,
        kind='line',
        facet_kws={'sharex': True, 'legend_out':True},
        )
    g.set(yscale="log")
    g.set_ylabels("Pi diversity")
    g.set_xlabels("ET first sample (Years)")
    # g.set_xticklabels(rotation=45, fontsize=14)

    # extracting plot
    plt.show()
    # plt.savefig(fname= '/Users/omer/PycharmProjects/SternLab/RG_HIVC_analysis/figures/pi_trends_ET86_2.pdf')
    # g.savefig('/Users/omer/PycharmProjects/SternLab/RG_HIVC_analysis/figures/pi_rates_ET86_2.png')

def aggregation_attempts():
    """Deprecated"""

    summary_table = pd.read_csv('/Users/omer/PycharmProjects/SternLab/RG_data_analysis/output_tables/final_ZA04.csv',
                                               sep=',').set_index('sample_id')
    # a= summary_table.groupby('ind_id')['sample_date'].aggregate(min).unstack().reset_index()
    a = summary_table[['ind_id', 'sample_date']].groupby('ind_id').agg(lambda x: x.iloc[0]).set_index('ind_id')
    print(a)
    sfil= summary_table[['ind_id', 'sample_date', 'pi_diversity']].set_index('ind_id')
    join= sfil.join(a, lsuffix='sample_date', rsuffix='fisrt_date')
    print(join)

    # # print(summary_table)
    # # a = summary_table[summary_table["ind_id"] == 12796][['ind_id', 'sample_date', 'pi_diversity']]
    # a = summary_table[['ind_id', 'sample_date', 'pi_diversity']]
    # fd = a[['ind_id', 'sample_date']].groupby('ind_id').agg(lambda x: x.iloc[0])
    # print(a)
    # print(fd)


    # diversity_trends2 = summary_table.groupby('ind_id').apply(list)
    # diversity_trends_x = summary_table.groupby('ind_id')['sample_date'].apply(list)
    # # diversity_trends = summary_table.groupby('ind_id').agg({'sample_date':'list','pi_diversity':'sum'})


def pi_rates_generate_and_plot_v2():
    # get freqs raw data
    run_folder = 'orig_high'
    unified_freq_df = pd.read_csv('/Users/omer/PycharmProjects/SternLab/RG_HIVC_analysis/runs/{}/unified_freqs_filtered_verbose.csv'.format(run_folder))

    # generate pi rate values (global & regionals)
    unified_freq_df = apply_pi_related_filters(unified_freq_df, frequency_threshold=constants.freq_threshold, coverage_threshold=constants.coverage_threshold)
    pi_rates_by_sample = pi_diversity_calc(data=unified_freq_df, pivot_cols= ['ind_id', 'sample_id', 'years_since_infection'])

    # gag_pi_rates_by_sample = pi_diversity_calc(data=unified_freq_df, pivot_cols= ['sample_id'], interval=gag_ET86_interval)
    # pi_rates_by_sample = pi_rates_by_sample.merge(gag_pi_rates_by_sample, on='sample_id', how='left', sort=False, suffixes=('', '_gag'))
    pol_pi_rates_by_sample = pi_diversity_calc(data=unified_freq_df, pivot_cols= ['sample_id'], interval=pol_ET86_interval)
    pi_rates_by_sample = pi_rates_by_sample.merge(pol_pi_rates_by_sample, on='sample_id', how='left', sort=False, suffixes=('', '_pol'))
    # env_pi_rates_by_sample = pi_diversity_calc(data=unified_freq_df, pivot_cols= ['sample_id'], interval=env_ET86_interval)
    # pi_rates_by_sample = pi_rates_by_sample.merge(env_pi_rates_by_sample, on='sample_id', how='left', sort=False, suffixes=('', '_env'))

    # TODO- understand double plotting
    # input VL values
    # samples_vl = pd.read_csv('/Users/omer/PycharmProjects/SternLab/RG_HIVC_analysis/output_tables/sample_vl.csv')
    #merge
    # pi_rates_by_sample = pi_rates_by_sample.merge(samples_vl, on='sample_id', how='left', sort=False, suffixes=('', ''))

    # plot
    print('plotting')
    # Explanation: by definition of the pi measure- plotting mean value (weighted mean & median are irrelevant), with no interest in confidence interval.
    # TODO- check this understanding ^ with maoz

    pi_rates_by_sample = pi_rates_by_sample.melt(id_vars= ('ind_id', 'years_since_infection'),
                        # value_vars= ('Pi', 'Pi_gag', 'Pi_pol', 'Pi_env'),
                        value_vars= ('Pi', 'Pi_pol'),
                        var_name='regions',  value_name='pi_diversity'
                        )
    g = sns.relplot(
        x='years_since_infection',
        y='pi_diversity',
        col='ind_id',
        hue='regions',
        data=pi_rates_by_sample,
        col_wrap=4,
        kind='line',
        facet_kws={'sharex': True, 'legend_out': True},
    )
    g.set(yscale="log")
    plt.ylim(8*10**-4, 5*10**-2)
    g.set_ylabels("Pi diversity")
    g.set_xlabels("ET first sample (Years)")
    # g.set_xticklabels(rotation=45, fontsize=14)

    # extract plot
    plt.show()
    # g.savefig('/Users/omer/PycharmProjects/SternLab/RG_HIVC_analysis/figures/pi_rates_ET86_{}_minimized2.png'.format(run_folder))

def plot_vl_rates():
    # input VL values
    samples_vl_dsi = pd.read_csv('/Users/omer/PycharmProjects/SternLab/RG_HIVC_analysis/output_tables/samples_to_patient_and_dsi.csv', sep='\t')

    # plot
    g = sns.relplot(
        x='days_since_infection',
        y='viral_load',
        col='patient',
        data=samples_vl_dsi,
        col_wrap=5,
        kind='line',
        facet_kws={'sharex': True, 'legend_out': True},
    )
    g.set(yscale="log")
    # plt.ylim(8*10**-4, 3*10**-2)
    g.set_ylabels("Viral load")
    g.set_xlabels("ET first sample (Days)")
    # g.set_xticklabels(rotation=45, fontsize=14)

    # extract plot
    plt.show()
    # g.savefig('/Users/omer/PycharmProjects/SternLab/RG_HIVC_analysis/figures/vl_rates_ET86_4s.png')


def pi_rates_generate_and_plot_v1():
    generate_pi_rates_summary()
    plot_diversity_by_time()


def add_pi_rates_to_unified_freqs_df():
    run_folder = 'orig_high'
    unified_freq_df = pd.read_csv('/Users/omer/PycharmProjects/SternLab/RG_HIVC_analysis/runs/{}/unified_freqs_filtered_verbose.csv'.format(run_folder))

    # filters?
    unified_freq_df = apply_pi_related_filters(unified_freq_df, frequency_threshold=0, coverage_threshold=0)
    samples = ['83476_S27','504193_S37','100888_S14','TASPX100711_S29','504212_S56','X160138_S81','504224_S68'] # 26892 samples
    unified_freq_df = unified_freq_df[unified_freq_df['sample_id'].isin(samples)]

    pi_rates_by_pos = pi_diversity_calc(data=unified_freq_df, pivot_cols=['ind_id', 'sample_id'])
    print(pi_rates_by_pos)

    print('plotting freq_plot')
    g = sns.relplot(x="Pos",
                    y="Pi",
                    col='sample_id',
                    # col_order='years_since_infection', # chronological presentation
                    # # hue='sample_id',
                    # hue='mutation_type',
                    col_wrap=7,
                    # kind="line",
                    # join=True,
                    data=pi_rates_by_pos)

    # plot adjustments
    # g.set(yscale="log")
    plt.ylim(5e-4, 1)
    # g.fig.suptitle(plot_header, y=0.1)
    # g.set_ylabels("mutation_rate")
    # g.set_xlabels("ET first sample (Years)")
    # g.set_xticklabels(rotation=45, fontsize=11)

    # extract plot
    plt.show()
    # plt.savefig(fname=fname)
    # g.savefig('')


if __name__ == "__main__":
    # pi
    # pi_rates_generate_and_plot_v1()
    pi_rates_generate_and_plot_v2()
    # add_pi_rates_to_unified_freqs_df()

    # vl
    # plot_vl_rates()
