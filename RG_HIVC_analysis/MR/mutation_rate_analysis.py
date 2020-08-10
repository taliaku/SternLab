import glob

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns;
from scipy import stats

from RG_HIVC_analysis.constants import orig_patients_ordered_by_mr, zn_patients_ordered_by_mr, get_ET86_region
from RG_HIVC_analysis.data_adaptations import generate_unified_filtered_verbose_freqs_df

sns.set_context("poster")

mr_dir = '/Users/omer/PycharmProjects/SternLab/RG_HIVC_analysis/MR/'
fits_summaries_dir = mr_dir + 'fits_runs_summaries/'
MR_MEAN_OF_MEANS = 1.81E-05 # arbitrary threshold


def plot_mutation_rate_distribution(unified_freq_df):
    pd.set_option('display.width', 600)  # TODO- remove
    pd.set_option('display.max_columns', 16)  # TODO- remove

    # g = sns.catplot(
    #     x='years_since_infection',
    #     y='Freq',
    #     col='ind_id',
    #     hue='Mutation_type',
    #     data=unified_freq_df,
    #     col_wrap=5,
    #     kind='box',
    #     facet_kws={'sharex': False, 'legend_out':True},
    #     )

    g = sns.boxplot(x='ind_id', y='Freq',hue='Mutation_type', data=unified_freq_df)

    # plot adjustments
    # g.set(yscale="log")
    plt.yscale('log')
    # plt.ylim(10**-6, 1)
    # g.set_ylabels("Muts rate")
    # g.set_xlabels("ET first sample (Years)")
    # g.set_xticklabels(rotation=45, fontsize=14)

    # extract plot
    plt.show()
    # plt.savefig(fname= '/Users/omer/PycharmProjects/SternLab/RG_HIVC_analysis/figures/muts_rates_distribution.png')
    # g.savefig('')


def plot_mutation_rate_conf_interval(unified_freq_df):
    add_weighted_freq(unified_freq_df)

    def weighted_varaint(x, **kws):
        var, count = map(np.asarray, zip(*x))
        return var.sum() / count.sum()

    # basic presentation
    # g = sns.factorplot(x="sample_id",
    #                    y="count_and_weight",
    #                    data=unified_freq_df,
    #                    hue="Mutation_type",
    #                    hue_order=["missense", "synonymous", "stop"],
    #                    # palette="tab20",
    #                    join=False,
    #                    orient="v",
    #                    estimator=weighted_varaint,
    #                    dodge=1)

    # advanced, chronological presentation
    g = sns.catplot(x= "years_since_infection",
                       y= "count_and_weight",
                       hue="Mutation_type",
                       hue_order=["missense", "synonymous", "stop"],
                       # palette="tab20",
                       col='ind_id',
                       col_wrap=5,
                       n_boot=1000,
                       # join=True,
                       orient="v",
                       kind='point',
                       estimator=weighted_varaint,
                       dodge=1,
                       data=unified_freq_df)

    # TODO- optional:
    # gateway to value extraction
    # ax = sns.pointplot(x= "ind_id",
    #                    y= "count_and_weight",
    #                    data=unified_freq_df,
    #                    hue="Mutation_type",
    #                    # palette="tab20",
    #                    join=False,
    #                    orient="v",
    #                    estimator=weighted_varaint,
    #                    dodge=1)


    # plot adjustments
    g.set(yscale="log")
    # plt.ylim(10**-6, 1)
    g.set_ylabels("mutation_rate")
    # g.set_xlabels("ET first sample (Years)")
    # g.set_xticklabels(rotation=45, fontsize=11)

    # extract plot
    # plt.savefig(fname= '/Users/omer/PycharmProjects/SternLab/RG_HIVC_analysis/figures/mutation_rates_ET86_4s.png')
    plt.show()
    # g.savefig('')


def add_weighted_freq(unified_freq_df):
    unified_freq_df["count_for_position"] = unified_freq_df["Freq"] * unified_freq_df["Read_count"]
    # unified_freq_df["count_for_position"] = np.where(unified_freq_df["Prob"] < 0.95, 0, unified_freq_df["no_variants"]) #TODO- relevant?
    unified_freq_df["count_and_weight"] = list(zip(unified_freq_df.count_for_position, unified_freq_df.Read_count))

    # unified_freq_df.to_csv('/Users/omer/PycharmProjects/SternLab/RG_HIVC_analysis/ET86_2s/with_muts_coordinated/examine_freqs_with_weights.csv', index=False)
    print(unified_freq_df.shape)


def calc_regression_lines(unified_freq_df):
    # calc slopes etc. for all patients
    regression_coeefs_summary = pd.DataFrame(columns=['ind_id','mut_type','slope', 'p_value'])

    # quick an dirty, but precise
    for mut_type in ["synonymous", "missense"]:
        for ind_id in ['12796','13003','15664','16207','17339','19937','22097','22763','22828','23271','26892','28545','28841','29219','29447','31254','34253','47939']:
            regressed_data = unified_freq_df[(unified_freq_df['ind_id'] == int(ind_id)) & (unified_freq_df['Mutation_type'] == mut_type)]

            regressed_data["count_for_position"] = regressed_data["Freq"] * regressed_data["Read_count"]
            regressed_data = regressed_data.groupby(['years_since_infection'])['count_for_position', 'Read_count'].agg(['sum'])
            regressed_data.columns = ['count_for_position_sum', 'Read_count_sum']
            regressed_data = regressed_data.reset_index()
            regressed_data['weighted_freq'] = regressed_data['count_for_position_sum'] / regressed_data['Read_count_sum']

            # get coeffs of linear fit
            slope, intercept, r_value, p_value, std_err = stats.linregress(regressed_data['years_since_infection'],
                                                                           regressed_data['weighted_freq'])

            regression_coeefs_summary = regression_coeefs_summary.append({'ind_id': ind_id, 'mut_type': mut_type, 'slope': "{:.2e}".format(slope), 'p_value': "{0:0.4f}".format(p_value)}, ignore_index=True)
            print(ind_id+','+mut_type+': '+ "{:.2e}".format(slope))

    # export to summary file
    regression_coeefs_summary.to_csv('/Users/omer/PycharmProjects/SternLab/RG_HIVC_analysis/output_tables/mutation_rates_regression_stats.csv')


    # regression attempt
    # g = sns.lmplot(x= "years_since_infection",
    #                y= "count_and_weight",
    #                hue="Mutation_type",
    #                hue_order=["missense", "synonymous", "stop"],
    #                col='ind_id',
    #                col_wrap=5,
    #                n_boot=100,
    #                # join=True,
    #                # kind='point',
    #                x_estimator=weighted_varaint,
    #                truncate=True,
    #                logx= True,
    #                data=unified_freq_df)

## DEPRECATED ##
def analyze_mutation_rates_by_fits():
    run_name = 'orig_high'
    patients = ['26892']

    # summarize all patients
    df = pd.DataFrame()
    for patient in patients:
        # get mr_per_pos_by_fits
        fits_input_files_dir = '/Users/omer/PycharmProjects/SternLab/RG_HIVC_analysis/runs/{}/fits/'.format(run_name)
        patient_dir = fits_input_files_dir + '{}/'.format(patient)

        p_mu_dict = {}
        for mut in ['GA', 'AG', 'CT', 'TC']:
            summary_average_per_mut = get_average_mr_after_quantile(patient_dir, mut, patient)
            print('Average MR for mutation - {}: {}'.format(mut, summary_average_per_mut))
            p_mu_dict[mut] = summary_average_per_mut

        if patient == patients[0]:
            df = pd.DataFrame(p_mu_dict)
        else:
            df = pd.merge(df, pd.DataFrame(p_mu_dict), how='left', left_index=True, right_index=True,
                          suffixes=('', '_' + patient))

    # extract
    df.to_csv(fits_input_files_dir + 'mutation_rate_summary_fits.csv')

    # TODO- run also on all syn positions (without conservation filter) & then do inter-quartile range filtering (instead of only high range)


def get_average_mr_after_quantile(input_files_path, mut, patient):
    mr_median_list = pd.read_csv(input_files_path + '{}_median_summary.txt'.format(mut), sep='\t')
    # filter out 4th quartile MR
    # visualize af distribution
    ax = sns.distplot(mr_median_list['mr_per_pos'])
    plot_header = '{} MR distribution across pos by fits - patient {}'.format(patient, mut)
    ax.set_title(plot_header)
    plt.show()
    # plt.savefig(input_files_path + '/' + plot_header + '.png')
    plt.cla()
    # filter
    q_thresh = mr_median_list.quantile(.75)['mr_per_pos']
    median_list_mid_low_range = mr_median_list[mr_median_list['mr_per_pos'] < q_thresh]
    # get mean MR per patient, per mutation
    summary_average_per_mut = median_list_mid_low_range.mean()['mr_per_pos']
    return summary_average_per_mut

def plot_fits_mr_genomewide():
    # run_name = 'orig_high'
    # plots_order = orig_patients_ordered_by_mr
    run_name = 'zn'
    plots_order = zn_patients_ordered_by_mr

    mr_data = pd.read_csv(fits_summaries_dir + '%s/fits_mr_summary_%s.csv' % (run_name, run_name) )
    mr_data = mr_data[mr_data['significance'] == 'significant']
    # mr_data = mr_data[mr_data['Patient'] == 26892]

    # analyse only suspicious positions
    analyse_only_high_pos = False
    if analyse_only_high_pos:
        mr_data = mr_data[mr_data['MR'] > MR_MEAN_OF_MEANS]

    g1 = sns.relplot(x="Pos",
                    y="MR",
                    col='Patient',
                    col_order= plots_order,
                    # hue='Mutation',
                    color='firebrick',
                    col_wrap=7,
                    kind='line',
                    data=mr_data)

    g1.set(yscale="log")
    plt.ylim(1e-8, 1e-1)
    plt.xlim(0, 9026)
    # plt.show()
    plt.savefig(fname=mr_dir + 'genomewide_mr_%s.png' % run_name)

def plot_mr_hotspots():
    run_name = 'orig_high'
    plots_order = orig_patients_ordered_by_mr
    # run_name = 'zn'
    # plots_order = zn_patients_ordered_by_mr
    generate_outputs = False

    mr_data = pd.read_csv(fits_summaries_dir + '%s/fits_mr_summary_%s.csv' % (run_name, run_name) )

    # filters
    mr_data = mr_data[mr_data['significance'] == 'significant']

    # threshold - deprecated (improved logic below)
    # mr_data = mr_data[mr_data['MR'] > mr_mean_of_means]

    # generate sliding window
    mr_data_rolling = mr_data[['Patient','Pos','MR']]
    # TODO- averaging on log instead of mean- makes more sense indeed?
    mr_data_rolling['log_MR'] = np.log(mr_data_rolling['MR'])
    mr_data_rolling['Pos_as_date'] = pd.to_datetime(mr_data_rolling['Pos'], unit='s')

    dfs = []
    window_size = '10'
    minimum_pos = '2'
    for p in mr_data_rolling.Patient.unique():
        df = mr_data_rolling[mr_data_rolling['Patient'] == p]
        df = df.drop(columns=['Patient']) # removing un-aggregatable columns
        df = df.rolling(window=('%ss' % window_size),
                        # center = True,
                        # win_type='triang',
                        min_periods=int(minimum_pos),
                        on='Pos_as_date')
        df = df.mean()
        df['Patient'] = p # restoring un-aggregatable columns
        dfs.append(df)
    mr_data_rolling = pd.concat(dfs)

    mr_data_rolling['Pos'] = mr_data_rolling['Pos_as_date'].astype(int)/ 10**9
    mr_data_rolling['MR_mean_rolling'] = np.exp(mr_data_rolling['log_MR'])
    mr_data_rolling = mr_data_rolling.drop(columns=['Pos_as_date','log_MR','MR'])
    mr_data_rolling = mr_data_rolling[mr_data_rolling['MR_mean_rolling'].notnull()]
    mr_data_rolling = mr_data_rolling.astype({"Pos": int})
    mr_data_rolling = mr_data_rolling.sort_values(['Patient','Pos'])
    # mr_data_rolling = mr_data_rolling.astype({"Pos": int, "Patient": int})

    # threshold
    # applying here- identification of non-hh in suspected hh zones- reduces FP in hotspot inference!
    mr_data_rolling = mr_data_rolling[mr_data_rolling['MR_mean_rolling'] > MR_MEAN_OF_MEANS]

    # basic summary
    print('Hotspot count: %s' % len(mr_data_rolling))
    per_patient = mr_data_rolling.groupby('Patient').agg(['count'])['MR_mean_rolling'].sort_values(['count'], ascending=False)
    print(per_patient)
    print(per_patient.mean())

    # additional analytics
    # v1- add regions only
    # mr_data_rolling['region'] = mr_data_rolling.apply(lambda row: get_ET86_region(row['Pos']), axis=1)
    # print(mr_data_rolling.groupby(['region']).agg(['count']))
    # print(mr_data_rolling.groupby(['Patient','region']).agg(['count']))

    # v2 - unify with freqs and get all data
    # TODO- if used- need to notice line multipication
    unified_freqs = pd.read_csv('/Users/omer/PycharmProjects/SternLab/RG_HIVC_analysis/runs/{}/unified_freqs_filtered_verbose.csv'.format(run_name))
    # mr_data_rolling = mr_data_rolling.rename(columns={'Patient': 'ind_id'})
    # merge 1- multiplying rows- rank, sample_id
    # freqs_with_hotspots = unified_freqs.merge(mr_data_rolling, on= ['ind_id', 'Pos'], how= 'left', sort=False)
    #
    # # merge 2- take one row from unified freqs
    unified_freqs_reduced = unified_freqs.groupby(['ind_id', 'Pos']).first()
    hotspots_with_freqs = mr_data_rolling.merge(unified_freqs_reduced, on= ['ind_id', 'Pos'], how= 'left', sort=False)

    if generate_outputs:
        fname = 'syn_high_mr_hotspots_v2_' + run_name + '_' + minimum_pos + '_from_' + window_size
        mr_data_rolling.to_csv(mr_dir + fname + '.csv', index= False)

    g2 = sns.relplot(x="Pos",
                     y="MR_mean_rolling",
                     col='Patient',
                     col_order= plots_order,
                     # hue='Mutation',
                     color='firebrick',
                     col_wrap=7,
                     # kind='line',
                     data=mr_data_rolling)

    # plot adjustments
    g2.set(yscale="log")
    plt.ylim(1e-8, 1e-1)
    plt.xlim(0, 9026)
    g2.fig.suptitle('Hotspot window definition: {}/{}\nfound: {}'.format(minimum_pos, window_size, len(mr_data_rolling)), y=0.2)

    # extract plot
    if generate_outputs:
        plt.savefig(fname= mr_dir + fname + '.png')
    else:
        plt.show()


def initial_naive_analysis():
    # get unified
    unified_freq_df = pd.read_csv('/Users/omer/PycharmProjects/SternLab/RG_HIVC_analysis/ET86_2s/with_muts_coordinated/unified_freqs_with_ids_ysi.csv')

    # freqs sanity check
    # print(unified_freq_df[(unified_freq_df['Freq'] > 0.5) & (unified_freq_df['Mutation_type'] == 'stop')].head().to_string())
    # print(len(unified_freq_df[(unified_freq_df['Freq'] > 0.5) & (unified_freq_df['Mutation_type'] == 'stop')]))

    # additional filters
    # G->A only
    # print(len(unified_freq_df))
    # unified_freq_df = unified_freq_df[unified_freq_df['Mutation'] == 'GA']
    # print(len(unified_freq_df))

    # manual verification analysis
    # verify_chosen_stats(unified_freq_df)
    # potential_analysis(unified_freq_df)

    # plot_mutation_rate_distribution(unified_freq_df)
    # plot_mutation_rate_conf_interval(unified_freq_df)

    # regression line extraction
    # calc_regression_lines(unified_freq_df)


if __name__ == "__main__":
    # plot_fits_mr_genomewide()
    plot_mr_hotspots()
