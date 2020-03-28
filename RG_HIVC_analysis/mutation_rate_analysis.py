import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns;
from scipy import stats

from RG_HIVC_analysis.data_adaptations import generate_unified_filtered_verbose_freqs_df

sns.set_context("poster")


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


if __name__ == "__main__":
    # generate unified
    # unified_freq_df = generate_unified_filtered_verbose_freqs_df()

    # get unified
    unified_freq_df = pd.read_csv('/Users/omer/PycharmProjects/SternLab/RG_HIVC_analysis/ET86_4s/with_muts_and_con_as_ref/unified_freqs_filtered_verbose.csv')
    # unified_freq_df = pd.read_csv('/Users/omer/PycharmProjects/SternLab/RG_HIVC_analysis/ET86_4s/with_muts_and_con_as_ref/TASPX119494_S74.freqs', sep='\t')
    # unified_freq_df = pd.read_csv('/Users/omer/PycharmProjects/SternLab/RG_HIVC_analysis/ET86_2s/with_muts_coordinated/unified_freqs_with_ids_ysi.csv')

    # freqs sanity check
    # print(unified_freq_df[(unified_freq_df['Freq'] > 0.5) & (unified_freq_df['Mutation_type'] == 'stop')].head().to_string())
    # print(len(unified_freq_df[(unified_freq_df['Freq'] > 0.5) & (unified_freq_df['Mutation_type'] == 'stop')]))

    # additional filters
    # G->A only
    # print(len(unified_freq_df))
    # unified_freq_df = unified_freq_df[unified_freq_df['Mutation'] == 'GA']
    # print(len(unified_freq_df))

    # plot_mutation_rate_distribution(unified_freq_df)
    # plot_mutation_rate_conf_interval(unified_freq_df)

    # manual verification analysis
    # verify_chosen_stats(unified_freq_df)
    # potential_analysis(unified_freq_df)

    # regression line extraction
    # calc_regression_lines(unified_freq_df)
