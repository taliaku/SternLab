import sys
import warnings

import numpy as np


def pi_diversity_calc(data, pivot_cols=[], interval = (0, sys.maxsize)): # start_pos=0, end_pos= sys.maxsize
    """
    Calculates PI diversity per position, than calculates mean per group according to pivot_vols. Assumes data is not indexed.

    # Example for using pivot:
    # pi_rates_by_sample = pi_diversity_calc(data=freqs, pivot_cols=['sample'])
    # pi_rates_by_sample = pi_rates_by_sample.groupby('sample').mean().reset_index()

    :param data: freqs file, with optional additional columns for pivoting (grouping)
    :param pivot_cols:
    :param min_read_count:
    :return:
    """
    def pairwise_differences_proportion(row):
        if row["Minor"] == 0:
            return 0
        total_part = (row["Total"] * (row["Total"] - 1))
        numerator = total_part - ((row["Major"] * (row["Major"] - 1)) + (row["Minor"] * (row["Minor"] - 1)))
        denominator = total_part
        return numerator * 1.0 / denominator

    # Filters
    # TODO- extract this filter too\ insert all others here
    # choose interval

    filtered_data = data[(data["ref_position"] >= interval[0]) & (data["ref_position"] <= interval[1])]
    filtered_data = data
    if filtered_data.empty:
        warnings.warn("No relevant data after filtering. Skipping")
        return None

    filtered_data['counts_for_position'] = np.round(filtered_data['coverage'] * filtered_data['frequency'])
    selecting_cols = pivot_cols[:]
    selecting_cols.extend(["ref_position", "base", "counts_for_position", "rank"])
    # selecting_cols = list(set(selecting_cols))
    filtered_data = filtered_data[selecting_cols]

    filtered_data["rank"] = np.where(filtered_data["rank"] == 0, "Major", "Minor")

    group_cols = pivot_cols[:]
    group_cols.extend(["ref_position", "rank"])
    # group_cols = list(set(group_cols))

    # selecting max on cfp, per Pos\Rank (Major\minor?)- than will be graded per Pos, and summed
    # TODO: use all minor variants (not only max)
    filtered_data = filtered_data.groupby(group_cols)['counts_for_position'].aggregate(max).unstack().reset_index()
    if 'Minor' not in filtered_data.columns:
        return 0

    filtered_data["Total"] = filtered_data["Major"] + filtered_data["Minor"]

    filtered_data["pdp"] = filtered_data.apply(lambda row: pairwise_differences_proportion(row), axis=1)

    pivot_cols.extend(["ref_position"])
    if any(pivot_cols):
        bysample_diversity = filtered_data.groupby(pivot_cols)['pdp'].agg(['count', 'sum']).reset_index()
        bysample_diversity["Pi"] = bysample_diversity["sum"] * 1.0 / bysample_diversity["count"]
        output_cols = pivot_cols[:]
        output_cols.append("Pi")
        pis = bysample_diversity[output_cols]
    else:
        pis = filtered_data['pdp'].mean()

    return pis


def apply_pi_related_filters(freqs_file,
                             coverage_threshold=0,
                             frequency_threshold=0,
                             base_count_threshold=0,
                             probability_threshold=0,
                             transitions_only= False,
                             remove_indels= False):

    filtered_freqs_file = freqs_file.copy()

    # base_count threshold
    filtered_freqs_file = filtered_freqs_file[filtered_freqs_file["base_count"] > base_count_threshold]

    # coverage threshold
    filtered_freqs_file = filtered_freqs_file[filtered_freqs_file["coverage"] > coverage_threshold]

    # set low frequencies to 0
    filtered_freqs_file["frequency"] = np.where(filtered_freqs_file["frequency"] >= frequency_threshold, filtered_freqs_file["frequency"], 0)

    # probability threshold
    filtered_freqs_file = filtered_freqs_file[filtered_freqs_file["probability"] > probability_threshold]

    # transitions only
    if transitions_only:
        filtered_freqs_file["mutation_type"] = filtered_freqs_file['read_base'] + filtered_freqs_file['ref_position']
        filtered_freqs_file = filtered_freqs_file[
            filtered_freqs_file["mutation_type"].isin(['GA', 'AG', 'GG', 'AA', 'CT', 'TC', 'CC', 'TT'])]

    # remove indels
    if remove_indels:
        filtered_freqs_file = filtered_freqs_file[(filtered_freqs_file["read_base"] != "-") & (filtered_freqs_file["ref_position"] != "-")]

    return filtered_freqs_file