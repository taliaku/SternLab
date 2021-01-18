import sys
import warnings
import numpy as np

def pi_diversity_calc(data, pivot_cols=[]):
    """
    Calculates PI diversity per position, than calculates mean per group according to pivot_vols.
    Assumes data is not indexed.

    :param data: freqs file + columns additional for pivoting
    :param pivot_cols:
    :return:
    """

    # TODO- call apply_pi_related_filters() here?

    # keep only pivot & data columns
    selecting_cols = pivot_cols[:]
    selecting_cols.extend(["ref_pos", "read_base", "base_count", "base_rank"]) # base_count = coverage * frequncy
    selecting_cols = list(set(selecting_cols))
    data = data[selecting_cols]

    # Mark major\minor variants per position
    data["base_rank"] = np.where(data["base_rank"] == 0, "Major", "Minor")

    group_cols = pivot_cols[:]
    group_cols.extend(["ref_pos", "base_rank"])
    group_cols = list(set(group_cols))

    # Getting minor base frequency
    # (selecting max of minor bases)
    data = data.groupby(group_cols)['base_count'].aggregate(max).unstack('base_rank').reset_index()
    # TODO: use all minor variants? (not only max?)

    if 'Minor' not in data.columns:
        # raise Exception("No minor frequencies found. check filtering of input data")
        warnings.warn("No minor frequencies found. check filtering of input data")
        return None

    data["Minor"] = data["Minor"].fillna(0)
    data["Total"] = data["Major"] + data["Minor"]

    # Pi diversity calculation by positions
    # (more on this in https://academic.oup.com/ve/article/5/1/vey041/5304643)
    def pairwise_differences_proportion(row):
        if row["Minor"] == 0:
            return 0
        else:
            total_part = (row["Total"] * (row["Total"] - 1))
            numerator = total_part - ((row["Major"] * (row["Major"] - 1)) + (row["Minor"] * (row["Minor"] - 1)))
            denominator = total_part
            return numerator * 1.0 / denominator

    data["pdp"] = data.apply(lambda row: pairwise_differences_proportion(row), axis=1)

    if len(pivot_cols) != 0:
        bysample_diversity = data.groupby(pivot_cols)['pdp'].agg(['count', 'sum']).reset_index()
        bysample_diversity["Pi"] = bysample_diversity["sum"] * 1.0 / bysample_diversity["count"]
        output_cols = pivot_cols[:]
        output_cols.append("Pi")
        pis = bysample_diversity[output_cols]
    else:
        pis = data['pdp'].mean()

    return pis


def apply_pi_related_filters(freqs_file,
                             coverage_threshold=0,
                             frequency_threshold=0,
                             base_count_threshold=10,
                             probability_threshold=0.8,
                             transitions_only= True,
                             remove_insertions= True,
                             positions_to_include = (0, sys.maxsize)  # start_pos=0, end_pos= sys.maxsize
                             ):

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
        filtered_freqs_file["mutation_type"] = filtered_freqs_file['read_base'].astype(str) + filtered_freqs_file['ref_base'].astype(str)
        filtered_freqs_file = filtered_freqs_file[
            filtered_freqs_file["mutation_type"].isin(['GA', 'AG', 'GG', 'AA', 'CT', 'TC', 'CC', 'TT'])]

    # remove insertions
    if remove_insertions:
        # filtered_freqs_file = filtered_freqs_file[filtered_freqs_file['ref_base'] != "-"]
        filtered_freqs_file = filtered_freqs_file[filtered_freqs_file['ref_pos'] % 1 == 0]

    # keep specified positions
    filtered_freqs_file = filtered_freqs_file[(filtered_freqs_file["ref_pos"] >= positions_to_include[0]) &
                                        (filtered_freqs_file["ref_pos"] <= positions_to_include[1])]

    if filtered_freqs_file.empty:
        warnings.warn("No relevant data after filtering")
        return None
    else:
        print("after filtering: {} rows ({}% from original)".format(len(filtered_freqs_file), len(filtered_freqs_file)*100/len(freqs_file)))
        return filtered_freqs_file