import os
import logging
import pandas as pd


def _get_last_end_value_location(string, substring, substring_end_loc, end_values):
    end_value_locations = []
    for end_value in end_values:
        end_value_loc = string.find(end_value, substring_end_loc)
        if end_value_loc != -1:
            end_value_locations.append(end_value_loc)
    if len(end_value_locations)==0:
        ret = None
    else:
        ret = min(end_value_locations)
    return ret


def find_all(string, substring):
    ret = []
    start_pos = string.find(substring)
    while start_pos != -1:
        ret.append(start_pos)
        start_pos = string.find(substring, start_pos+1)
    return ret


def get_value_after_substring(string, substring, start=None, end=None):
    end_values = (',', '\n', ':', ' ')
    substring_start_loc = string.find(substring, start, end)
    if substring_start_loc != -1:
        substring_end_loc = substring_start_loc + len(substring)
        return_value = string[substring_end_loc : _get_last_end_value_location(string, substring, substring_end_loc, end_values)]
    else: 
        return_value = None
    return return_value


def _get_repeats_contribution_data(summary_data):
    ret = {}
    for repeat_location in find_all(summary_data, ' repeats contributing '):
        repeat_number_location_start = summary_data.rfind('\n', repeat_location -5, repeat_location) + 1
        number_of_repeat = summary_data[repeat_number_location_start : repeat_location]
        ret[f'Reads contributed from {number_of_repeat} repeats'] = get_value_after_substring(summary_data, ' repeats contributing ', start=repeat_location)
        ret[f'Bases contributed from {number_of_repeat} repeats'] = get_value_after_substring(summary_data, ' reads and ', start=repeat_location)        
    return ret


def aggregate_summaries(project_dir_path):
    """
    Returns: a dataframe aggregating all the summaries in a project.
    """
    sub_dir_list = [name for name in os.listdir(project_dir_path) if os.path.isdir(os.path.join(project_dir_path,name))]
    data_categories_with_spaces = ['mode', 'task', '%id for blast', 'e-value', 'number of repeats used', 'q-score', 'protocol']
    data_categories_with_colon = ['Total number of reads', 'Number of reads mapped to reference', "% of mapped reads", 
                                  "Total number of reads contributing to frequency count", "% of reads contributing to frequency count",
                                  "Total number of bases contributing to frequency count", "Number of positions in reference genome",
                                  "Sum of Mutations"]
    #TODO: add the missing lines!
    data = {}
    for sub_dir_name in sub_dir_list:
        sub_dir_data = {}
        summary_path = os.path.join(project_dir_path, sub_dir_name, 'Summary.txt')
        with open(summary_path) as f:
            summary_data = f.read()
        # append variable columns
        value_after_x = get_value_after_substring(summary_data, "positions with min coverage x")
        if value_after_x:
            data_categories_with_colon.append(f"Number of positions with min coverage x{value_after_x}")
            data_categories_with_colon.append(f"% of positions with min coverage x{value_after_x}")
        # reads contributed from x repeats
        
        # get data
        for column in data_categories_with_colon:
            sub_dir_data[column] = get_value_after_substring(summary_data, f"{column}: ")
        for column in data_categories_with_spaces:
            sub_dir_data[column] = get_value_after_substring(summary_data, f"{column} = ")
        for column, value in _get_repeats_contribution_data(summary_data).items():
            sub_dir_data[column] = value
        data[sub_dir_name] = sub_dir_data
    return pd.DataFrame(data=data).T
