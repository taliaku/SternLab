

def _get_last_end_value_location(string, substring, substring_end_loc, end_values):
    end_value_locations = []
    for end_value in end_values:
        end_value_loc = string.find(end_value, start=substring_end_loc)
        if end_value_loc != -1:
            end_value_locations.append(end_value_loc)
    if len(end_value_locations)==0:
        ret = None
    else:
        ret = min(end_value_locations)
    return ret
            

def get_value_after_substring(string, substring, start=None, end=None):
    end_values = (',', '\n'. ':')
    substring_start_loc = string.find(substring, start, end)
    substring_end_loc = string_location + len(substring)
    return string[substring_end_loc : _get_last_end_value_location(string, substring, substring_end_loc, end_values)]


def aggregate_summaries(project_dir_path):
    """
    Returns: a dataframe aggregating all the summaries in a project.
    """
    sub_dir_list = [name for name in os.listdir(project_dir_path) if os.path.isdir(name)]
    data_categories_with_spaces = ['mode', 'task', 'for blast', 'e-value', 'number of repeats used', 'q-score', 'protocol']
    data_categories_with_colon = ['Total number of reads', 'Number of reads mapped to reference', "% of mapped reads", 
                                  "Total number of reads contributing to frequency count", "% of reads contributing to frequency count",
                                  "Total number of bases contributing to frequency count", "Number of positions in reference genome",
                                  "Sum of Mutations"]
    #TODO: add the missing lines!
    data = {}
    for sub_dir_name in sub_dir_list:
        sub_dir_data = {}
        summary_path = os.path.join(dir_path, sub_dir_name, 'Summary.txt')
        with open(summary_path) as f:
            summary_data = f.read()
        # append variable columns
        value_after_x = get_value_after_substring(summary_data, "positions with min coverage x")
        data_categories_with_colon.append(f"Number of positions with min coverage x"{value_after_x})
        data_categories_with_colon.append(f"% of positions with min coverage x"{value_after_x})
        # get data
        for column in data_categories_with_colon:
            sub_dir_data[column] = get_value_after_substring(summary_data, f"{column}: ")
        for column in data_categories_with_spaces:
            sub_dir_data[column] = get_value_after_substring(summary_data, f"{column} = ")
        data[sub_dir_name] = sub_dir_data
    return data