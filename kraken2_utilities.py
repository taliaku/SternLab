#! /usr/local/python_anaconda/bin/python3.4

import pandas as pd
import numpy as np
import glob
import os

def load_data(dir_path):
    """
    process kraken output, load all report data in dir_path and combine them to a single dataframe
    :param dir_path: a directory containing the report files
    :return: data frame with all sample names
    """
    dfs = []
    # column names as described in here:
    #https: // ccb.jhu.edu / software / kraken2 / index.shtml?t = manual
    names = ["percentage covered", "number covered", "number assigned", "rank code", "NCBI taxonomic ID", "name"]
    files = glob.glob(os.path.join(dir_path, "*report.tab"))
    for f in files:
        df = pd.read_table(f, names=names)
        df['sample'] = os.path.basename(f).split("_report")[0]
        dfs.append(df)
    df = pd.concat(dfs)
    df['name'] = df['name'].apply(lambda x: x.strip())
    return df

def map_index_to_kingdom(df):
    """
    map row index and taxonomy for later use in map taxonomy
    """
    d = {}
    df = df.reset_index()
    df['index'] = df.index
    for sample in df['sample'].unique():
        d[sample] = {}
        curr = df[df['sample'] == sample]
        curr = curr[curr['rank code'] == 'D']
        names = curr['name'].unique()
        for name in names:
            d[sample][str(curr[curr['name'] == name]['index'].values[0])] = name
    return d

def map_taxonomy(mapper, sample, row_loc):
    """
    add taxonomic kingdom hierarchy to each row
    """
    indices_to_tax = mapper[sample]
    indices = sorted([int(k) for k in indices_to_tax.keys()])
    # search the location of the idx and match to the relevant tax (the previous index in the list)
    idx = np.searchsorted(indices, row_loc)
    tax = indices_to_tax[str(indices[idx - 1])]
    return tax

def filter_df(df, thresh=1, rank='G', include_unclussified=True):
    """
    filter dataframe to contain only a specific taxonomy level
    :param df: a dataframe with all taxonomies
    :param thresh: threshold for read coverage percentage filter. default = 0.01
    :param rank: the taxonomic rank to filter by - default = G
    """
    mapper = map_index_to_kingdom(df)
    df = df.reset_index()
    df['index'] = df.index
    if include_unclussified:
        df = df[((df['rank code'] == rank) & (df['percentage covered'] >= thresh) | (df['rank code'] == "U"))]

    else:
        df = df[(df['rank code'] == rank) & (df['percentage covered'] >= thresh)]
    df['kingdom'] = df.apply(lambda row: map_taxonomy(mapper, row['sample'], int(row['index'])), axis=1)
    if include_unclussified:
        df.loc[df["name"] == "unclassified", "kingdom"] = "unclassified"
    return df


def merge_metadata(df, metadata_path = r"../data/all sequening run metadata.csv", gisaid=True):
    """
    merge kraken output dataframe per sample with metadata content
    """
    metadata = pd.read_csv(metadata_path)
    if gisaid:
        metadata['sample_num'] = metadata['sample']
        metadata = metadata[metadata['country'] == "Israel"]
        metadata = metadata[['age', 'sex', 'date','division', 'country', 'originating_lab', "sample_num", "averageCT"]]
        df['sample_num'] = df['sample'].apply(lambda x: x.split("_")[0])
        merged = df.merge(metadata, on="sample_num")
    else:   #SRA data
        try:
            metadata = metadata[["Run", "Instrument", "LibraryLayout", "LibrarySelection",	"LibrarySource", "Bases",
                                "Library Name", "Collection_Date", "isolate", "BioSampleModel", "geo_loc_name_country",
                                 "geo_loc_name_country_continent", "geo_loc_name", "Host", "collected_by", "host_disease"]]
        except:
            metadata = metadata
        df['Run'] = df['sample']
        merged = df.merge(metadata, on="Run")

    return merged


def main():
    data_files_path = r"/Users/daniellemiller/Google Drive/covid19/microbiome/kraken_outputs/merged"#
    output_dir = r"/Users/daniellemiller/Google Drive/covid19/microbiome/analysis"
    alias = "merged_0.1pct"
    data = kraken_utils.load_data(data_files_path)
    filtered_data = kraken_utils.filter_df(data, thresh=0.1)