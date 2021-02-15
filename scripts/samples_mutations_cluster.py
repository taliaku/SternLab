#! /usr/local/python_anaconda/bin/python3.4
import pandas as pd
import matplotlib
matplotlib.use('agg')
import seaborn as sns
import matplotlib.pyplot as plt
import os,sys,inspect
parentdir = os.path.dirname(os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe()))))
sys.path.insert(0,parentdir)
from optparse import OptionParser
from file_utilities import check_filename

def samples_mutations_cluster(df, output_excel, threshold_freq=0.1, threshold_read_count=10, sample_column='File'):
    '''
    This function gets a dataframe of multiple frequency files, and creates a clustered pivot
    table of mutations that appear at significant frequencies.
    Important! - clustering is run using only mutations whose frequency is known in all samples,
    the rest of the mutations (containing NAs in some columns) are appended at the bottom of the
    excel. The clustermap plot is also saved.
    @df - dataframe object of freqs
    @output_excel - file path
    @threshold_freq - a mutation needs to cross this threshold in at least one sample to be included.
    @threshold_read_count - rows that do not cross this are removed from analysis.
    @sample_column - column name that separates into samples (file, sample etc.)
.
    (Recommendation: use conditional formatting in excel afterwards to color according to frequency.)
    '''
    df['full_mutation'] = df.ref_base + df.ref_position.astype(int).astype(str) + df.base
    df = df[df.coverage > threshold_read_count]
    mutations_to_keep = df[(df.ref_base != df.base) & (df.ref_base != '-') & (df.frequency > threshold_freq)].full_mutation.drop_duplicates().tolist()
    mutations_to_keep = df[df.full_mutation.isin(mutations_to_keep)]
    to_pivot = mutations_to_keep.pivot_table(values='frequency', index=['full_mutation', 'ref_position'], columns=sample_column)
    to_pivot_na = to_pivot[to_pivot.isnull().any(axis=1)]
    to_pivot = to_pivot.dropna()
    clustergrid = sns.clustermap(to_pivot, figsize=(15,15), xticklabels=True, yticklabels=True, method='weighted')
    to_pivot['mutation_order'] = pd.Categorical(to_pivot.index, [to_pivot.index[s] for s in clustergrid.dendrogram_row.reordered_ind])
    to_pivot = pd.concat([to_pivot.sort_values('mutation_order'), to_pivot_na], sort=False)
    to_pivot = to_pivot[[to_pivot.columns[s] for s in clustergrid.dendrogram_col.reordered_ind]]
    to_pivot.to_excel(output_excel)
    clustergrid.savefig(output_excel.replace('.xlsx', '.png'), dpi=800)


def main():
    parser = OptionParser("usage: %prog [options]")
    parser.add_option("-i", "--input", dest="input", help="input dataframe csv path")
    parser.add_option("-o", "--output", dest="output", help="output excel path")
    parser.add_option("-f", "--threshold_freq", dest="threshold_freq", help="a mutation needs to cross this threshold in at least one sample to be included. default 0.1", default=0.1)
    parser.add_option("-r", "--threshold_read_count", dest="threshold_read_count", help="rows that do not cross this are removed from analysis. default 10", default=10)
    parser.add_option("-s", "--sample_column", dest="sample_column", help="column name that separates into samples (file, sample etc.). default 'File'", default='File')

    (options, args) = parser.parse_args()
    input = options.input
    output = options.output
    input = check_filename(input)
    output = check_filename(output, Truefile=False)

    samples_mutations_cluster(pd.read_csv(input, low_memory=False), output, options.threshold_freq, options.threshold_read_count, options.sample_column)

if __name__ == '__main__':
    main()
