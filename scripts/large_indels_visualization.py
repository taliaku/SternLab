import pandas as pd
import matplotlib.pyplot as plt
import glob


def plot_suspicious_reads(suspicious_reads_file, output_fig_path= None, sort=True):
    """
    visualize multi-mapped reads
    :param suspicious_reads_file: file generated from new python pipeline
    """

    print('plot_suspicious_reads: {}'.format(file))
    suspicious_reads = pd.read_csv(suspicious_reads_file, sep='\t')
    # suspicious_reads = suspicious_reads.head(100)

    # inverse minus mapping
    sp = suspicious_reads[suspicious_reads.plus_or_minus == 'plus']
    # TODO- inverse read_seq, ref_seq columns
    sm = suspicious_reads[suspicious_reads.plus_or_minus == 'minus'].rename(
        columns={'ref_start': 'ref_end', 'ref_end': 'ref_start'})
    df = pd.concat([sp, sm])
    # sort by position
    if sort:
        df = df.sort_values(by=['ref_start', 'ref_end']).reset_index()

    # remove duplicates
    df['interval'] = list(zip(df['ref_start'], df['ref_end']))
    df = df[['read_id', 'interval']].drop_duplicates()

    # add interval serial id (per read) for coloring
    reads = df['read_id'].unique().tolist()
    print(len(reads))
    df['interval_id'] = df.groupby('read_id').cumcount()

    # plotting
    print('plotting')
    plt.figure(figsize=(20, 10))
    df.apply(lambda row: plt.plot((row['interval'][0], row['interval'][1]),
                                  (reads.index(row['read_id']), reads.index(row['read_id'])),
                                  row['interval_id'], #TODO- color by id (dont know how it colors at the moment)
                                  marker='o'),
             axis=1)
    plt.xlabel('Position')
    plt.ylabel('read ID')

    if output_fig_path == None:
        output_fig_path = suspicious_reads_file.split(".")[0] + ".pdf"
    plt.savefig(output_fig_path)
    plt.show()

    # TODO- export table too (with seq columns)


if __name__ == '__main__':
    suspicious_reads_files = glob.glob('/sternadi/nobackup/volume1/HIVB_shafer_accungs_seq/new_pipeline_accungs_v2_15_iter/env*/suspicious_reads.tsv')
    print('len(suspicious_reads_files): {}'.format(len(suspicious_reads_files)))

    for file in suspicious_reads_files:
        sample = file.split('/')[-2]
        print(sample)
        plot_suspicious_reads(file, sort=True)
        plot_suspicious_reads(file,
                              output_fig_path= '/sternadi/home/volume1/shared/analysis/HIV_shafer/indels/potential_large_indels_{}.pdf'.format(sample))
