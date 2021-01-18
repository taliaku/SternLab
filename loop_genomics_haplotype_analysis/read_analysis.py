from Bio import SeqIO, Seq, SeqRecord
from Bio.SeqIO import FastaIO
import pandas as pd
import glob

def fastq_filter_by_qual(fastq_filepath, output_path, trim = False):
    with open(output_path, "w") as out_handle:
        # fasta_out = FastaIO.FastaWriter(out_handle, wrap=None)
        for record in SeqIO.parse(fastq_filepath, "fastq"):
            merged = zip(record.seq, record.letter_annotations["phred_quality"])
            if trim:
                new_seq = ''.join([base[0] if base[1]>10 else '' for base in merged ])
            else:
                new_seq = ''.join([base[0] if base[1]>10 else 'N' for base in merged ])

            rec = SeqRecord.SeqRecord(id=record.id, description='', seq=Seq.Seq(new_seq))

            SeqIO.write(rec, out_handle, "fasta")
            # fasta_out.write_record(rec)


def fastq_filter_by_qual_main():
    # input_file_template = '/sternadi/datasets/volume2/HIV_shafer_loop_genomics/fastq_only/s{}/s{}.fastq'
    input_file_template = '/sternadi/datasets/volume2/HIV_shafer_loop_genomics/fastq_long_only/s{}/s{}_long.fastq'
    # output_file_template = '/sternadi/datasets/volume2/HIV_shafer_loop_genomics/fasta_trimmed/env{}.trimmed.fasta'
    output_file_template = '/sternadi/datasets/volume2/HIV_shafer_loop_genomics/fasta_long_trimmed/env{}.trimmed.fasta'
    for i in range(3, 16):
        print('env{}'.format(i))
        print(input_file_template.format(i, i))
        print(output_file_template.format(i))

        fastq_filter_by_qual(input_file_template.format(i, i),
                             output_file_template.format(i),
                             trim=True)


def generate_blast_read_id_list():
    # loop_pipeline_root_dir = '/sternadi/home/volume1/shared/analysis/HIV_shafer/new_pipeline_loop_v2_q10'
    loop_pipeline_root_dir = '/Volumes/STERNADILABHOME$/volume1/shared/analysis/HIV_shafer/new_pipeline_loop_v2_q10'
    # files = glob.glob('{}/env*/blast.tsv'.format(loop_pipeline_root_dir))
    files = glob.glob('{}/env*/blast.tsv'.format(loop_pipeline_root_dir))

    for file in files:
        print(file)
        blast = pd.read_csv(file, sep='\t', header=None)

        # filter by length
        read_id_blast_mapped_long_plus = blast[(blast[7] > 1000) & (blast[6] == 'plus')]
        read_id_blast_mapped_long_minus = blast[(blast[7] > 1000) & (blast[6] == 'minus')]
        # print(read_id_blast_mapped_long_plus.head())

        # get read_ids & dedup
        read_id_blast_mapped_long_plus = read_id_blast_mapped_long_plus[0].drop_duplicates()
        read_id_blast_mapped_long_minus = read_id_blast_mapped_long_minus[0].drop_duplicates()

        # remove intersect of minus/plus
        # TODO- verify this code
        read_id_blast_mapped_long_minus = read_id_blast_mapped_long_minus.drop(read_id_blast_mapped_long_minus.intersection(read_id_blast_mapped_long_plus))

        # export
        read_id_blast_mapped_long_plus.to_csv(file + '.good_mapps_1000.dedup.plus', index=False, header=False, sep='\t')
        read_id_blast_mapped_long_minus.to_csv(file + '.good_mapps_1000.dedup.minus', index=False, header=False, sep='\t')

        # 3seq format
        # pd.DataFrame(read_id_blast_mapped_long_plus).T.to_csv(file + '.good_mapps_1800.dedup.row', index=False, header=False, sep=' ')

        # debug 'bad' reads (short etc.)
        # blast_bad = blast[blast[7]/blast [3] > 1000]
        # blast = pd.read_csv(file, header=None, sep='\t')
        # pd.DataFrame(blast).T.to_csv(file + '.read_id.dedup2', index=False, header=False, sep=' ')



if __name__ == '__main__':
    # fastq_filter_by_qual_main()
    generate_blast_read_id_list()