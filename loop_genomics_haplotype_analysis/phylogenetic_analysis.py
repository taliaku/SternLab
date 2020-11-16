from Bio import SeqIO,Seq
from Bio.Seq import MutableSeq

import pandas as pd

from loop_genomics_haplotype_analysis.summarize_visualize_strain_results import shafer_root_dir, accungs_root_dir


def strain_sequence_from_con_and_mutation_list(consensus_file, mutation_dict):
    with open(consensus_file, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            consensus_seq = record.seq
            tmp_record = record # TODO- must be a better way than this

    mutable_consensus = MutableSeq(str(consensus_seq))
    # print(len(consensus_seq))
    print('len(mutation_dict): {}'.format(len(mutation_dict)))

    for key, value in mutation_dict.items():
        # sanity
        if mutable_consensus[key-1] != value[0]:
            print('unexpected base in pos [{}]: [{}] instead of [{}]'.format(key, mutable_consensus[key-1], value[0]))
        else:
            mutable_consensus[key-1] = value[1]

    new_sequence = Seq.Seq(str(mutable_consensus))
    # print(len(new_sequence))
    tmp_record.seq = new_sequence # TODO- must be a better way than this

    return tmp_record


def generate_haplotype_sequences_v2():
    samples = ['env3', 'env4', 'env5', 'env6', 'env7', 'env8', 'env9', 'env10', 'env11', 'env12', 'env13', 'env14', 'env15']
    # samples = ['env10']
    mutation_data = pd.read_csv('{}/results_merged/strain_data_v2.csv'.format(shafer_root_dir))

    for sample in samples:
        print(sample)
        # accungs
        accungs_con_file = '{}/envs_output/{}/last_ref.fasta'.format(accungs_root_dir, sample)
        source = 'accungs'
        generate_fasta_strains(accungs_con_file, mutation_data, sample, source)

        # loop
        loop_con_file = '{}/loop_genomics_pipeline/envs_output/{}/last_ref.fasta'.format(shafer_root_dir, sample)
        source = 'loop'
        generate_fasta_strains(loop_con_file, mutation_data, sample, source)


def generate_fasta_strains(consensus_filename, mutation_data, sample, source):
    sample_source_muts = mutation_data[(mutation_data['source'] == source) & (mutation_data['sample'] == sample)]
    for strain in sample_source_muts['strain'].unique():
        print(strain)
        # get mutation list
        df = sample_source_muts[sample_source_muts['strain'] == strain]
        mutation_dict = dict(zip(df['Pos'], df['mutation']))

        strain_record = strain_sequence_from_con_and_mutation_list(consensus_filename, mutation_dict)
        # save fasta
        strain_name = '{}_{}_strain{}'.format(sample, source, strain)
        # TODO- improve name- '{}_strain{}_freqX_mutationsY.fasta'

        strain_record.id = strain_name
        strain_record.description = ''

        output_path = '{}/strain_phylogeny/{}.fasta'.format(shafer_root_dir, strain_name)
        with open(output_path, "w") as output_handle:
            SeqIO.write(strain_record, output_handle, "fasta")


if __name__ == '__main__':
    generate_haplotype_sequences_v2()

    # Generate tree & games with it-
    #   MAFFT online + FigTree


