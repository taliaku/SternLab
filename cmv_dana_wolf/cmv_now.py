# -*- coding: utf-8 -*-
"""
Created on Sat Dec  1 21:17:36 2018

@author: Noam
"""
import re
import pandas as pd
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
import numpy as np

def big_freqs_file():
    dfs = []
    freqs_list = ['Z:/volume1/noam/cmv/1_pipeline_output/1.freqs', 
                  'Z:/volume1/noam/cmv/2_pipeline_output/2.freqs', 
                  'Z:/volume1/noam/cmv/3_pipeline_output/3.freqs',
                  'Z:/volume1/noam/cmv/4_pipeline_output/4.freqs',
                  'Z:/volume1/noam/cmv/5_pipeline_output/5.freqs',
                  'Z:/volume1/noam/cmv/6_pipeline_output/6.freqs',
                  'Z:/volume1/noam/cmv/7_pipeline_output/7.freqs',
                  'Z:/volume1/noam/cmv/8_pipeline_output/8.freqs']
    for f in freqs_list:
        df = pd.read_csv(f, sep='\t')
        df['file'] = f.split('/')[-1].split('.')[0]
        dfs.append(df)
    df = pd.concat(dfs)
    return df



def compare_freqs(df, mutation_type_df, file1, file2, output_dir):
    # snps and deletions
    dfx = df[df.file == file1]
    dfy = df[df.file == file2]
    new_df = pd.merge(dfx, dfy, on=['Pos', 'Base', 'Ref'], suffixes=['_' + str(file1), '_' + str(file2)])
    deletions = new_df[(((new_df['Freq_' + str(file1)] - new_df['Freq_' + str(file2)]).abs() > 0.25) & (new_df.Base != new_df.Ref) & (new_df.Ref != '-') & (new_df.Base == '-'))]
    mutations = new_df[(((new_df['Freq_' + str(file1)] - new_df['Freq_' + str(file2)]).abs() > 0.25) & (new_df.Base != new_df.Ref) & (new_df.Ref != '-') & (new_df.Base != '-'))]
    mutations = pd.merge(mutations, mutation_type_df, on=['Pos', 'Ref', 'Base'], how='left')
    mutations = mutations[['Ref', 'Pos', 'Base', 'Freq_' + str(file1), 'Freq_' + str(file2), 'Read_count_' + str(file1), 'Read_count_' + str(file2), 'file_' + str(file1), 'file_' + str(file2), 'protein', 'mutation_type']]
    mutations.rename(columns={'mutation_type':'protein_mutation'}, inplace=True)
    conditions = [(mutations.protein_mutation == '[]'), (mutations.protein_mutation.isna()), ((mutations.protein_mutation != '[]') & ~(mutations.protein_mutation.isna()))]
    choices = ['synonymous', 'non_coding', 'non_synonymous']
    mutations['mutation_type'] = np.select(conditions, choices)
    
    deletions = pd.merge(deletions, mutation_type_df[['Pos', 'protein']].drop_duplicates(), on='Pos', how='left')
    deletions = deletions[['Ref', 'Pos', 'Base', 'Freq_' + str(file1), 'Freq_' + str(file2), 'Read_count_' + str(file1), 'Read_count_' + str(file2), 'file_' + str(file1), 'file_' + str(file2), 'protein']]
    # insertions
    #insertions = new_df[((new_df.Ref == '-') & ((new_df.Read_count_x - new_df.Read_count_y).abs() > 20))]    
    # positions appearing in one and not the other
    #weird_pos = pd.concat([dfx[(~(dfx.Pos.isin(dfy.Pos.tolist())))], dfy[(~(dfy.Pos.isin(dfx.Pos.tolist())))]])
    mutations.to_csv(output_dir + 'mutations_' + str(file1) + '_' + str(file2) + '.csv', index=False)
    deletions.to_csv(output_dir + 'deletions_' + str(file1) + '_' + str(file2) + '.csv', index=False)
    #insertions.to_csv(output_dir + 'insertions_' + file1 + '_' + file2 + '.csv', index=False)
    #weird_pos.to_csv(output_dir + 'weird_positions_' + file1 + '_' + file2 + '.csv', index=False)
    return
 
#for i in range(1,9):
#    for j in range(i+1,9):
#        if i != j:
#            compare_freqs(freqs, mutation_type, i, j, 'Z:/volume1/noam/cmv/comparisons2/')

def estimate_insertion_freq(df):
    read_counts = df[(df.Ref != '-')][['file', 'Pos', 'Read_count']].drop_duplicates()
    read_counts.rename(columns={'Read_count':'estimated_read_count', 'Pos':'rounded_pos'}, inplace=True)
    insertions = df[(df.Ref == '-')]
    not_insertions = df[(df.Ref != '-')]
    insertions['rounded_pos'] = insertions.Pos.astype(int).astype(float)
    insertions = pd.merge(insertions, read_counts, how='left', on=['file', 'rounded_pos'])
    insertions['estimated_freq'] = insertions.Read_count / insertions.estimated_read_count
    df = pd.concat([insertions, not_insertions])
    return df.sort_values(['file', 'Pos'])
    

#### mutation type #############3


def create_gene_dataframe():
    genes = []
    with open('Z:/volume1/noam/cmv/ncbi_genes.txt') as f:
        fi = f.read()
    for row in fi.splitlines():
        if not row.startswith('\t') and not row.startswith('>') and row != '':
            print('row' + row)
            protein = fi.split(row)[1].split('product\t')[1].split('\n')[0]
            start = int(row.split('\t')[0].replace('<', ''))
            end = int(row.split('\t')[1].replace('<', ''))
            complement = start > end
            genes.append((protein, start, end, complement))
        fi = fi.replace(row, '', 1)
    gene_df = pd.DataFrame(genes, columns=['protein', 'start', 'end', 'complement'])
    
    with open('Z:/volume1/noam/cmv/cmv_dana_ref.fasta') as f:
        reference_str = f.read()
        if reference_str.split('\n')[0].startswith('>'):
            reference_str = ''.join(reference_str.split('\n')[1:])
        else:
            reference_str = reference_str.replace('\n', '')
    gene_df['sequence'] = gene_df.apply(lambda x: get_rna_sequence(x.protein, x.start, x.end, x.complement, reference_str), axis=1)
    # not complement
    not_complement_df = gene_df[gene_df.complement == False].copy()
    not_complement_df['full_sequence'] = not_complement_df.sort_values('start', ascending=True).groupby('protein')['sequence'].transform(lambda x: x.sum())
    # is complement
    complement_df = gene_df[gene_df.complement == True].copy()
    complement_df['full_sequence'] = complement_df.sort_values('start', ascending=False).groupby('protein')['sequence'].transform(lambda x: x.sum())
    
    gene_df = pd.concat([complement_df, not_complement_df])
    gene_df['translated'] = gene_df.full_sequence.apply(translate_seq)
    gene_df.to_csv('Z:/volume1/noam/cmv/ncbi_genes_dataframe.csv', index=False)
    return gene_df

def get_rna_sequence(protein, start, end, complement, reference_str):
    if complement == False:
        return reference_str[start-1:end]
    else:
        return str(Seq(reference_str[end-1:start], IUPAC.unambiguous_dna).reverse_complement())
        
def translate_seq(seq):
    return str(Seq(seq, IUPAC.unambiguous_dna).translate())


def create_mutation_type_df():
    #gene_df = pd.read_csv('Z:/volume1/noam/cmv/ncbi_genes_dataframe.csv')
    #with open('Z:/volume1/noam/cmv/cmv_dana_ref.fasta') as f:
    gene_df = pd.read_csv('/sternadi/nobackup/volume1/noam/cmv/ncbi_genes_dataframe.csv')
    with open('/sternadi/nobackup/volume1/noam/cmv/cmv_dana_ref.fasta') as f:
        reference_str = f.read()
        if reference_str.split('\n')[0].startswith('>'):
            reference_str = ''.join(reference_str.split('\n')[1:])
        else:
            reference_str = reference_str.replace('\n', '')
    positions_dfs = []
    for i in range(len(reference_str)):
        positions_ref_bases = []
        relevant_genes_df = gene_df[(((gene_df.complement == False) & (i + 1 >= gene_df.start) & ( i + 1 <= gene_df.end)) | ((gene_df.complement == True) & (i + 1 >= gene_df.end) & ( i + 1 <= gene_df.start)))].copy()
        for base in ['A', 'C', 'T', 'G']:
            if base != reference_str[i]:
                positions_ref_bases.append([i+1, reference_str[i], base])
        temp_df = pd.DataFrame(positions_ref_bases, columns=['Pos', 'Ref', 'Base'])
        temp_df['temp'] = 1
        relevant_genes_df['temp'] = 1
        temp_df = pd.merge(temp_df, relevant_genes_df, how='right', on='temp')
        positions_dfs.append(temp_df)
    positions_df = pd.concat(positions_dfs)
    positions_df.Pos = positions_df.Pos.astype(float)
    
    # create mutations dict per mutation
    positions_df['mutation_type'] = positions_df.apply(mutation_type_per_row, axis=1)
    positions_df.to_csv('/sternadi/nobackup/volume1/noam/cmv/temp_positions_df_new.csv')
    return
    
def mutation_type_per_row(row):
    if row.complement == False:
        if row.sequence[int(row.Pos - row.start)] != row.Ref:
            return 'error 1'
        else:
            new_seq = list(row.sequence)
            new_seq[int(row.Pos - row.start)] = row.Base
            new_seq = "".join(new_seq)
            full_seq = row.full_sequence.replace(row.sequence, new_seq)
            new_translated = translate_seq(full_seq)
            diffs = [str(row.translated[i]) + str(i + 1) + str(new_translated[i]) for i in range(len(row.translated)) if row.translated[i] != new_translated[i]]
            return diffs
    else:
        if row.sequence[int(-1 * (row.Pos - row.end + 1))] != str(Seq(row.Ref, IUPAC.unambiguous_dna).reverse_complement()):
            return 'error 2'
        else:
            new_seq = list(row.sequence)
            new_seq[int(-1 * (row.Pos - row.end + 1))] = str(Seq(row.Base, IUPAC.unambiguous_dna).reverse_complement())
            new_seq = "".join(new_seq)
            full_seq = row.full_sequence.replace(row.sequence, new_seq)
            new_translated = translate_seq(full_seq)
            diffs = [str(row.translated[i]) + str(i + 1) + str(new_translated[i]) for i in range(len(row.translated)) if row.translated[i] != new_translated[i]]
            return diffs

create_mutation_type_df()