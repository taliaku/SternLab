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
import os
import sys
sys.path.append('/sternadi/home/volume2/noam/SternLab')
sys.path.append(r'X:\volume2\noam\Sternlab')
from blast_utilities import blast_to_mutations_list, blast_to_df

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
    mutations = pd.merge(mutations, f_df, on=['Pos', 'Ref', 'Base'], how='left')
    #mutations = mutations[['Ref', 'Pos', 'Base', 'Freq_' + str(file1), 'Freq_' + str(file2), 'Read_count_' + str(file1), 'Read_count_' + str(file2), 'file_' + str(file1), 'file_' + str(file2), 'protein', 'mutation_type']]
    mutations = mutations[['Ref', 'Pos', 'Base', 'Freq_' + str(file1), 'Freq_' + str(file2), 'Read_count_' + str(file1), 'Read_count_' + str(file2), 'protein', 'mutation_type']]
    mutations.rename(columns={'mutation_type':'protein_mutation'}, inplace=True)
    conditions = [(mutations.protein_mutation == '[]'), (mutations.protein_mutation.isna()), ((mutations.protein_mutation != '[]') & ~(mutations.protein_mutation.isna()))]
    choices = ['synonymous', 'non_coding', 'non_synonymous']
    mutations['mutation_type'] = np.select(conditions, choices)
    
    deletions = pd.merge(deletions, mutation_type_df[['Pos', 'protein']].drop_duplicates(), on='Pos', how='left')
    #deletions = deletions[['Ref', 'Pos', 'Base', 'Freq_' + str(file1), 'Freq_' + str(file2), 'Read_count_' + str(file1), 'Read_count_' + str(file2), 'file_' + str(file1), 'file_' + str(file2), 'protein']]
    deletions = deletions[['Ref', 'Pos', 'Base', 'Freq_' + str(file1), 'Freq_' + str(file2), 'Read_count_' + str(file1), 'Read_count_' + str(file2), 'protein']]
    # insertions
    #insertions = new_df[((new_df.Ref == '-') & ((new_df.Read_count_x - new_df.Read_count_y).abs() > 20))]    
    # positions appearing in one and not the other
    #weird_pos = pd.concat([dfx[(~(dfx.Pos.isin(dfy.Pos.tolist())))], dfy[(~(dfy.Pos.isin(dfx.Pos.tolist())))]])
    mutations.to_csv(output_dir + 'mutations_' + str(file1) + '_' + str(file2) + '.csv', index=False)
    deletions.to_csv(output_dir + 'deletions_' + str(file1) + '_' + str(file2) + '.csv', index=False)
    #insertions.to_csv(output_dir + 'insertions_' + file1 + '_' + file2 + '.csv', index=False)
    #weird_pos.to_csv(output_dir + 'weird_positions_' + file1 + '_' + file2 + '.csv', index=False)
    return
 
for i in range(1,9):
    for j in range(i+1,9):
        if i != j:
            compare_freqs(freqs, mutation_type, i, j, 'Z:/volume1/noam/cmv/comparisons2/')
    

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








#####
            
to_pivot = pd.merge(freqs, non_synonymous, how='right', on=['Ref', 'Pos', 'Base'])
to_pivot = pd.merge(to_pivot, mutation_type[['Base', 'Ref', 'Pos', 'protein', 'mutation_type']], how='left', on=['Ref', 'Pos', 'Base'])
to_pivot['mutation'] = to_pivot.Ref + to_pivot.Pos.astype(str) + to_pivot.Base
to_pivot = to_pivot.pivot_table(values='Freq', index=['protein', 'mutation'], columns='file')
to_pivot = to_pivot[[3,4,7]]
to_pivot = to_pivot[~(to_pivot.isnull().any(axis=1))]


## mutations
t = to_pivot[(to_pivot[7] > 0.8)].reset_index()
conditions = [(t[3] < 0.5), (t[3] >= 0.5) & (t[3] <= 0.7), (t[3] > 0.7) & (t[3] <= 0.8), (t[3] > 0.8)]
choices = ['less than 0.5', 'between 0.5 and 0.7', 'between 0.7 and 0.8', 'more than 0.8']
t['3_grouped'] = np.select(conditions, choices)

# deletions
freqs[(freqs.Pos.isin(freqs[(freqs.Base == '-') & (freqs.Freq > 0.5) & (freqs.Read_count > 5)].Pos.drop_duplicates().tolist())) & (freqs.Base == '-')].pivot_table(values='Freq', index='Pos', columns='file')


# mapped and unmapped
b = freqs[(freqs.Base == freqs.Ref)].pivot_table(values='Read_count', index='Pos', columns='file')
a = freqs[(freqs.Base == freqs.Ref) & (freqs.Ref != '-')].pivot_table(values='Read_count', index='mutation', columns='file')
a = a[(a).isnull().any(axis=1)]
a = a[a.max(axis=1) > 1]

c = pd.merge(a, mutation_type[['Pos', 'protein']].drop_duplicates(), left_on='mutation', right_on='Pos')


d = pd.merge(c.groupby('protein').Pos.count().reset_index().rename(columns={'Pos':'disappearance_count'}), mutation_type[['protein', 'Pos']].drop_duplicates().groupby('protein').count().reset_index().rename(columns={'Pos':'protein_length'}), on='protein')

# get all mutations for protein RL1
pd.merge(mutation_type[mutation_type.protein =='RL1'], non_synonymous, on=['Pos', 'Ref', 'Base']).sort_values('Pos')
# with frequencies from p16
pd.merge(pd.merge(mutation_type[mutation_type.protein =='RL1'], non_synonymous, on=['Pos', 'Ref', 'Base']).sort_values('Pos'), freqs[freqs.file==3], on=['Pos', 'Ref', 'Base'])







######## controls
nd_p16 = pd.read_csv('Z:/volume1/noam/cmv/cmv_no_medicine/P16ND/P16ND.freqs', '\t')
nd_p28 = pd.read_csv('Z:/volume1/noam/cmv/cmv_no_medicine/P28ND/P28ND.freqs', '\t')
lt_p16 = pd.read_csv('Z:/volume1/noam/cmv/cmv_no_medicine/P16Lt/P16LT.freqs', '\t')

def estimate_insertion_freq(df):
    #read_counts = df[(df.Ref != '-')][['Pos', 'Read_count']].drop_duplicates()
    read_counts = df[(df.Ref != '-')][['Pos', 'Read_count', 'file']].drop_duplicates()
    read_counts.rename(columns={'Read_count':'estimated_read_count', 'Pos':'rounded_pos'}, inplace=True)
    insertions = df[(df.Ref == '-')]
    not_insertions = df[(df.Ref != '-')]
    insertions['rounded_pos'] = insertions.Pos.astype(int).astype(float)
    insertions = pd.merge(insertions, read_counts, how='left', on=['rounded_pos', 'file'])
    insertions['estimated_freq'] = insertions.Freq * insertions.Read_count / insertions.estimated_read_count
    df = pd.concat([insertions, not_insertions])
    return df.sort_values(['Pos'])

a = pd.read_csv('Z:/volume1/noam/cmv/cmv_medicine/all_freqs.csv')
a['file'] = a.file.astype(str)

nd_p16['file'] = 'nd_p16'
nd_p28['file'] = 'nd_p28'
lt_p16['file'] = 'lt_p16'
a = pd.concat([a, nd_p16, nd_p28, lt_p16])
a = estimate_insertion_freq(a)
a['full_mutation'] = a.Ref + a.Pos.astype(str) + a.Base

a.to_csv('Z:/volume1/noam/cmv/all_freqs.csv', index=False)



######

def stats(df):
    df = df[df.Read_count > 5]
    substitutions = df[(df.Base != df.Ref) & (df.Ref != '-') & (df.Base != '-') & (df.Freq > 0.2)]
    deletions = df[(df.Base != df.Ref) & (df.Ref != '-') & (df.Base == '-') & (df.Freq > 0.2)]
    insertions = df[(df.Base != df.Ref) & (df.estimated_freq > 0.2)]
    return len(substitutions), len(deletions), len(insertions)

# substitutions in common, inner merge
# deletions in common, inner merge
# insertions in common, inner merge
pd.merge(a[(a.Read_count > 5) & (a.Base != a.Ref) & (a.estimated_freq > 0.2) & (a.file.str.contains('nd_p16'))][['Ref', 'Pos', 'Base']].drop_duplicates(), a[(a.Read_count > 5) & (a.Base != a.Ref) & (a.estimated_freq > 0.2) & (a.file.str.contains('lt_p16'))][['Ref', 'Pos', 'Base']].drop_duplicates(), how='inner', on=['Ref', 'Pos', 'Base'])    

######

a = pd.read_csv('Z:/volume1/noam/cmv/all_freqs.csv')
mutation_type = pd.read_csv('Z:/volume1/noam/cmv/mutation_type_position_df.csv')

# mutations exceeding 0.2 in one of the three samples and merging them with mutation type
pd.merge(a[(a.file.str.contains('p')) & (a.Base != a.Ref) & (a.Ref != '-') & (a.Base != '-') & (a.Freq > 0.2)][['Pos', 'Ref', 'Base', 'file']].groupby(['Pos', 'Base', 'Ref']).file.apply(','.join).reset_index(), mutation_type[mutation_type.mutation_type != '[]'][['Pos', 'Ref', 'Base', 'protein', 'mutation_type']], how='left', on=['Pos', 'Ref', 'Base'])

b = a[(a.Base != a.Ref) & (a.Ref != '-') & (a.Read_count > 5)].pivot_table(values='Freq', index='full_mutation', columns='file')

b[(b['nd_p16'] < 0.2) & (b.nd_p28 < 0.2) & (b.lt_p16 > 0.2)]


to_pivot = pd.merge(a, non_synonymous, how='right', on=['Ref', 'Pos', 'Base'])
to_pivot = pd.merge(to_pivot, mutation_type[['Base', 'Ref', 'Pos', 'protein', 'mutation_type']], how='left', on=['Ref', 'Pos', 'Base'])
to_pivot['mutation'] = to_pivot.Ref + to_pivot.Pos.astype(str) + to_pivot.Base
to_pivot = to_pivot.pivot_table(values='Freq', index=['protein', 'mutation'], columns='file')
to_pivot = to_pivot[[3,4,7]]
to_pivot = to_pivot[~(to_pivot.isnull().any(axis=1))]


## mutations
t = to_pivot[(to_pivot[7] > 0.8)].reset_index()
conditions = [(t[3] < 0.5), (t[3] >= 0.5) & (t[3] <= 0.7), (t[3] > 0.7) & (t[3] <= 0.8), (t[3] > 0.8)]
choices = ['less than 0.5', 'between 0.5 and 0.7', 'between 0.7 and 0.8', 'more than 0.8']
t['3_grouped'] = np.select(conditions, choices)






###### context - didn't find context dependency

refs = a[a.Ref != '-'][['Pos', 'Ref']].drop_duplicates().sort_values('Pos')
refs['ref_previous_1'] = refs.Ref.shift(1)
refs['ref_previous_2'] = refs.Ref.shift(2)
refs['ref_previous_3'] = refs.Ref.shift(3)
refs['ref_previous_4'] = refs.Ref.shift(4)
refs['ref_next_1'] = refs.Ref.shift(-1)
refs['ref_next_2'] = refs.Ref.shift(-2)
refs['ref_next_3'] = refs.Ref.shift(-3)


pd.merge(a[(a.file == 7) & (a.Base != a.Ref) & (a.Ref != '-') & (a.Base != '-') & (a.Freq > 0.8)], refs, on=['Ref', 'Pos']).groupby(['Ref', 'ref_previous_2']).Pos.count()
