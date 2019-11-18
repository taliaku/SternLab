# -*- coding: utf-8 -*-
"""
Created on Tue Jul 23 16:13:52 2019

@author: Noam
"""

import re
import pandas as pd
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
import numpy as np
import os
import sys


def get_rna_sequence(protein, start, end, complement, reference_str):
    if complement == False:
        return reference_str[start-1:end]
    else:
        return str(Seq(reference_str[end-1:start], IUPAC.unambiguous_dna).reverse_complement())
        
def translate_seq(seq):
    return str(Seq(seq, IUPAC.unambiguous_dna).translate())

 
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




def create_gene_dataframe():
    genes = []
    with open('X:/volume2/noam/ms2_reference/mutation_type_ncbi/ncbi_genes.txt') as f:
        fi = f.read()
    for row in fi.splitlines():
        if not row.startswith('\t') and not row.startswith('>') and row != '':
            print('row' + row)
            protein = fi.split(row)[1].split('gene\t')[1].split('\n')[0]
            start = int(row.split('\t')[0].replace('<', ''))
            end = int(row.split('\t')[1].replace('<', ''))
            complement = start > end
            genes.append((protein, start, end, complement))
        fi = fi.replace(row, '', 1)
    gene_df = pd.DataFrame(genes, columns=['protein', 'start', 'end', 'complement'])
    
    with open('X:/volume2/noam/ms2_reference/noam_concensus_ref.txt') as f:
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
    if len(complement_df) > 0:
        complement_df['full_sequence'] = complement_df.sort_values('start', ascending=False).groupby('protein')['sequence'].transform(lambda x: x.sum())
    
    gene_df = pd.concat([complement_df, not_complement_df])
    gene_df['translated'] = gene_df.full_sequence.apply(translate_seq)
    gene_df.to_csv('X:/volume2/noam/ms2_reference/mutation_type_ncbi/ncbi_genes.csv', index=False)
    return gene_df



def create_mutation_type_df():
    #gene_df = pd.read_csv('Z:/volume1/noam/cmv/ncbi_genes_dataframe.csv')
    #with open('Z:/volume1/noam/cmv/cmv_dana_ref.fasta') as f:
    gene_df = pd.read_csv('X:/volume2/noam/ms2_reference/mutation_type_ncbi/ncbi_genes.csv')
    with open('X:/volume2/noam/ms2_reference/noam_concensus_ref.txt') as f:
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
    positions_df.to_csv('X:/volume2/noam/ms2_reference/mutation_type_ncbi/mutation_types.csv')
    return
   