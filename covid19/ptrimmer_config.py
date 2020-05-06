import os
import sys
sys.path.append('/sternadi/home/volume1/shared/SternLab')
sys.path.append('X:/volume2/noam/SternLab')
from blast_utilities import blast_to_df
import pandas as pd


blasts = blast_to_df('X:/volume2/noam/covid/artic_amplicons/artic_fastas.blast')
seqs = pd.read_excel('X:/volume2/noam/covid/artic_amplicons/primersv3_arctic.xlsx')

blasts['amplicon'] = blasts.read.str.replace('nCoV-2019_', '').str.split('_').str[0]
pairs = pd.merge(blasts[blasts.strand == 'plus'], blasts[blasts.strand == 'minus'], on='amplicon')
pairs['amplicon_length'] = pairs.end_read_y - pairs.end_read_x

pairs = pd.merge(pairs, seqs[['name', 'seq']], left_on='read_x', right_on='name').rename(columns={'seq':'seq_x'})
pairs = pd.merge(pairs, seqs[['name', 'seq']], left_on='read_y', right_on='name').rename(columns={'seq':'seq_y'})
pairs['primer_pairs'] = pairs.read_x + '.' + pairs.read_y
pairs[['seq_x', 'seq_y', 'amplicon_length', 'primer_pairs']].to_csv('X:/volume2/noam/covid/artic_amplicons/ptrimmer_primers.csv', header=False, index=False, sep='\t')



#commands for running ptrimmer on all library when donwloaded from basespace
#for d in /sternadi/datasets/volume2/coronaTech1_20200415/*; do dir_name=$(basename $d); mkdir /sternadi/nobackup/volume1/noam/covid_data/coronaTech1_20200415/ptrimmer_cleanup/$dir_name; echo $d >> /sternadi/nobackup/volume1/noam/covid_data/coronaTech1_20200415/ptrimmer_cleanup/ptrimmer.log; /sternadi/home/volume1/shared/tools/pTrimmer/pTrimmer-1.3.1 -s pair -a /sternadi/home/volume2/noam/covid/artic_amplicons/ptrimmer_primers.txt --read1 $d/*R1* --read2 $d/*R2* -o /sternadi/nobackup/volume1/noam/covid_data/coronaTech1_20200415/ptrimmer_cleanup/$dir_name &>> /sternadi/nobackup/volume1/noam/covid_data/coronaTech1_20200415/ptrimmer_cleanup/ptrimmer.log; done
# for file in */*.fq; do mv "$file" "${file%.*}.fastq"; done

# pipeline
# python /sternadi/home/volume2/noam/SternLab/Python_pipeline/Project_Runner.py -o /sternadi/nobackup/volume1/noam/covid_data/coronaTech1_20200415/python_pipeline_x1_c0_after_ptrimmer/ -i /sternadi/nobackup/volume1/noam/covid_data/coronaTech1_20200415/ptrimmer_cleanup/ -r /sternadi/home/volume2/noam/covid/MN908947.fasta -x 1 -c 0



#commands for running ptrimmer on all library when donwloaded from ftp
#for d in /sternadi/datasets/volume2/coronaTech2_20200427/200427_M00654_0001_000000000-CNBNR/Raw_data/*R1*; do dir_name=$(basename $d); dir_name=$(echo $dir_name | cut -d"_" -f1); mkdir /sternadi/nobackup/volume1/noam/covid_data/coronaTech2_20200427/ptrimmer_cleanup/$dir_name; echo $d >> /sternadi/nobackup/volume1/noam/covid_data/coronaTech2_20200427/ptrimmer_cleanup/ptrimmer.log; /sternadi/home/volume1/shared/tools/pTrimmer/pTrimmer-1.3.1 -s pair -a /sternadi/home/volume2/noam/covid/artic_amplicons/ptrimmer_primers.txt --read1 $d --read2 "${d//_R1_/_R2_}" -o /sternadi/nobackup/volume1/noam/covid_data/coronaTech2_20200427/ptrimmer_cleanup/$dir_name &>> /sternadi/nobackup/volume1/noam/covid_data/coronaTech1_20200415/ptrimmer_cleanup/ptrimmer.log; done
# for file in */*.fq; do mv "$file" "${file%.*}.fastq"; done
# for file in *; do mv "$file" "${file%.*}_L001"; done

# pipeline
# python /sternadi/home/volume2/noam/SternLab/Python_pipeline/Project_Runner.py -o /sternadi/nobackup/volume1/noam/covid_data/coronaTech2_20200427/python_pipeline_x1_c0_after_ptrimmer -i /sternadi/nobackup/volume1/noam/covid_data/coronaTech2_20200427/ptrimmer_cleanup/ -r /sternadi/home/volume2/noam/covid/MN908947.fasta -x 1 -c 0




#technion3
#for d in /sternadi/datasets/volume2/coronaTech3_20200504/Raw_data/*R1*; do dir_name=$(basename $d); dir_name=$(echo $dir_name | cut -d"_" -f1); mkdir /sternadi/nobackup/volume1/noam/covid_data/coronaTech3_20200504/ptrimmer_cleanup/$dir_name; echo $d >> /sternadi/nobackup/volume1/noam/covid_data/coronaTech3_20200504/ptrimmer.log; /sternadi/home/volume1/shared/tools/pTrimmer/pTrimmer-1.3.1 -s pair -a /sternadi/home/volume2/noam/covid/artic_amplicons/ptrimmer_primers.txt --read1 $d --read2 "${d//_R1_/_R2_}" -o /sternadi/nobackup/volume1/noam/covid_data/coronaTech3_20200504/ptrimmer_cleanup/$dir_name &>> /sternadi/nobackup/volume1/noam/covid_data/coronaTech3_20200504/ptrimmer.log; done
# for file in */*.fq; do mv "$file" "${file%.*}.fastq"; done
# for file in *; do mv "$file" "${file%.*}_L001"; done

# pipeline
#python /sternadi/home/volume2/noam/SternLab/Python_pipeline/Project_Runner.py -o /sternadi/nobackup/volume1/noam/covid_data/coronaTech3_20200504/python_pipeline_x1_c0_after_ptrimmer/ -i /sternadi/nobackup/volume1/noam/covid_data/coronaTech3_20200504/ptrimmer_cleanup/ -r /sternadi/home/volume2/noam/covid/MN908947.fasta -x 1 -c 0
#python /sternadi/home/volume2/noam/SternLab/Python_pipeline/Project_Runner.py -o /sternadi/nobackup/volume1/noam/covid_data/coronaTech3_20200504/python_pipeline_x1_c0_v10-9_after_ptrimmer/ -i /sternadi/nobackup/volume1/noam/covid_data/coronaTech3_20200504/ptrimmer_cleanup/ -r /sternadi/home/volume2/noam/covid/MN908947.fasta -x 1 -c 0 -v 1e-9
