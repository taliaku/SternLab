import sys
sys.path.append('/sternadi/home/volume2/noam/SternLab/')
from pbs_runners import array_script_runner



#for f in ./*; do split $f $f -l 10000 --numeric-suffixes=1 -a 4 --additional-suffix .fastq; done

#pS2
array_script_runner('module load python/python-anaconda3.2019.7;\
             PYTHONPATH=":/sternadi/home/volume1/maozgelbart/barcode-utils";\
             printf -v p "%04d" $PBS_ARRAY_INDEX;\
             python -m barcode_utils /sternadi/nobackup/volume1/noam/hcv_data/180423_TMS2-74068001_R1_primer_id/HCV-PS2_S12_L001_R1_001.fastq$p.fastq /sternadi/home/volume2/noam/hcv/references/HCV_primer_id.txt -o /sternadi/nobackup/volume1/noam/hcv_data/180423_TMS2-74068001_R1_primer_id/HCV-PS2_S12_L001_R1_001.fastq$p.barcodes.txt',
             61, alias='barcoding_PS2', queue='hugemem')

#p1
array_script_runner('module load python/python-anaconda3.2019.7;\
             PYTHONPATH=":/sternadi/home/volume1/maozgelbart/barcode-utils";\
             printf -v p "%04d" $PBS_ARRAY_INDEX;\
             python -m barcode_utils /sternadi/nobackup/volume1/noam/hcv_data/180423_TMS2-74068001_R1_primer_id/HCV-P1_S1_L001_R1_001.fastq$p.fastq /sternadi/home/volume2/noam/hcv/references/HCV_primer_id.txt -o /sternadi/nobackup/volume1/noam/hcv_data/180423_TMS2-74068001_R1_primer_id/HCV-P1_S1_L001_R1_001.fastq$p.barcodes.txt',
             62, alias='barcoding_P1', queue='hugemem')

#p2
array_script_runner('module load python/python-anaconda3.2019.7;\
             PYTHONPATH=":/sternadi/home/volume1/maozgelbart/barcode-utils";\
             printf -v p "%04d" $PBS_ARRAY_INDEX;\
             python -m barcode_utils /sternadi/nobackup/volume1/noam/hcv_data/180423_TMS2-74068001_R1_primer_id/HCV-P2_S2_L001_R1_001.fastq$p.fastq /sternadi/home/volume2/noam/hcv/references/HCV_primer_id.txt -o /sternadi/nobackup/volume1/noam/hcv_data/180423_TMS2-74068001_R1_primer_id/HCV-P2_S2_L001_R1_001.fastq$p.barcodes.txt',
             61, alias='barcoding_P2', queue='hugemem')

#p3
array_script_runner('module load python/python-anaconda3.2019.7;\
             PYTHONPATH=":/sternadi/home/volume1/maozgelbart/barcode-utils";\
             printf -v p "%04d" $PBS_ARRAY_INDEX;\
             python -m barcode_utils /sternadi/nobackup/volume1/noam/hcv_data/180423_TMS2-74068001_R1_primer_id/HCV-P3_S3_L001_R1_001.fastq$p.fastq /sternadi/home/volume2/noam/hcv/references/HCV_primer_id.txt -o /sternadi/nobackup/volume1/noam/hcv_data/180423_TMS2-74068001_R1_primer_id/HCV-P3_S3_L001_R1_001.fastq$p.barcodes.txt',
             60, alias='barcoding_P3', queue='hugemem')

#p4
array_script_runner('module load python/python-anaconda3.2019.7;\
             PYTHONPATH=":/sternadi/home/volume1/maozgelbart/barcode-utils";\
             printf -v p "%04d" $PBS_ARRAY_INDEX;\
             python -m barcode_utils /sternadi/nobackup/volume1/noam/hcv_data/180423_TMS2-74068001_R1_primer_id/HCV-P4_S4_L001_R1_001.fastq$p.fastq /sternadi/home/volume2/noam/hcv/references/HCV_primer_id.txt -o /sternadi/nobackup/volume1/noam/hcv_data/180423_TMS2-74068001_R1_primer_id/HCV-P4_S4_L001_R1_001.fastq$p.barcodes.txt',
             60, alias='barcoding_P4', queue='hugemem')


#p5
array_script_runner('module load python/python-anaconda3.2019.7;\
             PYTHONPATH=":/sternadi/home/volume1/maozgelbart/barcode-utils";\
             printf -v p "%04d" $PBS_ARRAY_INDEX;\
             python -m barcode_utils /sternadi/nobackup/volume1/noam/hcv_data/180423_TMS2-74068001_R1_primer_id/HCV-P5_S5_L001_R1_001.fastq$p.fastq /sternadi/home/volume2/noam/hcv/references/HCV_primer_id.txt -o /sternadi/nobackup/volume1/noam/hcv_data/180423_TMS2-74068001_R1_primer_id/HCV-P5_S5_L001_R1_001.fastq$p.barcodes.txt',
             57, alias='barcoding_P5', queue='hugemem')

#p6 ##7
array_script_runner('module load python/python-anaconda3.2019.7;\
             PYTHONPATH=":/sternadi/home/volume1/maozgelbart/barcode-utils";\
             printf -v p "%04d" $PBS_ARRAY_INDEX;\
             python -m barcode_utils /sternadi/nobackup/volume1/noam/hcv_data/180423_TMS2-74068001_R1_primer_id/HCV-P6_S6_L001_R1_001.fastq$p.fastq /sternadi/home/volume2/noam/hcv/references/HCV_primer_id.txt -o /sternadi/nobackup/volume1/noam/hcv_data/180423_TMS2-74068001_R1_primer_id/HCV-P6_S6_L001_R1_001.fastq$p.barcodes.txt',
             60, alias='barcoding_P6', queue='hugemem')

#p7
array_script_runner('module load python/python-anaconda3.2019.7;\
             PYTHONPATH=":/sternadi/home/volume1/maozgelbart/barcode-utils";\
             printf -v p "%04d" $PBS_ARRAY_INDEX;\
             python -m barcode_utils /sternadi/nobackup/volume1/noam/hcv_data/180423_TMS2-74068001_R1_primer_id/HCV-P7_S7_L001_R1_001.fastq$p.fastq /sternadi/home/volume2/noam/hcv/references/HCV_primer_id.txt -o /sternadi/nobackup/volume1/noam/hcv_data/180423_TMS2-74068001_R1_primer_id/HCV-P7_S7_L001_R1_001.fastq$p.barcodes.txt',
             73, alias='barcoding_P7', queue='hugemem')

#p8
array_script_runner('module load python/python-anaconda3.2019.7;\
             PYTHONPATH=":/sternadi/home/volume1/maozgelbart/barcode-utils";\
             printf -v p "%04d" $PBS_ARRAY_INDEX;\
             python -m barcode_utils /sternadi/nobackup/volume1/noam/hcv_data/180423_TMS2-74068001_R1_primer_id/HCV-P8_S8_L001_R1_001.fastq$p.fastq /sternadi/home/volume2/noam/hcv/references/HCV_primer_id.txt -o /sternadi/nobackup/volume1/noam/hcv_data/180423_TMS2-74068001_R1_primer_id/HCV-P8_S8_L001_R1_001.fastq$p.barcodes.txt',
             65, alias='barcoding_P8', queue='hugemem')

#p9
array_script_runner('module load python/python-anaconda3.2019.7;\
             PYTHONPATH=":/sternadi/home/volume1/maozgelbart/barcode-utils";\
             printf -v p "%04d" $PBS_ARRAY_INDEX;\
             python -m barcode_utils /sternadi/nobackup/volume1/noam/hcv_data/180423_TMS2-74068001_R1_primer_id/HCV-P9_S9_L001_R1_001.fastq$p.fastq /sternadi/home/volume2/noam/hcv/references/HCV_primer_id.txt -o /sternadi/nobackup/volume1/noam/hcv_data/180423_TMS2-74068001_R1_primer_id/HCV-P9_S9_L001_R1_001.fastq$p.barcodes.txt',
             65, alias='barcoding_P9', queue='hugemem')

# p10
array_script_runner('module load python/python-anaconda3.2019.7;\
             PYTHONPATH=":/sternadi/home/volume1/maozgelbart/barcode-utils";\
             printf -v p "%04d" $PBS_ARRAY_INDEX;\
             python -m barcode_utils /sternadi/nobackup/volume1/noam/hcv_data/180423_TMS2-74068001_R1_primer_id/HCV-P10_S10_L001_R1_001.fastq$p.fastq /sternadi/home/volume2/noam/hcv/references/HCV_primer_id.txt -o /sternadi/nobackup/volume1/noam/hcv_data/180423_TMS2-74068001_R1_primer_id/HCV-P10_S10_L001_R1_001.fastq$p.barcodes.txt',
             56, alias='barcoding_P10', queue='hugemem')

#p11
array_script_runner('module load python/python-anaconda3.2019.7;\
             PYTHONPATH=":/sternadi/home/volume1/maozgelbart/barcode-utils";\
             printf -v p "%04d" $PBS_ARRAY_INDEX;\
             python -m barcode_utils /sternadi/nobackup/volume1/noam/hcv_data/180423_TMS2-74068001_R1_primer_id/HCV-P11_S11_L001_R1_001.fastq$p.fastq /sternadi/home/volume2/noam/hcv/references/HCV_primer_id.txt -o /sternadi/nobackup/volume1/noam/hcv_data/180423_TMS2-74068001_R1_primer_id/HCV-P11_S11_L001_R1_001.fastq$p.barcodes.txt',
             44, alias='barcoding_P11', queue='hugemem')

#p12
array_script_runner('module load python/python-anaconda3.2019.7;\
             PYTHONPATH=":/sternadi/home/volume1/maozgelbart/barcode-utils";\
             printf -v p "%04d" $PBS_ARRAY_INDEX;\
             python -m barcode_utils /sternadi/nobackup/volume1/noam/hcv_data/180423_TMS2-74068001_R1_primer_id/HCV-P12_S22_L001_R1_001.fastq$p.fastq /sternadi/home/volume2/noam/hcv/references/HCV_primer_id.txt -o /sternadi/nobackup/volume1/noam/hcv_data/180423_TMS2-74068001_R1_primer_id/HCV-P12_S22_L001_R1_001.fastq$p.barcodes.txt',
             2, alias='barcoding_P12', queue='hugemem')


#p2-1
array_script_runner('module load python/python-anaconda3.2019.7;\
             PYTHONPATH=":/sternadi/home/volume1/maozgelbart/barcode-utils";\
             printf -v p "%04d" $PBS_ARRAY_INDEX;\
             python -m barcode_utils /sternadi/nobackup/volume1/noam/hcv_data/180423_TMS2-74068001_R1_primer_id/HCV-P2-1_S13_L001_R1_001.fastq$p.fastq /sternadi/home/volume2/noam/hcv/references/HCV_primer_id.txt -o /sternadi/nobackup/volume1/noam/hcv_data/180423_TMS2-74068001_R1_primer_id/HCV-P2-1_S13_L001_R1_001.fastq$p.barcodes.txt',
             60, alias='barcoding_P2-1', queue='hugemem')
#p11-1
array_script_runner('module load python/python-anaconda3.2019.7;\
             PYTHONPATH=":/sternadi/home/volume1/maozgelbart/barcode-utils";\
             printf -v p "%04d" $PBS_ARRAY_INDEX;\
             python -m barcode_utils /sternadi/nobackup/volume1/noam/hcv_data/180423_TMS2-74068001_R1_primer_id/HCV-P11-1_S19_L001_R1_001.fastq$p.fastq /sternadi/home/volume2/noam/hcv/references/HCV_primer_id.txt -o /sternadi/nobackup/volume1/noam/hcv_data/180423_TMS2-74068001_R1_primer_id/HCV-P11-1_S19_L001_R1_001.fastq$p.barcodes.txt',
             48, alias='barcoding_P11-1', queue='hugemem')

#p3-1
array_script_runner('module load python/python-anaconda3.2019.7;\
             PYTHONPATH=":/sternadi/home/volume1/maozgelbart/barcode-utils";\
             printf -v p "%04d" $PBS_ARRAY_INDEX;\
             python -m barcode_utils /sternadi/nobackup/volume1/noam/hcv_data/180423_TMS2-74068001_R1_primer_id/HCV-P3-1_S14_L001_R1_001.fastq$p.fastq /sternadi/home/volume2/noam/hcv/references/HCV_primer_id.txt -o /sternadi/nobackup/volume1/noam/hcv_data/180423_TMS2-74068001_R1_primer_id/HCV-P3-1_S14_L001_R1_001.fastq$p.barcodes.txt',
             57, alias='barcoding_P3-1', queue='hugemem')

#p5-1
array_script_runner('module load python/python-anaconda3.2019.7;\
             PYTHONPATH=":/sternadi/home/volume1/maozgelbart/barcode-utils";\
             printf -v p "%04d" $PBS_ARRAY_INDEX;\
             python -m barcode_utils /sternadi/nobackup/volume1/noam/hcv_data/180423_TMS2-74068001_R1_primer_id/HCV-P5-1_S15_L001_R1_001.fastq$p.fastq /sternadi/home/volume2/noam/hcv/references/HCV_primer_id.txt -o /sternadi/nobackup/volume1/noam/hcv_data/180423_TMS2-74068001_R1_primer_id/HCV-P5-1_S15_L001_R1_001.fastq$p.barcodes.txt',
             59, alias='barcoding_P5-1', queue='hugemem')

#p6-1
array_script_runner('module load python/python-anaconda3.2019.7;\
             PYTHONPATH=":/sternadi/home/volume1/maozgelbart/barcode-utils";\
             printf -v p "%04d" $PBS_ARRAY_INDEX;\
             python -m barcode_utils /sternadi/nobackup/volume1/noam/hcv_data/180423_TMS2-74068001_R1_primer_id/HCV-P6-1_S16_L001_R1_001.fastq$p.fastq /sternadi/home/volume2/noam/hcv/references/HCV_primer_id.txt -o /sternadi/nobackup/volume1/noam/hcv_data/180423_TMS2-74068001_R1_primer_id/HCV-P6-1_S16_L001_R1_001.fastq$p.barcodes.txt',
             60, alias='barcoding_P6-1', queue='hugemem')

#p7-2
array_script_runner('module load python/python-anaconda3.2019.7;\
             PYTHONPATH=":/sternadi/home/volume1/maozgelbart/barcode-utils";\
             printf -v p "%04d" $PBS_ARRAY_INDEX;\
             python -m barcode_utils /sternadi/nobackup/volume1/noam/hcv_data/180423_TMS2-74068001_R1_primer_id/HCV-P7-2_S17_L001_R1_001.fastq$p.fastq /sternadi/home/volume2/noam/hcv/references/HCV_primer_id.txt -o /sternadi/nobackup/volume1/noam/hcv_data/180423_TMS2-74068001_R1_primer_id/HCV-P7-2_S17_L001_R1_001.fastq$p.barcodes.txt',
             57, alias='barcoding_P7-2', queue='hugemem')

#p8-1
array_script_runner('module load python/python-anaconda3.2019.7;\
             PYTHONPATH=":/sternadi/home/volume1/maozgelbart/barcode-utils";\
             printf -v p "%04d" $PBS_ARRAY_INDEX;\
             python -m barcode_utils /sternadi/nobackup/volume1/noam/hcv_data/180423_TMS2-74068001_R1_primer_id/HCV-P8-1_S20_L001_R1_001.fastq$p.fastq /sternadi/home/volume2/noam/hcv/references/HCV_primer_id.txt -o /sternadi/nobackup/volume1/noam/hcv_data/180423_TMS2-74068001_R1_primer_id/HCV-P8-1_S20_L001_R1_001.fastq$p.barcodes.txt',
             60, alias='barcoding_P8-1', queue='hugemem')

#p8-2
array_script_runner('module load python/python-anaconda3.2019.7;\
             PYTHONPATH=":/sternadi/home/volume1/maozgelbart/barcode-utils";\
             printf -v p "%04d" $PBS_ARRAY_INDEX;\
             python -m barcode_utils /sternadi/nobackup/volume1/noam/hcv_data/180423_TMS2-74068001_R1_primer_id/HCV-P8-2_S21_L001_R1_001.fastq$p.fastq /sternadi/home/volume2/noam/hcv/references/HCV_primer_id.txt -o /sternadi/nobackup/volume1/noam/hcv_data/180423_TMS2-74068001_R1_primer_id/HCV-P8-2_S21_L001_R1_001.fastq$p.barcodes.txt',
             54, alias='barcoding_P8-2', queue='hugemem')

#p9-2
array_script_runner('module load python/python-anaconda3.2019.7;\
             PYTHONPATH=":/sternadi/home/volume1/maozgelbart/barcode-utils";\
             printf -v p "%04d" $PBS_ARRAY_INDEX;\
             python -m barcode_utils /sternadi/nobackup/volume1/noam/hcv_data/180423_TMS2-74068001_R1_primer_id/HCV-P9-2_S18_L001_R1_001.fastq$p.fastq /sternadi/home/volume2/noam/hcv/references/HCV_primer_id.txt -o /sternadi/nobackup/volume1/noam/hcv_data/180423_TMS2-74068001_R1_primer_id/HCV-P9-2_S18_L001_R1_001.fastq$p.barcodes.txt',
             65, alias='barcoding_P9-2', queue='hugemem')



# filter problematic samples by aligned
#cat /sternadi/nobackup/volume1/noam/hcv_data/180423_TMS2-74068001_pipeline_optimized4/HCV-P9/tmp/HCV-P9_S9_L001_001.part*.fasta.blast | cut -f1 > /sternadi/nobackup/volume1/noam/hcv_data/180423_TMS2-74068001_R1_primer_id_filtered/HCV-P9.read_ids.txt
#cat /sternadi/nobackup/volume1/noam/hcv_data/180423_TMS2-74068001_pipeline_optimized4/HCV-P7/tmp/HCV-P7*.fasta.blast | cut -f1 > /sternadi/nobackup/volume1/noam/hcv_data/180423_TMS2-74068001_R1_primer_id_filtered/HCV-P7.read_ids.txt
#cat /sternadi/nobackup/volume1/noam/hcv_data/180423_TMS2-74068001_pipeline_optimized4/HCV-P8-2/tmp/HCV-P8-2*.fasta.blast | cut -f1 > /sternadi/nobackup/volume1/noam/hcv_data/180423_TMS2-74068001_R1_primer_id_filtered/HCV-P8-2.read_ids.txt
#cat /sternadi/nobackup/volume1/noam/hcv_data/180423_TMS2-74068001_pipeline_optimized4/HCV-P7-2/tmp/HCV-P7-2*.fasta.blast | cut -f1 > /sternadi/nobackup/volume1/noam/hcv_data/180423_TMS2-74068001_R1_primer_id_filtered/HCV-P7-2.read_ids.txt
#cat /sternadi/nobackup/volume1/noam/hcv_data/180423_TMS2-74068001_pipeline_optimized4/HCV-P8-1/tmp/HCV-P8-1*.fasta.blast | cut -f1 > /sternadi/nobackup/volume1/noam/hcv_data/180423_TMS2-74068001_R1_primer_id_filtered/HCV-P8-1.read_ids.txt

# filterbyname.sh in=/sternadi/nobackup/volume1/noam/hcv_data/180423_TMS2-74068001_R1_primer_id/HCV-P9_S9_L001_R1_001.fastq out=/sternadi/nobackup/volume1/noam/hcv_data/180423_TMS2-74068001_R1_primer_id_filtered/HCV-P9_R1.filtered.fastq names=/sternadi/nobackup/volume1/noam/hcv_data/180423_TMS2-74068001_R1_primer_id_filtered/HCV-P9.read_ids.txt include=t
# filterbyname.sh in=/sternadi/nobackup/volume1/noam/hcv_data/180423_TMS2-74068001_R1_primer_id/HCV-P7_S7_L001_R1_001.fastq out=/sternadi/nobackup/volume1/noam/hcv_data/180423_TMS2-74068001_R1_primer_id_filtered/HCV-P7_R1.filtered.fastq names=/sternadi/nobackup/volume1/noam/hcv_data/180423_TMS2-74068001_R1_primer_id_filtered/HCV-P7.read_ids.txt include=t
# filterbyname.sh in=/sternadi/nobackup/volume1/noam/hcv_data/180423_TMS2-74068001_R1_primer_id/HCV-P7-2_S17_L001_R1_001.fastq out=/sternadi/nobackup/volume1/noam/hcv_data/180423_TMS2-74068001_R1_primer_id_filtered/HCV-P7-2_R1.filtered.fastq names=/sternadi/nobackup/volume1/noam/hcv_data/180423_TMS2-74068001_R1_primer_id_filtered/HCV-P7-2.read_ids.txt include=t
# filterbyname.sh in=/sternadi/nobackup/volume1/noam/hcv_data/180423_TMS2-74068001_R1_primer_id/HCV-P8-1_S20_L001_R1_001.fastq out=/sternadi/nobackup/volume1/noam/hcv_data/180423_TMS2-74068001_R1_primer_id_filtered/HCV-P8-1_R1.filtered.fastq names=/sternadi/nobackup/volume1/noam/hcv_data/180423_TMS2-74068001_R1_primer_id_filtered/HCV-P8-1.read_ids.txt include=t
# filterbyname.sh in=/sternadi/nobackup/volume1/noam/hcv_data/180423_TMS2-74068001_R1_primer_id/HCV-P8-2_S21_L001_R1_001.fastq out=/sternadi/nobackup/volume1/noam/hcv_data/180423_TMS2-74068001_R1_primer_id_filtered/HCV-P8-2_R1.filtered.fastq names=/sternadi/nobackup/volume1/noam/hcv_data/180423_TMS2-74068001_R1_primer_id_filtered/HCV-P8-2.read_ids.txt include=t


#for f in ./*.fastq; do split $f $f -l 10000 --numeric-suffixes=1 -a 4 --additional-suffix .fastq; done

#p7
array_script_runner('module load python/python-anaconda3.2019.7;\
             PYTHONPATH=":/sternadi/home/volume1/maozgelbart/barcode-utils";\
             printf -v p "%04d" $PBS_ARRAY_INDEX;\
             python -m barcode_utils /sternadi/nobackup/volume1/noam/hcv_data/180423_TMS2-74068001_R1_primer_id_filtered/HCV-P7_R1.filtered.fastq$p.fastq /sternadi/home/volume2/noam/hcv/references/HCV_primer_id.txt -o /sternadi/nobackup/volume1/noam/hcv_data/180423_TMS2-74068001_R1_primer_id_filtered/HCV-P7_R1.filtered.fastq$p.barcodes.txt',
             60, alias='barcoding_P7', queue='hugemem')

#p9
array_script_runner('module load python/python-anaconda3.2019.7;\
             PYTHONPATH=":/sternadi/home/volume1/maozgelbart/barcode-utils";\
             printf -v p "%04d" $PBS_ARRAY_INDEX;\
             python -m barcode_utils /sternadi/nobackup/volume1/noam/hcv_data/180423_TMS2-74068001_R1_primer_id_filtered/HCV-P9_R1.filtered.fastq$p.fastq /sternadi/home/volume2/noam/hcv/references/HCV_primer_id.txt -o /sternadi/nobackup/volume1/noam/hcv_data/180423_TMS2-74068001_R1_primer_id_filtered/HCV-P9_R1.filtered.fastq$p.barcodes.txt',
             49, alias='barcoding_P9', queue='hugemem')

#p7-2
array_script_runner('module load python/python-anaconda3.2019.7;\
             PYTHONPATH=":/sternadi/home/volume1/maozgelbart/barcode-utils";\
             printf -v p "%04d" $PBS_ARRAY_INDEX;\
             python -m barcode_utils /sternadi/nobackup/volume1/noam/hcv_data/180423_TMS2-74068001_R1_primer_id_filtered/HCV-P7-2_R1.filtered.fastq$p.fastq /sternadi/home/volume2/noam/hcv/references/HCV_primer_id.txt -o /sternadi/nobackup/volume1/noam/hcv_data/180423_TMS2-74068001_R1_primer_id_filtered/HCV-P7-2_R1.filtered.fastq$p.barcodes.txt',
             44, alias='barcoding_P7-2', queue='hugemem')

#p8-1
array_script_runner('module load python/python-anaconda3.2019.7;\
             PYTHONPATH=":/sternadi/home/volume1/maozgelbart/barcode-utils";\
             printf -v p "%04d" $PBS_ARRAY_INDEX;\
             python -m barcode_utils /sternadi/nobackup/volume1/noam/hcv_data/180423_TMS2-74068001_R1_primer_id_filtered/HCV-P8-1_R1.filtered.fastq$p.fastq /sternadi/home/volume2/noam/hcv/references/HCV_primer_id.txt -o /sternadi/nobackup/volume1/noam/hcv_data/180423_TMS2-74068001_R1_primer_id_filtered/HCV-P8-1_R1.filtered.fastq$p.barcodes.txt',
             2, alias='barcoding_P8-1', queue='hugemem')

#p8-2
array_script_runner('module load python/python-anaconda3.2019.7;\
             PYTHONPATH=":/sternadi/home/volume1/maozgelbart/barcode-utils";\
             printf -v p "%04d" $PBS_ARRAY_INDEX;\
             python -m barcode_utils /sternadi/nobackup/volume1/noam/hcv_data/180423_TMS2-74068001_R1_primer_id_filtered/HCV-P8-2_R1.filtered.fastq$p.fastq /sternadi/home/volume2/noam/hcv/references/HCV_primer_id.txt -o /sternadi/nobackup/volume1/noam/hcv_data/180423_TMS2-74068001_R1_primer_id_filtered/HCV-P8-2_R1.filtered.fastq$p.barcodes.txt',
             31, alias='barcoding_P8-2', queue='hugemem')