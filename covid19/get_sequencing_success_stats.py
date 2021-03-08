import pandas as pd
import re
import glob
from Bio import SeqIO


p1 = re.compile('Total number of reads: (\d*)')
p2 = re.compile('Number of reads mapped to reference: (\d*)')

results = []

for sample in glob.glob('/sternadi/nobackup/volume1/covid/israel_artic_pipeline/*/python_pipeline/*/Summary.txt'):
    with open(sample) as fi:
        fi = fi.read()
    name = sample.replace('/Summary.txt', '').split('/')[-1]
    library = sample.replace('/sternadi/nobackup/volume1/covid/israel_artic_pipeline/', '').split('/')[0]
    try:
        read_count = p1.findall(fi)[0]
        mapped_read_count = p2.findall(fi)[0]
        freqs = pd.read_csv(glob.glob(sample.replace('Summary.txt', '*.freqs.csv'))[0])
        covered_positions_count = len(freqs[(freqs.base == freqs.ref_base) & (freqs.ref_base != '-') & (freqs.coverage > 10)])
        suspicious_30_70_positions_percent = 100 * len(freqs[(freqs.base != freqs.ref_base) & (freqs.ref_base != '-') & (freqs.frequency >= 0.3) & (freqs.frequency <= 0.7) & (freqs.ref_position.isin(range(300,29603)))]) / 29903
        spike_covered_positions_count = len(freqs[(freqs.base == freqs.ref_base) & (freqs.ref_base != '-') & (freqs.coverage > 10) & (freqs.ref_position.isin(range(21563,25385)))])
        median_coverage = freqs[(freqs.base == freqs.ref_base) & (freqs.ref_base != '-') & (freqs.coverage > 10)].coverage.median()
        spike_median_coverage = freqs[(freqs.base == freqs.ref_base) & (freqs.ref_base != '-') & (freqs.coverage > 10) & (freqs.ref_position.isin(range(21563, 25385)))].coverage.median()
        amp74_median_coverage = freqs[(freqs.base == freqs.ref_base) & (freqs.ref_base != '-') & (freqs.coverage > 10) & (freqs.ref_position.isin(range(22346, 22516)))].coverage.median()
        amp74_covered_positions_count = len(freqs[(freqs.base == freqs.ref_base) & (freqs.ref_base != '-') & (freqs.coverage > 10) & (freqs.ref_position.isin(range(22346, 22516)))])
        amp76_median_coverage = freqs[(freqs.base == freqs.ref_base) & (freqs.ref_base != '-') & (freqs.coverage > 10) & (freqs.ref_position.isin(range(22878, 23144)))].coverage.median()
        amp76_covered_positions_count = len(freqs[(freqs.base == freqs.ref_base) & (freqs.ref_base != '-') & (freqs.coverage > 10) & (freqs.ref_position.isin(range(22878, 23144)))])
    except:
        read_count = None
        mapped_read_count = None
        freqs = None
        covered_positions_count = None
        suspicious_30_70_positions_percent = None
        spike_covered_positions_count = None
        median_coverage = None
        spike_median_coverage = None
        amp74_median_coverage = None
        amp74_covered_positions_count = None
        amp76_median_coverage = None
        amp76_covered_positions_count = None
    results.append((name, library, read_count, mapped_read_count, covered_positions_count, suspicious_30_70_positions_percent, spike_covered_positions_count, median_coverage, spike_median_coverage, amp74_median_coverage, amp74_covered_positions_count, amp76_median_coverage, amp76_covered_positions_count))
results = pd.DataFrame(results, columns=['name', 'library', 'read_count', 'mapped_read_count', 'covered_positions_count', 'suspicious_30_70_positions_percent', 'spike_covered_positions_count', 'median_coverage', 'spike_median_coverage', 'amp74_median_coverage', 'amp74_covered_positions_count', 'amp76_median_coverage', 'amp76_covered_positions_count'])
results.fillna(value=pd.np.nan, inplace=True)
results['read_count'] = results.read_count.astype(float)
results['mapped_read_count'] = results.mapped_read_count.astype(float)
results['mapped_read_precent'] = 100 * results.mapped_read_count / results.read_count
results['covered_positions_percent'] = 100 * results.covered_positions_count / 29903
results['name'] = results.name.astype(str)

fasta_strict = list(SeqIO.parse('/sternadi/nobackup/volume1/covid/israel_artic_pipeline/all_results/strict/israel_sequences.fasta', "fasta"))
fasta_majority = list(SeqIO.parse('/sternadi/nobackup/volume1/covid/israel_artic_pipeline/all_results/majority/israel_sequences.fasta', "fasta"))

n_counts_strict = []
for f in fasta_strict:
    non_n_precent = 100 - (f.seq.count('N') * 100 / 29903)
    name = f.id.replace('Israel/', '').replace('/2020', '').replace('/2021', '')
    n_counts_strict.append((name, non_n_precent))
n_counts_strict = pd.DataFrame(n_counts_strict, columns=['name', 'strict_non_n_precent'])

n_counts_majority = []
for f in fasta_majority:
    non_n_precent = 100 - (f.seq.count('N') * 100 / 29903)
    name = f.id.replace('Israel/', '').replace('/2020', '').replace('/2021', '')
    n_counts_majority.append((name, non_n_precent))
n_counts_majority = pd.DataFrame(n_counts_majority, columns=['name', 'majority_non_n_precent'])

n_counts = pd.merge(n_counts_strict, n_counts_majority, on='name')
results = pd.merge(results, n_counts, on='name', how='left')
results = results[~results.library.str.contains('Tech6')]
results = results[['name', 'library', 'read_count', 'mapped_read_precent', 'covered_positions_percent', 'strict_non_n_precent', 'majority_non_n_precent', 'suspicious_30_70_positions_percent', 'median_coverage', 'spike_median_coverage', 'spike_covered_positions_count', 'amp74_median_coverage', 'amp74_covered_positions_count', 'amp76_median_coverage', 'amp76_covered_positions_count']]
results.to_csv('/sternadi/nobackup/volume1/covid/israel_artic_pipeline/sequencing_success_stats.csv', index=False)

