import random
import pandas as pd
import math

# known
PCR_CYCLES = 10
FINAL_SAMPLING = 0.1
# unknown
AMPLIFICATION_PERCENT_PER_CYCLE = 0.5
RECOMBINATION_PERCENT_PER_CYCLE = 0.1
#depends on sample
MUTANT_PERCENT = 0.1
START_VIRUS_COUNT = 1000
# analysis
CONCENSUS_CUTOFF = 0.9
CONCENSUS_MIN_READS_PER_PID = 0

def run_model():
    viruses = create_viruses(START_VIRUS_COUNT, MUTANT_PERCENT)
    print(analyze_viruses(viruses, CONCENSUS_CUTOFF, CONCENSUS_MIN_READS_PER_PID))
    print('pcr cycles:')
    for i in range(PCR_CYCLES):
        viruses = pcr_cycle(viruses, AMPLIFICATION_PERCENT_PER_CYCLE)
        viruses = recombination_cycle(viruses, RECOMBINATION_PERCENT_PER_CYCLE)
        print(analyze_viruses(viruses, CONCENSUS_CUTOFF, CONCENSUS_MIN_READS_PER_PID))        
    viruses = sample(viruses, FINAL_SAMPLING)
    print('after sampling:')
    print(analyze_viruses(viruses, CONCENSUS_CUTOFF, CONCENSUS_MIN_READS_PER_PID))
    
def create_viruses(start_virus_count, mutant_percent):
    pids = list(range(start_virus_count))
    mutants = pids[:int(start_virus_count * mutant_percent)]
    wt = pids[int(start_virus_count * mutant_percent):]
    mutants = [(i, 'm') for i in mutants]
    wt = [(i, 'w') for i in wt]
    viruses = mutants + wt
    return viruses

def sample(viruses, factor):
    random.shuffle(viruses)
    return viruses[:int(len(viruses) * factor)]
    
def pcr_cycle(viruses, amplification_percent_per_cycle):
    random.shuffle(viruses)
    amplified = viruses[:int(len(viruses) * amplification_percent_per_cycle)]
    not_amplified = viruses[int(len(viruses) * amplification_percent_per_cycle):]
    return (amplified * 2) + not_amplified

def recombination_cycle(viruses, recombination_percent_per_cycle):
    random.shuffle(viruses)
    recombinants = viruses[:int(len(viruses) * recombination_percent_per_cycle)]
    random.shuffle(recombinants)
    non_recombinants = viruses[int(len(viruses) * recombination_percent_per_cycle):]
    after_recombination = []
    for i in range(0, math.floor(len(recombinants) / 2), 2):
        after_recombination.append((recombinants[i][0], recombinants[i+1][1]))
        after_recombination.append((recombinants[i+1][0], recombinants[i][1]))
    return after_recombination + non_recombinants

def analyze_viruses(viruses, concensus_cutoff, concensus_min_reads_per_pid):
    viruses = pd.DataFrame(viruses)
    viruses.columns = ['pid', 'type']
    viruses['read_count_per_pid'] = viruses.groupby('pid')['type'].transform('count')
    viruses['read_count_per_pid_type'] = viruses.groupby(['pid', 'type'])['type'].transform('count')
    viruses['pid_concensus'] = viruses.read_count_per_pid_type / viruses.read_count_per_pid
    viruses = viruses[['pid', 'type', 'pid_concensus', 'read_count_per_pid']].drop_duplicates()
    all_pids = len(viruses[['pid']].drop_duplicates())
    viruses = viruses[(viruses.pid_concensus >= concensus_cutoff)]
    if concensus_min_reads_per_pid:
        viruses = viruses[(viruses.read_count_per_pid >= concensus_min_reads_per_pid)]
    kept_pids = len(viruses[['pid']].drop_duplicates())
    percent_kept_pids = kept_pids / all_pids * 100
    viruses['mutant_percent'] = viruses.groupby('type')['pid'].transform('count') / len(viruses)
    mutant_freq = viruses[['type', 'mutant_percent']].drop_duplicates()
    mutant_freq = float(mutant_freq[(viruses.type == 'm')]['mutant_percent'])
    return mutant_freq, percent_kept_pids
    