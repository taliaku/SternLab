
import os,sys,inspect
parentdir = os.path.dirname(os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe()))))
sys.path.insert(0,parentdir)
import pandas as pd
import os
from freqs_utilities import estimate_insertion_freq

UNITED_FREQS_FILE = '/sternadi/nobackup/volume1/covid/israel_artic_pipeline/all_results/israel_freqs.csv'
UNITED_CONS_FILE = '/sternadi/nobackup/volume1/covid/israel_artic_pipeline/all_results/israel_sequences.mutations.csv'
ALERTS_OUTPUT = '/sternadi/nobackup/volume1/covid/israel_artic_pipeline/all_results/alerts.txt'
TEMP_FREQS_FILE = '/sternadi/nobackup/volume1/covid/israel_artic_pipeline/all_results/temp_freqs.csv'

########### add mutations here for alerts ##################
# mutations to search for in freqs file as keys, values are frequencies to alert over
freqs_alerts = {'-6696.1T':0.4, '-6696.1A':0.4, '-6696.1C':0.4, '-6696.1G':0.4}
# mutations to search for on the concensus level, values are comments for why this is added
consensus_nuc_alerts = {'-6696.001T':'correlated with negative serology',
                        'C22879.0A':'spike N439K variant',
                        'C22879.0G':'spike N439K variant',
                        'A22920.0T':'mink mutation Y453F',
                        'A21766.0-':'mink deletion, mutation from B.1.1.7',
                        'C21767.0-':'mink deletion, mutation from B.1.1.7',
                        'A21768.0-':'mink deletion, mutation from B.1.1.7',
                        'T21769.0-':'mink deletion, mutation from B.1.1.7',
                        'G21770.0-':'mink deletion, mutation from B.1.1.7',
                        'T21771.0-':'mink deletion',
                        'A23063.0T':'N501Y spike RBD',
                        'C23604.0A':'P681H spike furin cleavage site',
                        'G5230.0T':'mutation from S.501Y.V2',
                        'A10323.0G':'mutation from S.501Y.V2',
                        'C21614.0T':'mutation from S.501Y.V2',
                        'A21801.0C':'mutation from S.501Y.V2',
                        'A22206.0G':'mutation from S.501Y.V2',
                        'T22287.0A':'mutation from S.501Y.V2',
                        'C22286.0-':'mutation from S.501Y.V2',
                        'T22287.0-':'mutation from S.501Y.V2',
                        'T22288.0-':'mutation from S.501Y.V2',
                        'G22289.0-':'mutation from S.501Y.V2',
                        'C22290.0-':'mutation from S.501Y.V2',
                        'T22291.0-':'mutation from S.501Y.V2',
                        'T22292.0-':'mutation from S.501Y.V2',
                        'T22293.0-':'mutation from S.501Y.V2',
                        'A22294.0-':'mutation from S.501Y.V2',
                        'G22299.0T':'mutation from S.501Y.V2',
                        'G22813.0T':'mutation from S.501Y.V2',
                        'G23012.0A':'mutation from S.501Y.V2, E484K reinfection',
                        'A23063.0T':'mutation from S.501Y.V2, B.1.1.7',
                        'C23664.0T':'mutation from S.501Y.V2',
                        #'G25563.0T':'mutation from S.501Y.V2',
                        'C25904.0T':'mutation from S.501Y.V2',
                        'C26456.0T':'mutation from S.501Y.V2',
                        'C28887.0T':'mutation from S.501Y.V2',
                        'C3267.0T':'mutation from B.1.1.7',
                        'C5388.0A':'mutation from B.1.1.7',
                        'T6954.0C':'mutation from B.1.1.7',
                        'T11288.0-':'mutation from B.1.1.7',
                        'C11289.0-':'mutation from B.1.1.7',
                        'T11290.0-':'mutation from B.1.1.7',
                        'G11291.0-':'mutation from B.1.1.7',
                        'G11292.0-':'mutation from B.1.1.7',
                        'T11293.0-':'mutation from B.1.1.7',
                        'T11294.0-':'mutation from B.1.1.7',
                        'T11295.0-':'mutation from B.1.1.7',
                        'T11296.0-':'mutation from B.1.1.7',
                        'T21765.0-':'mutation from B.1.1.7',
                        'T21991.0-':'mutation from B.1.1.7',
                        'T21992.0-':'mutation from B.1.1.7',
                        'A21993.0-':'mutation from B.1.1.7',
                        'C23271.0A':'mutation from B.1.1.7',
                        'C23604.0A':'mutation from B.1.1.7',
                        'C23709.0T':'mutation from B.1.1.7',
                        'T24506.0G':'mutation from B.1.1.7',
                        'G24914.0C':'mutation from B.1.1.7',
                        'C27972.0T':'mutation from B.1.1.7',
                        'G28048.0T':'mutation from B.1.1.7',
                        'A28111.0G':'mutation from B.1.1.7',
                        'G28280.0C':'mutation from B.1.1.7',
                        'A28281.0T':'mutation from B.1.1.7',
                        'T28282.0A':'mutation from B.1.1.7',
                        'C28977.0T':'mutation from B.1.1.7',}

############################################################


if __name__ == '__main__':
    # filter freqs data using grep for memory reasons
    for i in freqs_alerts:
        os.system(f'grep ^{i.split(".")[0][1:]} {UNITED_FREQS_FILE} >> {TEMP_FREQS_FILE}')
    os.system(f'head -n 1 {UNITED_FREQS_FILE} > {TEMP_FREQS_FILE}2')
    os.system(f'cat {TEMP_FREQS_FILE} | sort | uniq >> {TEMP_FREQS_FILE}2')
    os.system(f'mv {TEMP_FREQS_FILE}2 {TEMP_FREQS_FILE}')

    # read data
    #freqs level
    f = pd.read_csv(TEMP_FREQS_FILE, low_memory=False)
    f['ref_position'] = f.ref_position.astype(float)
    f = f.sort_values(['File', 'ref_position'])
    f = estimate_insertion_freq(f, extra_columns=['File'])
    f['full_mutation'] = f.ref_base + f.ref_position.astype(str) + f.base
    #f.to_csv(TEMP_FREQS_FILE + '.temp.csv')
    # cons level
    df = pd.read_csv(UNITED_CONS_FILE)

    # run alerts
    # freqs
    freqs_output = []
    for m in freqs_alerts:
        if not m.startswith('-'): # substitutions and deletions
            freqs_output.append(f[(f.full_mutation == m) & (f.frequency > freqs_alerts[m])])
        else: # insertion
            freqs_output.append(f[(f.full_mutation == m) & (f.estimated_freq > freqs_alerts[m])])
    freqs_output = pd.concat(freqs_output)

    # consensuses
    consensus_output = []
    for m in consensus_nuc_alerts:
        o = df[(df.position == float(m[1:-1])) & (df.ref_base.str.upper() == m[0]) & (df.base.str.upper() == m[-1])]
        o['reason'] = consensus_nuc_alerts[m]
        consensus_output.append(o)
    consensus_output = pd.concat(consensus_output)

    # pretty print output:
    output = f'''consensus level events:
{consensus_output.sort_values(['sample', 'position']).to_string()}
freqs level events:
{freqs_output.to_string()}'''
    with open(ALERTS_OUTPUT, 'w') as f:
        f.write(output)
    os.system(f'rm {TEMP_FREQS_FILE}')
