

import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import sys
sys.path.append('/sternadi/home/volume1/shared/SternLab')
sys.path.append('/Volumes/STERNADILABHOME$/volume2/noam/SternLab/')
from blast_utilities import blast_to_df
from freqs_utilities import estimate_insertion_freq, unite_all_freq_files

COLORS = {'T1764.0-':'#F50202', 'A1664.0G':'#F49D09', 'A535.0G':'#5EC0D2', 'T1440.0C':'#F1F87A', 'T1440.0G':'#C4A0E4', 'A1443.0G':'#FE1D1D', 'A1611.0G':'#327CFD', 'C1724.0T':'#8FD95A', 'A1744.0G':'#FBB3DD', 'G1560.0A':'#A7F0A2', 'G1906.0A':'#A3A3A3', 'C3358.0T':'#26451C', 'G3114.0A':'#B37A42', 'A1770.0G':'#7603BA', 'G2310.0A':'#033E86', 'A2626.0G':'#8FD95A', 'C3299.0T':'#211785', 'C1718.0T':'#DFC236', 'T862.0C':'#880E05', 'A2790.0T':'#DF36C6', 'G1736.0A':'#CFFD2F', 'C1549.0T':'#2CA403', 'G531.0A':'#972FFE', 'C1050.0T':'#13B908', 'G1560.0A':'#B5FE84', 'G1688.0T':'#1B0398','A2356.0G':'#830276', 'T170.0A':'#C60DC3', 'A1673.0G':'#E2D492', 'C2859.0T':'#B3FD04', 'G20.0-':'black', 'T21.0C':'grey'}


####### data ###########
df = pd.read_csv('/Volumes/STERNADILABHOME$/volume2/noam/ms2_km/new_pipeline_1e-9/all.freqs.csv')
# only passages:
df = df[~(df.File.str.startswith('GEL')) & ~(df.File.str.startswith('Plasmid'))]
df['Passage'] = df.File.str.replace('-', '').str.split('KM').str[0].str.replace('p', '').astype(int)
df['Replicate'] = df.File.str.replace('-', '').str.split('KM').str[-1].str.replace('p', '').astype(int)

########################

# check coverages
fig, ax = plt.subplots(4,6)
ax = ax.flatten()
i = 1
for f in df[['File']].drop_duplicates().File.tolist():
    df[(df.File == f) & (df.ref_base == df.base) & (df.ref_base != '-')].plot(x='ref_position', y='coverage', ax=ax[i], kind='scatter')
    ax[i].set_title(f)
    #ax[i].set_yscale('log')
    i += 1
fig.set_size_inches(10,10)
fig.subplots_adjust(hspace=0.3)

####### mutations

def create_mutations_graph(df, mutations_list, output_path, title=None):
    """
    This function gets a df of freq files, a list of mutations and a path to save graph to.
    """
    plt.style.use('ggplot')
    fig, axes = plt.subplots(nrows=1, ncols=3)
    axes = axes.flatten()
    for sample, a in zip([1,2,3], axes):
                df_line = df[df.Replicate == sample]
                df_line = df_line.sort_values('Passage')
                for m in mutations_list:
                    df_line_mutation = df_line[df.Full_mutation == m].sort_values('Passage')
                    if m in COLORS:
                        a.plot('Passage', 'frequency', data = df_line_mutation, linestyle='-', marker='.', label = m, color=COLORS[m])
                    else:
                        a.plot('Passage', 'frequency', data = df_line_mutation, linestyle='-', marker='.', label = m)
                a.set_title(sample, fontsize=16)
                a.set_xlabel('Time (Passages)', fontsize=14, color='black')
                a.set_ylim(0,1)
    axes[0].set_ylabel('Mutation Frequency', fontsize=14, color='black')
    fig.set_size_inches(10, 3)
    if title:
        fig.set_title(title)
    #plt.subplots_adjust(wspace=0.3, hspace=0.4)
    plt.legend(bbox_to_anchor=(1.04,0.5), loc="center left", borderaxespad=0, facecolor='white', edgecolor='white')
    plt.savefig(output_path, bbox_inches='tight', dpi=800)
    plt.show()
    return

df['Full_mutation'] = df.ref_base + df.ref_position.astype(str) + df.base
mutations = df[(df.base != df.ref_base) & (df.ref_base != '-') & (df.frequency >= 0.1) & (df.ref_position > 30) & (df.ref_position < 3540)].sort_values('frequency', ascending=False).Full_mutation.unique().tolist()
create_mutations_graph(df, mutations, '/Volumes/STERNADILABHOME$/volume2/noam/ms2_km/new_pipeline_1e-9/mutations_over_time.png')

mutations = df[(df.base != df.ref_base) & (df.ref_base != '-') & (df.frequency >= 0.1) & ((df.ref_position < 30) | (df.ref_position > 3540))].sort_values('frequency', ascending=False).Full_mutation.unique().tolist()
mutations.remove('C3566.0G')
mutations.remove('C3566.0A')
create_mutations_graph(df, mutations, '/Volumes/STERNADILABHOME$/volume2/noam/ms2_km/new_pipeline_1e-9/mutations_over_time_edges.png')




###### mutations around edges on cheater passages

def create_mutations_graph2(df, mutations_list, output_file, title=None):
    """
    This function gets a df of freq files, a list of mutations and a path to save graph to.
    """
    plt.style.use('ggplot')
    fig, axes = plt.subplots(nrows=1, ncols=2)
    axes = axes.flatten()
    for sample, a in zip(['37A', '37B'], axes):
                df_line = df[(df.Replica==sample[-1]) & (df.Degree==int(sample[:2]))]
                df_line = df_line.sort_values('Time')
                for m in mutations_list:
                    df_line_mutation = df_line[(df_line.Pos == float(m[1:-1])) & (df_line.Base == m[-1])].sort_values('Time')
                    if m in COLORS:
                        a.plot('Time', 'Freq', data = df_line_mutation, marker='.', linestyle='-', label = m.replace('.0', '').replace('T', 'U').replace('U1764-', '$\Delta$1764'), color=COLORS[m])
                    else:
                        a.plot('Time', 'Freq', data = df_line_mutation, marker='.', linestyle='-', label = m)
                #a.set_ylim(top=0.8, bottom=0)
                a.set_yticks([0.0,0.2,0.4,0.6,0.8])
                a.set_title('Line ' + sample.replace('37', ''), fontsize=16)
                a.set_xlabel('Time (Passages)', fontsize=14, color='black')
                a.set_ylabel('Mutation Frequency', fontsize=14, color='black')
                #a.minorticks_on()
                #a.grid(which='minor', alpha=0.2)
                #a.grid(which='major', alpha=0.7)
                #a.set_xticks([1,3,5,7,9,11,13,15,17])
                #a.set_xticks(range(0,25,3))
    a.set_ylabel('', fontsize=14)
    fig.set_size_inches(16, 6)
    if title:
        fig.set_title(title)
    #plt.subplots_adjust(wspace=0.3, hspace=0.4)
    plt.legend(bbox_to_anchor=(1.04,0.5), loc="center left", borderaxespad=0, facecolor='white', edgecolor='white')
    plt.savefig(output_file, bbox_inches='tight', dpi=800)
    plt.show()
    return

df = pd.read_csv('/Volumes/STERNADILABHOME$/volume2/noam/passages/201909/all_freqs.csv')

df = df[~(df.Time.isin([11,12,14]))]
mutations = df[(df.Base != df.Ref) & (df.Ref != '-') & (df.Freq >= 0.1) & ((df.Pos < 30) | (df.Pos > 3540))].sort_values('Freq', ascending=False).Full_mutation.unique().tolist()

mutations.remove('G5.0C')
mutations.remove('C3566.0G')

#mutations.remove('G1554.0A')
#mutations.remove('C224.0T')
#mutations.remove('C3299.0T')
#create_mutations_graph2(df, mutations, 'X:/volume2/noam/passages/201909/mutation_over_time.png')
create_mutations_graph2(df, mutations, '/Volumes/STERNADILABHOME$/volume2/noam/passages/201909/mutation_over_time_edges.png')


