################################################################
# This script creates all the graphs appearing in the cheaters paper
# (will rewrite the graphs though, so don't run all of it as is)
################################################################



import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
import os
import math

COLORS = {'T1764.0-':'#F50202', 'A1664.0G':'#F49D09', 'A535.0G':'#5EC0D2', 'T1440.0C':'#F1F87A', 'T1440.0G':'#C4A0E4', 'A1443.0G':'#FE1D1D', 'A1611.0G':'#327CFD', 'C1724.0T':'#8FD95A', 'A1744.0G':'#FBB3DD', 'G1560.0A':'#A7F0A2', 'G1906.0A':'#A3A3A3', 'C3358.0T':'#26451C', 'G3114.0A':'#B37A42', 'A1770.0G':'#7603BA', 'G2310.0A':'#033E86', 'A2626.0G':'#8FD95A', 'C3299.0T':'#211785', 'C1718.0T':'#DFC236', 'T862.0C':'#880E05', 'A2790.0T':'#DF36C6', 'G1736.0A':'#CFFD2F', 'C1549.0T':'#2CA403', 'G531.0A':'#972FFE', 'C1050.0T':'#13B908', 'G1560.0A':'#B5FE84', 'G1688.0T':'#1B0398','A2356.0G':'#830276', 'T170.0A':'#C60DC3', 'A1673.0G':'#E2D492', 'C2859.0T':'#B3FD04'}


n = [5, 17, 18, 19, 20, 21, 22, 23, 3462, 3463, 3524, 3542, 3543, 3544, 3545, 3546, 3547, 3548, 3564, 3566]
        
n = list(range(30)) + list(range(3539, 3570))



def create_mutations_barplot(df, out_file):
    fig, axes = plt.subplots(nrows=1, ncols=5)
    axes = axes.flatten()
    for sample, a in zip([5, 8, 10, 13, 15], axes):
        df2 = df[(df.Full_mutation.isin(['A1664.0G', 'T1764.0-'])) & (df.passage == sample)]
        df2 = df2.sort_values('Full_mutation', ascending=False)
        df2['Freq'] = df2.Freq.round(2)
        wt = (1 - df2.groupby('time').Freq.sum()).reset_index()
        wt['Full_mutation'] = 'other'
        df2 = pd.concat([wt, df2], sort=True)
        df3 = df2.pivot_table(values='Freq', index=['time'], columns='Full_mutation')
        #column_order = ['T1764.0-', 'A1664.0G', 'other']
        column_order = ['A1664.0G', 'T1764.0-','other']
        df3 = df3.reindex(column_order, axis=1)
        df3.plot.bar(ax=a, stacked=True, color=['#F49D09', '#F50202', 'white'], legend=False)
        for p in a.patches[-4:]:
            p.set_hatch('///')
            p.set_edgecolor('#5EC0D2')
        a.set_title('Passage ' + str(sample), fontsize=16)
        a.set_ylim(top=1, bottom=0)
        a.set_xlabel('')
        a.set_ylabel('')
        a.minorticks_on()
        a.tick_params(axis='both', which='major', labelsize=12)
        if sample != 5:
            a.set_yticklabels(['']*6)
        #a.grid(which='minor', alpha=0.2)
        #a.grid(which='major', alpha=0.7)
        for p in a.patches:
            width, height = p.get_width(), p.get_height()
            x, y = p.get_xy() 
            if height >= 0.02:
                a.text(x+width/2.0, y+height/2.0, str(round(height, 2)), horizontalalignment='center', verticalalignment='center', fontsize=12, rotation=30)
        plt.setp(a.get_xticklabels(), rotation=360)
    #plt.legend(bbox_to_anchor=(1.04,0.5), loc="center left", borderaxespad=0, labels=['A1664G', '$\Delta$1764', 'WT/A535G/other'])
    plt.legend(bbox_to_anchor=(-1.95,-0.3), loc='upper center', borderaxespad=0, labels=['A1664G', '$\Delta$1764', 'WT/A535G/other'], ncol=3, fontsize=12)
    fig = plt.gcf()
    fig.text(0.5, -0.02, 'Time Post Infection', ha='center', va='center', fontsize=16)
    fig.text(0.07, 0.5, 'Frequency', ha='center', va='center', rotation='vertical', fontsize=16)
    fig.set_size_inches(12, 4)
    plt.subplots_adjust(wspace=0.15)
    plt.savefig(out_file , bbox_inches='tight', dpi=800)
    plt.show()
    plt.clf()
    
zoom_in = pd.read_csv('X:/volume2/noam/zoom_in_passages/zoom_in_freqs.csv')
create_mutations_barplot(zoom_in, 'X:/volume2/noam/zoom_in_passages/zoom_in_barplot.png')



########################3
# plaques
f = pd.read_csv('Z:/volume1/noam/ms2_data/Moran_Yaara_Oded-119067150/pipeline_runs/StokeT30/StokeT30.freqs', '\t')
f = estimate_insertion_freq(f, extra_columns=[])
f[((f.Base != f.Ref) & (f.Ref != '-') & (f.Freq > 0.1)) | ((f.Ref =='-') & (f.Base != '-') & (f.estimated_freq > 0.1))].to_clipboard()


f = pd.read_csv('Z:/volume1/noam/ms2_data/190415_M04473_0037_000000000-G3PG4/12_pipeline/12.freqs', '\t')
f = estimate_insertion_freq(f, extra_columns=[])
f[((f.Base != f.Ref) & (f.Ref != '-') & (f.Freq > 0.1)) | ((f.Ref =='-') & (f.Base != '-') & (f.estimated_freq > 0.1))].to_clipboard()


###### lysis plate reader

def lysis(input_plate_reader_excel, input_mapping_excel, output_directory) :      
    if not os.path.isdir(output_directory + '/triplicates'):
        os.mkdir(output_directory + '/triplicates')
    if not os.path.isdir(output_directory + '/samples'):
        os.mkdir(output_directory + '/samples')
    lysis = pd.read_excel(input_plate_reader_excel, header=4)
    lysis = lysis[~(lysis.apply(lambda row: row.astype(str).str.contains('\?').any(), axis=1))]
    lysis['kinetic_read'] = pd.to_datetime(lysis['Kinetic read'], format='%H:%M:%S').dt.time
    mapping = pd.read_excel(input_mapping_excel)
    mapping['plate_reader_columns'] = mapping.plate_reader_columns.str.split(',')
    for s in mapping.sample_name.tolist():
        sample_columns = mapping[mapping.sample_name == s].iloc[0]['plate_reader_columns']
        data = lysis[['kinetic_read'] + sample_columns]
        #fig, axes = plt.subplots(nrows=1, ncols=1)
        #data.plot(x='kinetic_read', ax=axes, xticks=range(0, 3600*6 + 1, 3600), legend=False)
        #axes.set_title(s)
        #axes.xaxis.set_tick_params(rotation=45)
        #fig.legend(loc='center right')
        #fig.savefig(output_directory + '/triplicates/' + s + '.png', dpi=800, bbox_inches='tight')
        lysis[s] = data[sample_columns].mean(axis=1)
    for gn in mapping.graph_name.drop_duplicates().tolist():
        fig, axes = plt.subplots(nrows=1, ncols=1)
        lysis2 = lysis[['kinetic_read'] + mapping[mapping.graph_name == gn].sample_name.tolist()]
        lysis2.plot(x='kinetic_read', ax=axes, xticks=range(0, 3600*6 + 1, 3600), legend=False, cmap='viridis')
        #axes.set_title('The Timing of Host Cell Lysis')
        fig.legend(loc='center left', bbox_to_anchor=(1.0, 0.7))
        axes.set_xlabel('Time Post Infection (hrs)', fontsize=14)
        axes.set_ylabel('O.D.600', fontsize=14)
        axes.xaxis.set_tick_params(rotation=45)
        fig.set_size_inches(3, 2)
        fig.savefig(output_directory + '/samples/lysis_' + gn + '.png', dpi=800, bbox_inches='tight')

lysis('X:/volume2/noam/lysis/Lysis 1.4.19 (0001).xls', 'X:/volume2/noam/lysis/map_c3000_1.4.xlsx', 'X:/volume2/noam/lysis/c3000/')
lysis('X:/volume2/noam/lysis/Helper13.5.19 (0001).xls', 'X:/volume2/noam/lysis/map_helper.xlsx', 'X:/volume2/noam/lysis/helper/')


def c3000_helper(input_plate_reader_excel_c3000, c3000_columns, input_plate_reader_excel_helper, helper_columns, output_file) :
    lysis_c3000 = pd.read_excel(input_plate_reader_excel_c3000, header=4)
    lysis_c3000 = lysis_c3000[~(lysis_c3000.apply(lambda row: row.astype(str).str.contains('\?').any(), axis=1))]
    lysis_c3000['kinetic_read'] = pd.to_datetime(lysis_c3000['Kinetic read'], format='%H:%M:%S').dt.time
    lysis_c3000['c3000'] = lysis_c3000[c3000_columns].mean(axis=1)
    
    lysis_helper = pd.read_excel(input_plate_reader_excel_helper, header=4)
    lysis_helper = lysis_helper[~(lysis_helper.apply(lambda row: row.astype(str).str.contains('\?').any(), axis=1))]
    lysis_helper['kinetic_read'] = pd.to_datetime(lysis_helper['Kinetic read'], format='%H:%M:%S').dt.time
    lysis_helper['kinetic_read'] = pd.to_datetime(lysis_helper['Kinetic read'], format='%H:%M:%S').dt.time
    lysis_helper['helper'] = lysis_helper[helper_columns].mean(axis=1)
    
    fig, axes = plt.subplots(nrows=1, ncols=1)
    lysis_c3000.plot(x='kinetic_read', y='c3000', kind='line', ax=axes, xticks=range(0, 3600*6 + 1, 3600), legend=False, label='C3000', color='#940786')
    lysis_helper.plot(x='kinetic_read', y='helper', kind='line', ax=axes, xticks=range(0, 3600*6 + 1, 3600), legend=False, label='XL1 + pBAD + Replicase', color='#9CCC65')
    fig.legend(loc='center left', bbox_to_anchor=(1.0, 0.7))
    axes.set_xlabel('Time Post Infection (hrs)', fontsize=14)
    axes.set_ylabel('O.D.600', fontsize=14)
    axes.xaxis.set_tick_params(rotation=45)
    fig.set_size_inches(3, 2)
    fig.savefig(output_file, dpi=800, bbox_inches='tight')
    return

c3000_helper('X:/volume2/noam/lysis/Lysis 1.4.19 (0001).xls', ['B1', 'B2', 'B3'], 'X:/volume2/noam/lysis/Helper13.5.19 (0001).xls', ['B1', 'B2', 'B3'], 'X:/volume2/noam/lysis/c3000_vs_helper.png')    


  

# lysis main text

def lysis_main_text(input_plate_reader_excel, input_mapping_excel, output) :    
    plt.style.use('default')
    #if not os.path.isdir(output_directory + '/triplicates'):
    #    os.mkdir(output_directory + '/triplicates')
    #if not os.path.isdir(output_directory + '/samples'):
    #    os.mkdir(output_directory + '/samples')
    lysis = pd.read_excel(input_plate_reader_excel, header=4)
    lysis = lysis[~(lysis.apply(lambda row: row.astype(str).str.contains('\?').any(), axis=1))]
    lysis['kinetic_read'] = pd.to_datetime(lysis['Kinetic read'], format='%H:%M:%S').dt.time
    mapping = pd.read_excel(input_mapping_excel)
    mapping['plate_reader_columns'] = mapping.plate_reader_columns.str.split(',')
    for s in mapping.sample_name.tolist():
        sample_columns = mapping[mapping.sample_name == s].iloc[0]['plate_reader_columns']
        data = lysis[['kinetic_read'] + sample_columns]
        #fig, axes = plt.subplots(nrows=1, ncols=1)
        #data.plot(x='kinetic_read', ax=axes, xticks=range(0, 3600*6 + 1, 3600), legend=False)
        #axes.set_title(s)
        #axes.xaxis.set_tick_params(rotation=45)
        #fig.legend(loc='center right')
        #fig.savefig(output_directory + '/triplicates/' + s + '.png', dpi=800, bbox_inches='tight')
        lysis[s] = data[sample_columns].mean(axis=1)
    fig, axes = plt.subplots(nrows=2, ncols=1)
    for gn, a in zip(mapping.graph_name.drop_duplicates().tolist(), axes):
        lysis2 = lysis[['kinetic_read'] + mapping[mapping.graph_name == gn].sample_name.tolist()]
        lysis2.plot(x='kinetic_read', ax=a, xticks=range(0, 3600*6 + 1, 3600), legend=False, cmap='Dark2')        
        a.set_ylabel('', fontsize=14)
        a.xaxis.set_tick_params(rotation=45)
        #a.set_xticklabels(labels=['0','1','2','3','4','5','6'])
        a.set_title(gn.replace('37', 'Line '))
        a.set_xlabel('', fontsize=14)
    fig.text(-0.05, 0.5, 'O.D.600', va='center', rotation='vertical', fontsize = 14)
    a.set_xlabel('Time Post Infection (hrs)', fontsize=14)
    plt.subplots_adjust(hspace=0.7)
    
    handles, labels = a.get_legend_handles_labels()
    labels = [x.split('-')[0] for x in labels]
    a.legend(labels=labels, loc='center left', bbox_to_anchor=(1.05, 1.2))
    
    #fig.legend(loc='center left', bbox_to_anchor=(1.05, 0.5))
    fig.set_size_inches(3, 4)
    fig.savefig(output, dpi=800, bbox_inches='tight')

lysis_main_text('X:/volume2/noam/lysis/Lysis 1.4.19 (0001).xls', 'X:/volume2/noam/lysis/map_c3000_1.4.xlsx', 'X:/volume2/noam/lysis/c3000/lysis_main_text.png')



####### mutation trajectories
# original passages
def create_mutations_graph(df, mutations_list, output_file, title=None):
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
                        a.plot('Time', 'Freq', data = df_line_mutation, linestyle='-', marker='.', label = m.replace('.0', '').replace('T', 'U').replace('U1764-', '$\Delta$1764'), color=COLORS[m])
                    else:
                        a.plot('Time', 'Freq', data = df_line_mutation, linestyle='-', marker='.', label = m)
                a.set_ylim(top=0.8, bottom=0)
                a.set_yticks([0.0,0.2,0.4,0.6,0.8])
                a.set_title('Line ' + sample.replace('37', ''), fontsize=16)
                a.set_xlabel('Time (Passages)', fontsize=14, color='black')
                a.set_ylabel('Mutation Frequency', fontsize=14, color='black')
                #a.minorticks_on()
                #a.grid(which='minor', alpha=0.2)
                #a.grid(which='major', alpha=0.7)
                #a.set_xticks([1,3,5,7,9,11,13,15])
    a.set_ylabel('', fontsize=14)
    fig.set_size_inches(8, 3)
    if title:
        fig.set_title(title)
    #plt.subplots_adjust(wspace=0.3, hspace=0.4)
    plt.legend(bbox_to_anchor=(1.04,0.5), loc="center left", borderaxespad=0, facecolor='white', edgecolor='white')
    plt.savefig(output_file, bbox_inches='tight', dpi=800)
    plt.show()
    return

joined = pd.read_csv('X:/volume2/noam/passages/201908_w_2019_passages/old_passages/all_freqs.csv')
joined['Full_mutation'] = joined.Ref + joined.Pos.astype(str) + joined.Base
#
joined = joined[(joined.Time <= 15) & (joined.Degree == 37)]
joined = joined[~joined.Time.isin([11, 12, 14])]

mutations = joined[(joined.Base != joined.Ref) & (joined.Ref != '-') & ~(joined.Pos.isin(n)) & (joined.Freq >= 0.1)].sort_values('Freq', ascending=False).Full_mutation.unique().tolist()
mutations.remove('C3299.0T')
mutations.remove('G1554.0A')
mutations.remove('C224.0T')
create_mutations_graph(joined, mutations, 'X:/volume2/noam/passages/201908_w_2019_passages/old_passages/mutations_over_time_p15_no11,12,14.png')

### 41C 15 passages for review

def create_mutations_graph(df, mutations_list, output_file, title=None):
    """
    This function gets a df of freq files, a list of mutations and a path to save graph to.
    """
    plt.style.use('ggplot')
    fig, axes = plt.subplots(nrows=1, ncols=2)
    axes = axes.flatten()
    for sample, a in zip(['41A', '41B'], axes):
                df_line = df[(df.Replica==sample[-1]) & (df.Degree==int(sample[:2]))]
                df_line = df_line.sort_values('Time')
                for m in mutations_list:
                    df_line_mutation = df_line[(df_line.Pos == float(m[1:-1])) & (df_line.Base == m[-1])].sort_values('Time')
                    if m in COLORS:
                        a.plot('Time', 'Freq', data = df_line_mutation, linestyle='-', marker='.', label = m.replace('.0', '').replace('T', 'U').replace('U1764-', '$\Delta$1764'), color=COLORS[m])
                    else:
                        a.plot('Time', 'Freq', data = df_line_mutation, linestyle='-', marker='.', label = m)
                a.set_ylim(top=0.8, bottom=0)
                a.set_yticks([0.0,0.2,0.4,0.6,0.8])
                a.set_title('Line ' + sample.replace('37', ''), fontsize=16)
                a.set_xlabel('Time (Passages)', fontsize=14, color='black')
                a.set_ylabel('Mutation Frequency', fontsize=14, color='black')
                #a.minorticks_on()
                #a.grid(which='minor', alpha=0.2)
                #a.grid(which='major', alpha=0.7)
                #a.set_xticks([1,3,5,7,9,11,13,15])
    a.set_ylabel('', fontsize=14)
    fig.set_size_inches(8, 3)
    if title:
        fig.set_title(title)
    #plt.subplots_adjust(wspace=0.3, hspace=0.4)
    plt.legend(bbox_to_anchor=(1.04,0.5), loc="center left", borderaxespad=0, facecolor='white', edgecolor='white')
    plt.savefig(output_file, bbox_inches='tight', dpi=800)
    plt.show()
    return

joined = pd.read_csv('X:/volume2/noam/passages/201908_w_2019_passages/old_passages/all_freqs.csv')
joined['Full_mutation'] = joined.Ref + joined.Pos.astype(str) + joined.Base
#
joined = joined[(joined.Time <= 15) & (joined.Degree == 41)]
joined = joined[~joined.Time.isin([11, 12, 14])]

mutations = joined[(joined.Base != joined.Ref) & (joined.Ref != '-') & ~(joined.Pos.isin(n)) & (joined.Freq >= 0.1)].sort_values('Freq', ascending=False).Full_mutation.unique().tolist()
mutations.remove('C3299.0T')
mutations.remove('G1554.0A')
mutations.remove('C224.0T')
create_mutations_graph(joined, mutations, 'X:/volume2/noam/passages/201908_w_2019_passages/old_passages/mutations_over_time_41C_p15_no11,12,14.png')



# new passages
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
                a.set_ylim(top=0.8, bottom=0)
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
    fig.set_size_inches(8, 3)
    if title:
        fig.set_title(title)
    #plt.subplots_adjust(wspace=0.3, hspace=0.4)
    plt.legend(bbox_to_anchor=(1.04,0.5), loc="center left", borderaxespad=0, facecolor='white', edgecolor='white')
    plt.savefig(output_file, bbox_inches='tight', dpi=800)
    plt.show()
    return

df = pd.read_csv('X:/volume2/noam/passages/201909/all_freqs.csv')
df = df[~(df.Time.isin([11,12,14]))]

mutations = df[(df.Base != df.Ref) & (df.Ref != '-') & ~(df.Pos.isin(n)) & (df.Freq >= 0.1)].sort_values('Freq', ascending=False).Full_mutation.unique().tolist()
mutations.remove('G1554.0A')
mutations.remove('C224.0T')
mutations.remove('C3299.0T')
create_mutations_graph2(df, mutations, 'X:/volume2/noam/passages/201909/mutation_over_time.png')

# only original mutations
joined = pd.read_csv('X:/volume2/noam/passages/201908_w_2019_passages/old_passages/all_freqs.csv')
joined['Full_mutation'] = joined.Ref + joined.Pos.astype(str) + joined.Base
joined = joined[(joined.Time <= 15) & (joined.Degree == 37)]
joined = joined[~joined.Time.isin([11, 12, 14])]
mutations = joined[(joined.Base != joined.Ref) & (joined.Ref != '-') & ~(joined.Pos.isin(n)) & (joined.Freq >= 0.1)].sort_values('Freq', ascending=False).Full_mutation.unique().tolist()
mutations.remove('C3299.0T')
mutations.remove('G1554.0A')
mutations.remove('C224.0T')
create_mutations_graph2(df, mutations, 'X:/volume2/noam/passages/201909/mutation_over_time_simplified.png')

##########
# diverging_mois
def create_mutations_graph3(df, mutations_list, output_file, moi, title=None):
    """
    This function gets a df of freq files, a list of mutations and a path to save graph to.
    """
    plt.style.use('ggplot')
    fig, axes = plt.subplots(nrows=1, ncols=2)
    axes = axes.flatten()
    for sample, a in zip(['A', 'B'], axes):
                df_line = df[(df.Replica==sample)]
                df_line = df_line.sort_values('Time')
                for m in mutations_list:
                    df_line_mutation = df_line[(df_line.Pos == float(m[1:-1])) & (df_line.Base == m[-1])].sort_values('Time')
                    if m in COLORS:
                        a.plot('Time', 'Freq', data = df_line_mutation, marker='.', linestyle='-', label = m.replace('.0', '').replace('T', 'U').replace('U1764-', '$\Delta$1764'), color=COLORS[m])
                    else:
                        a.plot('Time', 'Freq', data = df_line_mutation, marker='.', linestyle='-', label = m)
                a.set_ylim(top=0.6, bottom=0)
                a.set_yticks([0.0,0.2,0.4,0.6])
                a.set_title('Line ' + sample.replace('37', ''), fontsize=16)
                a.set_xlabel('Time (Passages)', fontsize=14, color='black')
                a.set_ylabel('Mutation Frequency', fontsize=14, color='black')
                #a.minorticks_on()
                #a.grid(which='minor', alpha=0.2)
                #a.grid(which='major', alpha=0.7)
                a.set_xticks([5,10,15,16])
                a.xaxis.set_tick_params(rotation=45)
                a.axvspan(15, 17, facecolor='#B5DF8B', alpha=0.5)
                a.set_xlim(0,17)
                
    
    patch2 = mpatches.Patch(color='#9BBD78', label='MOI ' + moi, alpha=0.7)
    patch1 = mpatches.Patch(color='#ECEBEB', label='MOI 1')
    lgd2 = axes[0].legend(bbox_to_anchor=(1.1, -0.48), handles=[patch1, patch2], loc='lower center', borderaxespad=0, facecolor='white', edgecolor='white', ncol=2)
    
    a.set_ylabel('', fontsize=14)
    fig.set_size_inches(8, 2.4)
    if title:
        fig.set_title(title)
    #plt.subplots_adjust(wspace=0.3, hspace=0.4)
    lgd = plt.legend(bbox_to_anchor=(1.04,0.5), loc="center left", borderaxespad=0, facecolor='white', edgecolor='white')
    plt.savefig(output_file, bbox_inches='tight', dpi=800)#, bbox_extra_artists=[lgd,lgd2])
    plt.show()
    return

joined = pd.read_csv('X:/volume2/noam/passages/201908_w_2019_passages/new_2019/all_freqs.csv')
joined[(joined.Degree == 37)]
joined = joined[~(joined.Time.isin([11, 12, 14])) & (joined.Time <= 15)]
joined['Full_mutation'] = joined.Ref + joined.Pos.astype(str) + joined.Base

moi_01_16_37A = pd.read_csv('X:/volume2/noam/passages/201908_w_2019_passages/freqs/P16-01-37A-2019.freqs', '\t')
moi_01_16_37B = pd.read_csv('X:/volume2/noam/passages/201908_w_2019_passages/freqs/P16-01-37B-2019.freqs', '\t')
moi_01_16_37A['Replica'] = 'A'
moi_01_16_37B['Replica'] = 'B'
moi_01_16_37A['Time'] = 16
moi_01_16_37B['Time'] = 16
moi_01_16_37A['Degree'] = 37
moi_01_16_37B['Degree'] = 37

moi_01_16 = pd.concat([moi_01_16_37A, moi_01_16_37B, joined])
moi_01_16['Full_mutation'] = moi_01_16.Ref + moi_01_16.Pos.astype(str) + moi_01_16.Base

mutations = moi_01_16[(moi_01_16.Base != moi_01_16.Ref) & (moi_01_16.Ref != '-') & ~(moi_01_16.Pos.isin(n)) & (moi_01_16.Freq >= 0.1)].sort_values('Freq', ascending=False).Full_mutation.unique().tolist()
mutations.remove('G1554.0A')
mutations.remove('C224.0T')
mutations.remove('C3299.0T')
create_mutations_graph3(moi_01_16, mutations, 'X:/volume2/noam/passages/201908_w_2019_passages/diverging_mois/p16_moi_0.1_no11,12,14.png', '0.1')



moi_001_16_37A = pd.read_csv('X:/volume2/noam/passages/201908_w_2019_passages/freqs/P16-001-37A-2019.freqs', '\t')
moi_001_16_37B = pd.read_csv('X:/volume2/noam/passages/201908_w_2019_passages/freqs/P16-001-37B-2019.freqs', '\t')
moi_001_16_37A['Replica'] = 'A'
moi_001_16_37B['Replica'] = 'B'
moi_001_16_37A['Time'] = 16
moi_001_16_37B['Time'] = 16
moi_001_16_37A['Degree'] = 37
moi_001_16_37B['Degree'] = 37

moi_001_16 = pd.concat([moi_001_16_37A, moi_001_16_37B, joined])
moi_001_16['Full_mutation'] = moi_001_16.Ref + moi_001_16.Pos.astype(str) + moi_001_16.Base

mutations = moi_001_16[(moi_001_16.Base != moi_001_16.Ref) & (moi_001_16.Ref != '-') & ~(moi_001_16.Pos.isin(n)) & (moi_001_16.Freq >= 0.1)].sort_values('Freq', ascending=False).Full_mutation.unique().tolist()
mutations.remove('G1554.0A')
mutations.remove('C224.0T')
mutations.remove('C3299.0T')
create_mutations_graph3(moi_001_16, mutations, 'X:/volume2/noam/passages/201908_w_2019_passages/diverging_mois/p16_moi_0.01_no11,12,14.png', '0.01')

############# 20 plaques
# get original mutations
joined = pd.read_csv('X:/volume2/noam/passages/201908_w_2019_passages/old_passages/all_freqs.csv')
joined['Full_mutation'] = joined.Ref + joined.Pos.astype(int).astype(str) + joined.Base
joined = joined[(joined.Time <= 15) & (joined.Degree == 37)]
joined = joined[~joined.Time.isin([11, 12, 14])]
mutations = joined[(joined.Base != joined.Ref) & (joined.Ref != '-') & ~(joined.Pos.isin(n)) & (joined.Freq >= 0.1)].sort_values('Freq', ascending=False)[['Full_mutation']].drop_duplicates()
mutations = mutations[~(mutations.Full_mutation.isin(['C3299T', 'G1554A', 'C224T']))]

df = pd.read_csv('X:/volume2/noam/plaques_40/all_freqs.csv')
df['Full_mutation'] = df.Ref + df.Pos.astype(int).astype(str) + df.Base
df['Line'] = df.Sample.str.split('_').str[1]
df = df[df.Line.str.contains('37')]
df = pd.merge(mutations, df, on='Full_mutation')

df = df[(df.Freq > 0.8)].groupby(['Line', 'Full_mutation']).Sample.count().reset_index()
df = df.rename(columns={'Sample':'Plaque_freq'})
df['Plaque_freq'] = df.Plaque_freq / 10

mutations37A = mutations.copy()
mutations37A['Line'] = '37A'
mutations37B = mutations.copy()
mutations37B['Line'] = '37B'
df = pd.merge(df, pd.concat([mutations37A, mutations37B]), how='right', on=['Line', 'Full_mutation'])
df = df.fillna(0)
df['Line'] = 'Line ' + df['Line'].str.replace('37', '')

df['Full_mutation'] = df.Full_mutation.str.replace('T', 'U').str.replace('U1764-', '$\Delta$1764')

plt.style.use('default')
fig, axes = plt.subplots(nrows=1, ncols=2)
for sample, a in zip(['Line A', 'Line B'], axes):
    sns.stripplot(x='Line', y='Plaque_freq', hue='Full_mutation', data=df[(df.Line==sample) & (df.Plaque_freq != 0)], jitter=0, size=7, palette={i.replace('.0', '').replace('T', 'U').replace('U1764-', '$\Delta$1764'):COLORS[i] for i in COLORS}, ax=a)
    sns.stripplot(x='Line', y='Plaque_freq', hue='Full_mutation', data=df[(df.Line==sample) & (df.Plaque_freq == 0)], jitter=0.2, size=7, palette={i.replace('.0', '').replace('T', 'U').replace('U1764-', '$\Delta$1764'):COLORS[i] for i in COLORS}, ax=a, alpha=0.3)
    a.set_ylabel('Plaque Frequency', fontsize=14)
    a.set_xlabel('')
    a.minorticks_on()
    a.grid(which='major', alpha=0.7, axis='y')
    a.get_legend().remove()
    a.set_ylim(-0.05,0.35)
    a.tick_params(axis='x', which='major', labelsize=12)
    #a.set_title(sample.replace('_', ' '))
a.legend(bbox_to_anchor=(1.04,0.5), loc="center left", borderaxespad=0)
a.set_ylabel('')
fig.set_size_inches(3.5, 2.5)
fig.subplots_adjust(wspace=0.4)
fig.savefig('X:/volume2/noam/plaques_40/37_plaques_new3.png', bbox_inches='tight', dpi=800)




################# coverage plots
joined = pd.read_csv('X:/volume2/noam/passages/201908_w_2019_passages/old_passages/all_freqs.csv')
joined = joined[(joined.Time <= 15) & (joined.Degree == 37)]
joined = joined[~joined.Time.isin([11, 12, 14])]
joined['Sample'] = joined.Sample.str.replace('P', 'p')
joined = joined.sort_values(['Time', 'Replica', 'Pos'])
plasmid = pd.read_csv('X:/volume2/noam/passages/201907/freqs/plasmid_for_p5,13,15,18,20.freqs', '\t')
plasmid['Sample'] = 'plasmid'
joined = pd.concat([joined, plasmid])

fig, ax = plt.subplots(nrows=1, ncols=1)
colors = list(COLORS.values())
colors = colors[:-5]
for name, group in joined.groupby('Sample', sort=False):
    group[(group.Base == group.Ref) & (group.Ref != '-')].plot(x='Pos', y='Read_count', ax=ax, label=name.replace('-37', ''), color = colors[-1])
    colors.pop(-1)
ax.set_yscale('log')
ax.set_ylabel('Read count')
ax.set_xlabel('Position in the genome')
ax.legend(bbox_to_anchor=(1.04,0.5), loc="center left", borderaxespad=0)
fig.set_size_inches(4,3)
fig.savefig('X:/volume2/noam/passages/201908_w_2019_passages/old_passages/coverage_p15_no11,12,14.png', bbox_inches='tight', dpi=800)


df = pd.read_csv('X:/volume2/noam/passages/201909/all_freqs.csv')
df = df[~(df.Time.isin([11,12,14]))]
df = df[df.Time >= 15]

fig, ax = plt.subplots(nrows=1, ncols=1)
colors = list(COLORS.values())
colors = colors[:-5]
for name, group in df.groupby('Sample', sort=False):
    group[(group.Base == group.Ref) & (group.Ref != '-')].plot(x='Pos', y='Read_count', ax=ax, label=name.replace('-37', '').replace('P', 'p'), color = colors[0])
    colors.pop(0)
ax.set_yscale('log')
ax.set_ylabel('Read count')
ax.set_xlabel('Position in the genome')
ax.legend(bbox_to_anchor=(1.04,0.5), loc="center left", borderaxespad=0)
fig.set_size_inches(4,3)
fig.savefig('X:/volume2/noam/passages/201909/ccoverage_plots_2019.png', bbox_inches='tight', dpi=800)


########################## inflection point


a = pd.read_excel('X:/volume2/noam/inflection_point/inflection_WT.xlsx')
X = a.minutes_after_infection.values
Y = np.log10(a.PFU_ml.values)


def jmp_4pl(x, A, B, C, D):
    return (D + ((A-D) / (1 + math.exp((x - C)*B))))

A = 6.85 # minimum asymptote
B = 0.044 # hill's slope, steepness
C = 66.79 # inflection point
D = 11.27 # maximum asymptote

fig, ax = plt.subplots(nrows=1, ncols=1)
xs = np.linspace(X.min() + 0.00000001, X.max(), 10000)
ax.scatter(X, Y, color='#830276')
ax.plot(xs, [jmp_4pl(x, A, B, C, D) for x in xs], color='#830276')
ax.set_xlabel('Time after infection (min)')
ax.set_ylabel('Log10(PFU/ml)')
fig.set_size_inches(4,3)
fig.savefig('X:/volume2/noam/inflection_point/inflection_WT.png', bbox_inches='tight', dpi=800)


####### Real time graph ######
fig, ax = plt.subplots(nrows=1, ncols=1)
data = pd.read_excel('X:/volume2/noam/ms2_real_time/qPCR.graph.xlsx')
data = data[data.time != 0]
data.plot.bar(x='time', y='norm data 7', ax=ax, legend=False, color='purple')
fig.set_size_inches(4,3)
ax.set_yscale('log')
ax.set_xlabel('Time (min)')
ax.set_ylabel('Relative Normalized MS2 Genomes')
ax.xaxis.set_tick_params(rotation=0)
fig.savefig('X:/volume2/noam/ms2_real_time/qPCR.graph.png', bbox_inches='tight', dpi=800)


########## p15 old and p15 new the same ######
COLORS = {'$\\Delta$1764':'#F50202', 'A1664G':'#F49D09', 'A535G':'#5EC0D2', 'U1440C':'#F1F87A', 'U1440G':'#C4A0E4','A1611G':'#327CFD', 'A1744G':'#FBB3DD', 'G1906A':'#A3A3A3', 'G3114A':'#B37A42', 'G531A':'#972FFE', 'C1050U':'#13B908'}


new_A = pd.read_csv('X:/volume2/noam/passages/201908_w_2019_passages/freqs/P15-37A-2019.freqs', '\t')
new_B = pd.read_csv('X:/volume2/noam/passages/201908_w_2019_passages/freqs/P15-37B-2019.freqs', '\t')
old_A = pd.read_csv('X:/volume2/noam/passages/freqs/p15-37A.freqs', '\t')
old_B = pd.read_csv('X:/volume2/noam/passages/freqs/p15-37B.freqs', '\t')
new_A['Sample'], new_A['Type'] = 'p15_A', 'new'
new_B['Sample'], new_B['Type'] = 'p15_B', 'new'
old_A['Sample'], old_A['Type'] = 'p15_A', 'original'
old_B['Sample'], old_B['Type'] = 'p15_B', 'original'
df = pd.concat([new_A, new_B, old_A, old_B])
df['Full_mutation'] = df.Ref + df.Pos.astype(int).astype(str) + df.Base
mutations = df[(df.Base != df.Ref) & (df.Ref != '-') & ~(df.Pos.isin(n)) & (df.Freq >= 0.1)].sort_values('Freq', ascending=False).Full_mutation.drop_duplicates().tolist()
mutations.remove('C3299T')
df = df[df.Full_mutation.isin(mutations)]
df['Full_mutation'] = df.Full_mutation.str.replace('T', 'U').str.replace('U1764-', '$\Delta$1764')
fig, axes = plt.subplots(nrows=1, ncols=2)
for sample, a in zip(['p15_A', 'p15_B'], axes):
    sns.stripplot(x='Type', y='Freq', hue='Full_mutation', data=df[df.Sample==sample].sort_values('Freq', ascending=False), jitter=0.2, size=10, palette={i.replace('.0', ''):COLORS[i] for i in COLORS}, ax=a)
    a.set_ylabel('Frequency')
    a.set_xlabel('')
    a.minorticks_on()
    a.grid(which='major', alpha=0.7, axis='y')
    a.get_legend().remove()
    a.set_ylim(-0.02,0.6)
    a.set_title(sample.replace('_', ' '))
a.legend(bbox_to_anchor=(1.04,0.5), loc="center left", borderaxespad=0)
fig.set_size_inches(5, 3)
fig.subplots_adjust(wspace=0.4)
fig.savefig('X:/volume2/noam/passages/201908_w_2019_passages/new_2019/compare_new_original_p15s.png', bbox_inches='tight', dpi=800)


############## GFP 1664 tests ########
results = pd.read_excel('X:/volume2/noam/gfp_1664/results.xlsx')
results['normalized_fluorescence'] = results.fluorescence - results.baseline
#results['normalized_fluorescence'] = results.fluorescence / results.baseline
import seaborn as sns
fig, ax = plt.subplots(nrows=1, ncols=1)
fig.set_size_inches(4,3)
sns.barplot(x='sample', y='normalized_fluorescence', data=results, ax=ax, palette=['white', '#F49D09'])
ax.set_ylabel('GFP (A.U.)')
ax.set_xlabel('RBS')
ax.patches[0].set_hatch('///')
ax.patches[0].set_edgecolor('#5EC0D2')
fig.savefig('X:/volume2/noam/gfp_1664/gfp_results_minus.png', bbox_inches='tight', dpi=800)

# other version
results = pd.read_excel('X:/volume2/noam/gfp_1664/results2.xlsx')
import seaborn as sns
fig, ax = plt.subplots(nrows=1, ncols=1)
fig.set_size_inches(5,2.5)
sns.barplot(x='sample', y='fluorescence', data=results, ax=ax, palette=['white', '#F49D09', 'white'])
ax.set_ylabel('GFP (A.U.)')
ax.set_xlabel('Sample')
ax.patches[0].set_hatch('///')
ax.patches[0].set_edgecolor('#5EC0D2')
ax.patches[2].set_hatch('')
ax.patches[2].set_edgecolor('black')
fig.savefig('X:/volume2/noam/gfp_1664/gfp_results.png', bbox_inches='tight', dpi=800)



##################### abc smc ############

def get_best_parameters(astroabc_output, n):
    data = pd.read_csv(astroabc_output, '\t', index_col=False)
    data = data.tail(n)
    data['wt_with_wt_p'] = 1
    data['del_with_del_p'] = 0
    data = data.rename(columns={
     'param#0 ':'wt_with_del_p',
     ' param#1 ':'wt_with_syn_p',
     ' param#2 ':'del_with_wt_p',
     ' param#3 ':'del_with_syn_p',
     ' param#4 ':'syn_with_wt_p',
     ' param#5 ':'syn_with_del_p',
     ' param#6 ':'syn_with_syn_p',
     ' param#7 ':'triple_wt',
     ' param#8 ':'triple_del',
     ' param#9 ':'triple_syn',
     ' param#10 ':'starting_n_syn'})
    return data

# multi panel graph
#
#def prior_posterior_graph(data, outpath):
#    data['del_with_wt_syn_with_wt_ratio'] = data.del_with_wt_p / data.syn_with_wt_p
#    data['wt_with_del_wt_with_syn_ratio'] = data.wt_with_del_p / data.wt_with_syn_p
#    #data['del_with_syn_syn_with_del_dratio'] = np.clip(((data.del_with_syn_p + 0.05) / (data.syn_with_del_p + 0.05)), 0, 4.000001)
#    data['del_with_syn_syn_with_del_dratio'] = (data.del_with_syn_p + 0.05) / (data.syn_with_del_p + 0.05)
#    data.loc[data['del_with_syn_syn_with_del_dratio'] > 3, 'del_with_syn_syn_with_del_dratio'] = 4.5
#    
#    fig, axes = plt.subplots(nrows=5, ncols=3)
#    fig.set_size_inches(6,7)
#    fig.tight_layout()
#    axes = axes.flatten()
#    fig.subplots_adjust(hspace=0.5)
#    fig.text(0.5, -0.02, 'Relationship', ha='center', fontsize = 14)
#    fig.text(0.5, 0.38, 'Payoff', ha='center', fontsize = 14)
#    fig.text(-0.05, 0.72, 'Frequency', va='center', rotation='vertical', fontsize = 14)
#    fig.text(-0.05, 0.12, 'Frequency', va='center', rotation='vertical', fontsize = 14)
#    for i, column in enumerate(['wt_with_wt_p', 'wt_with_del_p', 'wt_with_syn_p', 'del_with_wt_p', 'del_with_del_p', 'del_with_syn_p', 'syn_with_wt_p', 'syn_with_del_p', 'syn_with_syn_p', None, None, None, 'del_with_wt_syn_with_wt_ratio', 'wt_with_del_wt_with_syn_ratio', 'del_with_syn_syn_with_del_dratio']):
#        if column:
#            if column in ['wt_with_del_p', 'syn_with_syn_p']:
#                axes[i].hist(np.linspace(0,1,len(data)), alpha=0.2, color='#267BB8', bins=np.linspace(0, 4, 21))
#            elif column in ['del_with_wt_p', 'del_with_syn_p', 'syn_with_wt_p', 'syn_with_del_p', 'wt_with_syn_p']:
#                axes[i].hist(np.linspace(0,4,len(data)), alpha=0.2, color='#267BB8', bins=np.linspace(0, 4, 21))
#            
#            if column not in ['del_with_del_p', 'wt_with_wt_p']:
#                if i < 9:
#                    data[column].plot.hist(ax=axes[i], color='#267BB8', bins=np.linspace(0, 4, 21))
#                    axes[i].set_xticks([0,1,2,3,4])
#                    axes[i].set_xlim(-0.1,4.1)
#                else:
#                    data[column].plot.hist(ax=axes[i], color='green')
#                axes[i].set_xlabel(None)
#                axes[i].set_ylabel(None)
#                
#            else:   
#                axes[i].set_visible(False)
#                            
#            if 'with' in column and i < 9:
#                title = 'W$_\mathrm{' + '|'.join(column.replace('_p', '').split('_with_')).replace('wt', 'WT').replace('del', '\Delta1764').replace('syn', 'A1664G') + '}$'
#                axes[i].set_title(fix_symbols(title), fontsize=12)
#            else:
#                #titles = {'del_with_wt_syn_with_wt_ratio':r'$\frac{\mathrm{W_\mathrm{\Delta1764|WT}}}{\mathrm{W_\mathrm{A1664G|WT}}}$', 'wt_with_del_wt_with_syn_ratio':r'$\frac{\mathrm{W_\mathrm{WT|\Delta1764}}}{\mathrm{W_\mathrm{WT|A1664G}}}$', 'del_with_syn_syn_with_del_dratio':r'$\frac{\mathrm{W_\mathrm{\Delta1764|A1664G}} + 0.05}{\mathrm{W_\mathrm{A1664G|\Delta1764}} + 0.05}$'}
#                titles = {'del_with_wt_syn_with_wt_ratio':r'$\frac{\mathrm{W_\mathrm{\Delta1764|WT}}}{\mathrm{W_\mathrm{A1664G|WT}}}$', 'wt_with_del_wt_with_syn_ratio':r'$\frac{\mathrm{W_\mathrm{WT|\Delta1764}}}{\mathrm{W_\mathrm{WT|A1664G}}}$', 'del_with_syn_syn_with_del_dratio':r'$\frac{\mathrm{W_\mathrm{\Delta1764|A1664G}}}{\mathrm{W_\mathrm{A1664G|\Delta1764}}}$'}
#                axes[i].set_title(fix_symbols(titles[column]), fontsize=14)
#            axes[i].grid(which='major', alpha=0.2, linestyle='-')
#        else:
#            axes[i].set_visible(False)
#    axes[-3].set_xticks([1.0, 1.25, 1.5, 1.75])
#    axes[-1].set_xticks([0, 1, 2, 3])
#    fig.text(0.93, 0.028, '>3', fontsize=9.7)
#    fig.text(0.9, 0.05, '/')
#    fig.text(0.06, 0.25, 'Which cheater has\na greater advantage?', fontsize=10, ma='center', style='italic')
#    fig.text(0.4, 0.25, 'Which cheater hurts\nthe WT more?', fontsize=10, ma='center', style='italic')
#    fig.text(0.7, 0.25, 'Which cheater can grow\nwith which cheater?', fontsize=10, ma='center', style='italic')
#    fig.text(0.18, 0.9, fix_symbols('W$_\mathrm{WT|WT}$\nfixed at 1'), ha='center', va='center', fontsize=12)
#    fig.text(0.5, 0.72, fix_symbols('W$_\mathrm{\Delta1764|\Delta1764}$\nfixed at 0'), ha='center', va='center', fontsize=12)
#    prior_patch = mpatches.Patch(color='#267BB8', label='Prior', alpha=0.2)
#    posterior_patch = mpatches.Patch(color='#267BB8', label='Posterior')
#    fig.legend(handles=[posterior_patch, prior_patch], bbox_to_anchor=(0.56,0.45), loc='upper center', borderaxespad=0., ncol=2)
#    line = plt.Line2D([0.1,0.95],[0.31, 0.31], transform=fig.transFigure, color="black", linestyle=':')
#    fig.add_artist(line)
#    posterior_patch2 = mpatches.Patch(color='green', label='Posterior')
#    axes[13].legend(handles=[posterior_patch2], bbox_to_anchor=(0.42,-0.7), loc='upper center', borderaxespad=0.)
#    fig.savefig(outpath, dpi=800, bbox_inches='tight')
#    return
    

#def prior_posterior_graph2(data, outpath):
#    data['del_with_wt_syn_with_wt_ratio'] = data.del_with_wt_p / data.syn_with_wt_p
#    #data['wt_with_del_wt_with_syn_ratio'] = data.wt_with_del_p / data.wt_with_syn_p
#    #data['del_with_syn_syn_with_del_dratio'] = (data.del_with_syn_p + 0.05) / (data.syn_with_del_p + 0.05)
#    data['del_with_syn_syn_with_del_dratio'] = (data.del_with_syn_p) / (data.syn_with_del_p)
#    data.loc[data['del_with_syn_syn_with_del_dratio'] > 3, 'del_with_syn_syn_with_del_dratio'] = 4.5
#    #data['triple_del_triple_syn_dratio'] = (data.triple_del + 0.05) / (data.triple_syn + 0.05)
#    data['triple_del_triple_syn_dratio'] = (data.triple_del) / (data.triple_syn)
#    data.loc[data['triple_del_triple_syn_dratio'] > 3, 'triple_del_triple_syn_dratio'] = 4.5
#
#    fig, axes = plt.subplots(nrows=6, ncols=3)
#    fig.set_size_inches(6,8.5)
#    fig.tight_layout()
#    axes = axes.flatten()
#    fig.subplots_adjust(hspace=0.6)
#    fig.text(0.5, -0.02, 'Relationship', ha='center', fontsize = 14)
#    fig.text(0.5, 0.32, 'Payoff', ha='center', fontsize = 14)
#    fig.text(-0.09, 0.65, 'Frequency', va='center', rotation='vertical', fontsize = 14)
#    fig.text(-0.09, 0.12, 'Frequency', va='center', rotation='vertical', fontsize = 14)
#    for i, column in enumerate(['wt_with_wt_p', 'wt_with_del_p', 'wt_with_syn_p', 'del_with_wt_p', 'del_with_del_p', 'del_with_syn_p', 'syn_with_wt_p', 'syn_with_del_p', 'syn_with_syn_p', 'triple_wt', 'triple_del', 'triple_syn', None, None, None, 'del_with_wt_syn_with_wt_ratio', 'triple_del_triple_syn_dratio', 'del_with_syn_syn_with_del_dratio']):
#        if column:
#            if column in ['wt_with_del_p', 'syn_with_syn_p']:
#                axes[i].hist(np.linspace(0,1,len(data)), alpha=0.2, color='#267BB8', bins=np.linspace(0, 4, 21))
#            elif column in ['del_with_wt_p', 'del_with_syn_p', 'syn_with_wt_p', 'syn_with_del_p', 'wt_with_syn_p', 'triple_wt', 'triple_syn', 'triple_del']:
#                axes[i].hist(np.linspace(0,4,len(data)), alpha=0.2, color='#267BB8', bins=np.linspace(0, 4, 21))
#            
#            if column not in ['del_with_del_p', 'wt_with_wt_p']:
#                if i < 12:
#                    data[column].plot.hist(ax=axes[i], color='#267BB8', bins=np.linspace(0, 4, 21))
#                    axes[i].set_xticks([0,1,2,3,4])
#                    axes[i].set_xlim(-0.1,4.1)
#                    axes[i].axvline(data[column].median(), color='black', linestyle='dashed', linewidth=1)
#                else:
#                    data[column].plot.hist(ax=axes[i], color='green')
#                axes[i].set_xlabel(None)
#                axes[i].set_ylabel(None)
#                
#            else:   
#                axes[i].set_visible(False)
#                            
#            if 'with' in column and i < 9:
#                title = 'W$_\mathrm{' + '|'.join(column.replace('_p', '').split('_with_')).replace('wt', 'WT').replace('del', '\Delta1764').replace('syn', 'A1664G') + '}$'
#                axes[i].set_title(fix_symbols(title), fontsize=12)
#            else:
#                #titles = {'del_with_wt_syn_with_wt_ratio':r'$\frac{\mathrm{W_\mathrm{\Delta1764|WT}}}{\mathrm{W_\mathrm{A1664G|WT}}}$', 'wt_with_del_wt_with_syn_ratio':r'$\frac{\mathrm{W_\mathrm{WT|\Delta1764}}}{\mathrm{W_\mathrm{WT|A1664G}}}$', 'del_with_syn_syn_with_del_dratio':r'$\frac{\mathrm{W_\mathrm{\Delta1764|A1664G}} + 0.05}{\mathrm{W_\mathrm{A1664G|\Delta1764}} + 0.05}$'}
#                titles = {'del_with_wt_syn_with_wt_ratio':r'$\frac{\mathrm{W_\mathrm{\Delta1764|WT}}}{\mathrm{W_\mathrm{A1664G|WT}}}$', 'wt_with_del_wt_with_syn_ratio':r'$\frac{\mathrm{W_\mathrm{WT|\Delta1764}}}{\mathrm{W_\mathrm{WT|A1664G}}}$', 'del_with_syn_syn_with_del_dratio':r'$\frac{\mathrm{W_\mathrm{\Delta1764|A1664G}}}{\mathrm{W_\mathrm{A1664G|\Delta1764}}}$', 'triple_wt':'$\mathrm{W_\mathrm{WT|\Delta1764,A1664G}}$', 'triple_del':'$\mathrm{W_\mathrm{\Delta1764|WT,A1664G}}$', 'triple_syn':'$\mathrm{W_\mathrm{A1664G|WT,\Delta1764}}$', 'triple_del_triple_syn_dratio':r'$\frac{\mathrm{W_\mathrm{\Delta1764|WT,A1664G}}}{\mathrm{W_\mathrm{A1664G|WT,\Delta1764}}}$'}
#                axes[i].set_title(fix_symbols(titles[column]), fontsize=14)
#            axes[i].grid(which='major', alpha=0.2, linestyle='-')
#        else:
#            axes[i].set_visible(False)
#    axes[-3].set_xticks([1.0, 1.25, 1.5, 1.75])
#    axes[-1].set_xticks([0, 1, 2, 3])
#    axes[-2].set_xticks([0, 1, 2, 3])
#    fig.text(0.91, 0.023, '>3', fontsize=9.5)
#    fig.text(0.9, 0.039, '/')
#    fig.text(0.59, 0.023, '>3', fontsize=9.5)
#    fig.text(0.58, 0.039, '/')
#    fig.text(0.06, 0.21, 'Which cheater has\na greater advantage?', fontsize=10, ma='center', style='italic')
#    fig.text(0.53, 0.21, 'Which cheater can grow\nwith which cheater?', fontsize=10, ma='center', style='italic')
#    fig.text(0.18, 0.94, fix_symbols('W$_\mathrm{WT|WT}$\nfixed at 1'), ha='center', va='center', fontsize=12)
#    fig.text(0.5, 0.76, fix_symbols('W$_\mathrm{\Delta1764|\Delta1764}$\nfixed at 0'), ha='center', va='center', fontsize=12)
#    prior_patch = mpatches.Patch(color='#267BB8', label='Prior', alpha=0.2)
#    posterior_patch = mpatches.Patch(color='#267BB8', label='Posterior')
#    fig.legend(handles=[posterior_patch, prior_patch], bbox_to_anchor=(0.62,0.34), loc='upper center', borderaxespad=0., ncol=2)
#    line = plt.Line2D([0.1,0.95],[0.26, 0.26], transform=fig.transFigure, color="black", linestyle=':')
#    fig.add_artist(line)
#    posterior_patch2 = mpatches.Patch(color='green', label='Posterior')
#    axes[13].legend(handles=[posterior_patch2], bbox_to_anchor=(0.42,-0.7), loc='upper center', borderaxespad=0.)
#    fig.savefig(outpath, dpi=800, bbox_inches='tight')
#    return
#
#def prior_posterior_graph3(data, outpath):
#    data['del_with_wt_syn_with_wt_ratio'] = data.del_with_wt_p / data.syn_with_wt_p
#    #data['wt_with_del_wt_with_syn_ratio'] = data.wt_with_del_p / data.wt_with_syn_p
#    #data['del_with_syn_syn_with_del_dratio'] = (data.del_with_syn_p + 0.05) / (data.syn_with_del_p + 0.05)
#    data['del_with_syn_syn_with_del_dratio'] = (data.del_with_syn_p) / (data.syn_with_del_p)
#    data.loc[data['del_with_syn_syn_with_del_dratio'] > 3, 'del_with_syn_syn_with_del_dratio'] = 4.5
#    #data['triple_del_triple_syn_dratio'] = (data.triple_del + 0.05) / (data.triple_syn + 0.05)
#    data['triple_del_triple_syn_dratio'] = (data.triple_del) / (data.triple_syn)
#    data.loc[data['triple_del_triple_syn_dratio'] > 3, 'triple_del_triple_syn_dratio'] = 4.5
#
#    fig, axes = plt.subplots(nrows=6, ncols=3)
#    fig.set_size_inches(6,8.5)
#    fig.tight_layout()
#    axes = axes.flatten()
#    fig.subplots_adjust(hspace=0.6)
#    fig.text(0.5, -0.02, 'Relationship', ha='center', fontsize = 14)
#    fig.text(0.5, 0.32, 'Payoff', ha='center', fontsize = 14)
#    fig.text(-0.09, 0.65, 'Frequency', va='center', rotation='vertical', fontsize = 14)
#    fig.text(-0.09, 0.12, 'Frequency', va='center', rotation='vertical', fontsize = 14)
#    
#    temp_fig, temp_ax = plt.subplots(nrows=1, ncols=1)
#    
#    for i, column in enumerate(['wt_with_wt_p', 'wt_with_del_p', 'wt_with_syn_p', 'del_with_wt_p', 'del_with_del_p', 'del_with_syn_p', 'syn_with_wt_p', 'syn_with_del_p', 'syn_with_syn_p', 'triple_wt', 'triple_del', 'triple_syn', None, None, None, 'del_with_wt_syn_with_wt_ratio', 'triple_del_triple_syn_dratio', 'del_with_syn_syn_with_del_dratio']):
#        if column:
#            
#            if column not in ['del_with_del_p', 'wt_with_wt_p']:
#                if i < 12:
#                    data[column].plot.hist(ax=axes[i], color='#267BB8', bins=np.linspace(0, 4, 21))
#                    axes[i].set_xticks([0,1,2,3,4])
#                    axes[i].set_xlim(-0.1,4.1)
#                    axes[i].axvline(data[column].median(), color='black', linestyle='dashed', linewidth=1)
#                else:
#                    data[column].plot.hist(ax=axes[i], color='green')
#                axes[i].set_xlabel(None)
#                axes[i].set_ylabel(None)
#                axes[i].set_yticks([])
#                
#            else:   
#                axes[i].set_visible(False)
#            
#            if column in ['wt_with_del_p', 'syn_with_syn_p']:
#                a,b,c = temp_ax.hist(data[column], bins=np.linspace(0, 4, 21))
#                axes[i].hist(np.linspace(0,1,max(a)/2), alpha=0.4, color='#C5C6C6', bins=1)
#            elif column in ['del_with_wt_p', 'del_with_syn_p', 'syn_with_wt_p', 'syn_with_del_p', 'wt_with_syn_p', 'triple_wt', 'triple_syn', 'triple_del']:
#                a,b,c = temp_ax.hist(data[column], bins=np.linspace(0, 4, 21))
#                axes[i].hist(np.linspace(0,4,max(a)/2), alpha=0.4, color='#C5C6C6', bins=1)
#                
#            if 'with' in column and i < 9:
#                title = 'W$_\mathrm{' + '|'.join(column.replace('_p', '').split('_with_')).replace('wt', 'WT').replace('del', '\Delta1764').replace('syn', 'A1664G') + '}$'
#                axes[i].set_title(fix_symbols(title), fontsize=12)
#            else:
#                #titles = {'del_with_wt_syn_with_wt_ratio':r'$\frac{\mathrm{W_\mathrm{\Delta1764|WT}}}{\mathrm{W_\mathrm{A1664G|WT}}}$', 'wt_with_del_wt_with_syn_ratio':r'$\frac{\mathrm{W_\mathrm{WT|\Delta1764}}}{\mathrm{W_\mathrm{WT|A1664G}}}$', 'del_with_syn_syn_with_del_dratio':r'$\frac{\mathrm{W_\mathrm{\Delta1764|A1664G}} + 0.05}{\mathrm{W_\mathrm{A1664G|\Delta1764}} + 0.05}$'}
#                titles = {'del_with_wt_syn_with_wt_ratio':r'$\frac{\mathrm{W_\mathrm{\Delta1764|WT}}}{\mathrm{W_\mathrm{A1664G|WT}}}$', 'wt_with_del_wt_with_syn_ratio':r'$\frac{\mathrm{W_\mathrm{WT|\Delta1764}}}{\mathrm{W_\mathrm{WT|A1664G}}}$', 'del_with_syn_syn_with_del_dratio':r'$\frac{\mathrm{W_\mathrm{\Delta1764|A1664G}}}{\mathrm{W_\mathrm{A1664G|\Delta1764}}}$', 'triple_wt':'$\mathrm{W_\mathrm{WT|\Delta1764,A1664G}}$', 'triple_del':'$\mathrm{W_\mathrm{\Delta1764|WT,A1664G}}$', 'triple_syn':'$\mathrm{W_\mathrm{A1664G|WT,\Delta1764}}$', 'triple_del_triple_syn_dratio':r'$\frac{\mathrm{W_\mathrm{\Delta1764|WT,A1664G}}}{\mathrm{W_\mathrm{A1664G|WT,\Delta1764}}}$'}
#                axes[i].set_title(fix_symbols(titles[column]), fontsize=14)
#            axes[i].grid(which='major', alpha=0.2, linestyle='-')
#        else:
#            axes[i].set_visible(False)
#    axes[-3].set_xticks([1.0, 1.25, 1.5, 1.75])
#    axes[-1].set_xticks([0, 1, 2, 3])
#    axes[-2].set_xticks([0, 1, 2, 3])
#    fig.text(0.91, 0.023, '>3', fontsize=9.5)
#    fig.text(0.9, 0.039, '/')
#    fig.text(0.59, 0.023, '>3', fontsize=9.5)
#    fig.text(0.58, 0.039, '/')
#    fig.text(0.06, 0.21, 'Which cheater has\na greater advantage?', fontsize=10, ma='center', style='italic')
#    fig.text(0.53, 0.21, 'Which cheater can grow\nwith which cheater?', fontsize=10, ma='center', style='italic')
#    fig.text(0.18, 0.94, fix_symbols('W$_\mathrm{WT|WT}$\nfixed at 1'), ha='center', va='center', fontsize=12)
#    fig.text(0.5, 0.76, fix_symbols('W$_\mathrm{\Delta1764|\Delta1764}$\nfixed at 0'), ha='center', va='center', fontsize=12)
#    prior_patch = mpatches.Patch(color='#C5C6C6', label='Prior', alpha=0.4)
#    posterior_patch = mpatches.Patch(color='#267BB8', label='Posterior')
#    fig.legend(handles=[posterior_patch, prior_patch], bbox_to_anchor=(0.62,0.34), loc='upper center', borderaxespad=0., ncol=2)
#    line = plt.Line2D([0.1,0.95],[0.26, 0.26], transform=fig.transFigure, color="black", linestyle=':')
#    fig.add_artist(line)
#    posterior_patch2 = mpatches.Patch(color='green', label='Posterior')
#    axes[13].legend(handles=[posterior_patch2], bbox_to_anchor=(0.42,-0.7), loc='upper center', borderaxespad=0.)
#    #fig.savefig(outpath, dpi=800, bbox_inches='tight')
#    return

def prior_posterior_graph4(data, outpath):
    data['del_with_wt_syn_with_wt_ratio'] = data.del_with_wt_p / data.syn_with_wt_p
    #data['wt_with_del_wt_with_syn_ratio'] = data.wt_with_del_p / data.wt_with_syn_p
    #data['del_with_syn_syn_with_del_dratio'] = (data.del_with_syn_p + 0.05) / (data.syn_with_del_p + 0.05)
    data['del_with_syn_syn_with_del_dratio'] = (data.del_with_syn_p) / (data.syn_with_del_p)
    data.loc[data['del_with_syn_syn_with_del_dratio'] > 3, 'del_with_syn_syn_with_del_dratio'] = 4.5
    #data['triple_del_triple_syn_dratio'] = (data.triple_del + 0.05) / (data.triple_syn + 0.05)
    data['triple_del_triple_syn_dratio'] = (data.triple_del) / (data.triple_syn)
    data.loc[data['triple_del_triple_syn_dratio'] > 3, 'triple_del_triple_syn_dratio'] = 4.5
    data['syn_with_syn_syn_with_del_dratio'] = (data.syn_with_syn_p / data.syn_with_del_p)
    data.loc[data['syn_with_syn_syn_with_del_dratio'] > 3, 'syn_with_syn_syn_with_del_dratio'] = 4.5

    fig, axes = plt.subplots(nrows=4, ncols=3)
    fig.set_size_inches(6,6)
    fig.tight_layout()
    axes = axes.flatten()
    fig.subplots_adjust(wspace=0.2)
    fig.subplots_adjust(hspace=0.6)
    fig.text(0.51, -0.02, 'Fitness', ha='center', fontsize = 14)
    fig.text(0, 0.5, 'Frequency', va='center', rotation='vertical', fontsize = 14)
    
    temp_fig, temp_ax = plt.subplots(nrows=1, ncols=1)

    
    for i, column in enumerate(['wt_with_wt_p', 'wt_with_del_p', 'wt_with_syn_p', 'del_with_wt_p', 'del_with_del_p', 'del_with_syn_p', 'syn_with_wt_p', 'syn_with_del_p', 'syn_with_syn_p', 'triple_wt', 'triple_del', 'triple_syn']):
        if column:            
            if column not in ['del_with_del_p', 'wt_with_wt_p']:
                data[column].plot.hist(ax=axes[i], color='#267BB8', bins=np.linspace(0, 4, 21))
                axes[i].set_xticks([0,1,2,3,4])
                axes[i].set_xlim(-0.1,4.1)
                axes[i].axvline(data[column].median(), color='black', linestyle='dashed', linewidth=1)
                axes[i].set_xlabel(None)
                axes[i].set_ylabel(None)
                axes[i].set_yticks([])
            else:
                axes[i].set_visible(False)
                            
            if column in ['wt_with_del_p', 'syn_with_syn_p']:
                a,b,c = temp_ax.hist(data[column], bins=np.linspace(0, 4, 21))
                axes[i].hist(np.linspace(0,1,max(a)/2), alpha=0.4, color='#C5C6C6', bins=1)
            elif column in ['del_with_wt_p', 'del_with_syn_p', 'syn_with_wt_p', 'syn_with_del_p', 'wt_with_syn_p', 'triple_wt', 'triple_syn', 'triple_del']:
                a,b,c = temp_ax.hist(data[column], bins=np.linspace(0, 4, 21))
                axes[i].hist(np.linspace(0,4,max(a)/2), alpha=0.4, color='#C5C6C6', bins=1)
                
            if 'with' in column:
                title = 'W$_\mathrm{' + '|'.join(column.replace('_p', '').split('_with_')).replace('wt', 'WT').replace('del', '\Delta1764').replace('syn', 'A1664G') + '}$'
                axes[i].set_title(fix_symbols(title), fontsize=14)
            else:
                titles = {'del_with_wt_syn_with_wt_ratio':r'$\frac{\mathrm{W_\mathrm{\Delta1764|WT}}}{\mathrm{W_\mathrm{A1664G|WT}}}$', 'wt_with_del_wt_with_syn_ratio':r'$\frac{\mathrm{W_\mathrm{WT|\Delta1764}}}{\mathrm{W_\mathrm{WT|A1664G}}}$', 'del_with_syn_syn_with_del_dratio':r'$\frac{\mathrm{W_\mathrm{\Delta1764|A1664G}}}{\mathrm{W_\mathrm{A1664G|\Delta1764}}}$', 'triple_wt':'$\mathrm{W_\mathrm{WT|\Delta1764,A1664G}}$', 'triple_del':'$\mathrm{W_\mathrm{\Delta1764|WT,A1664G}}$', 'triple_syn':'$\mathrm{W_\mathrm{A1664G|WT,\Delta1764}}$', 'triple_del_triple_syn_dratio':r'$\frac{\mathrm{W_\mathrm{\Delta1764|WT,A1664G}}}{\mathrm{W_\mathrm{A1664G|WT,\Delta1764}}}$'}
                axes[i].set_title(fix_symbols(titles[column]), fontsize=14)
            axes[i].grid(which='major', alpha=0.3, linestyle='-')


    fig.text(0.2, 0.89, fix_symbols('W$_\mathrm{WT|WT}$\nfixed at 1'), ha='center', va='center', fontsize=12)
    fig.text(0.51, 0.65, fix_symbols('W$_\mathrm{\Delta1764|\Delta1764}$\nfixed at 0'), ha='center', va='center', fontsize=12)
    prior_patch = mpatches.Patch(color='#C5C6C6', label='Prior', alpha=0.4)
    posterior_patch = mpatches.Patch(color='#267BB8', label='Posterior')
    axes[10].legend(handles=[posterior_patch, prior_patch], bbox_to_anchor=(0.5,-0.7), loc='upper center', borderaxespad=0., ncol=2)
    #########
    fig2, axes2 = plt.subplots(nrows=1, ncols=4)
    fig2.set_size_inches(7,1.4)
    fig2.tight_layout()
    fig2.subplots_adjust(wspace=0.35)
    axes2 = axes2.flatten()
    fig2.text(0.52, -0.07, 'Posterior Relationship', ha='center', fontsize = 14)
    fig2.text(-.0, 0.6, 'Frequency', va='center', rotation='vertical', fontsize = 14)
    
    for i,column in enumerate(['del_with_wt_syn_with_wt_ratio', 'del_with_syn_syn_with_del_dratio', 'triple_del_triple_syn_dratio', 'syn_with_syn_syn_with_del_dratio']):
        data[column].plot.hist(ax=axes2[i], color='green')
        axes2[i].axvline(data[column].median(), color='black', linestyle='dashed', linewidth=1)
        axes2[i].set_xlabel(None)
        axes2[i].set_ylabel(None)
        axes2[i].set_yticks([])
        axes2[i].grid(which='major', alpha=0.3, linestyle='-')
        titles = {'del_with_wt_syn_with_wt_ratio':r'$\frac{\mathrm{W_\mathrm{\Delta1764|WT}}}{\mathrm{W_\mathrm{A1664G|WT}}}$', 'wt_with_del_wt_with_syn_ratio':r'$\frac{\mathrm{W_\mathrm{WT|\Delta1764}}}{\mathrm{W_\mathrm{WT|A1664G}}}$', 'del_with_syn_syn_with_del_dratio':r'$\frac{\mathrm{W_\mathrm{\Delta1764|A1664G}}}{\mathrm{W_\mathrm{A1664G|\Delta1764}}}$', 'triple_wt':'$\mathrm{W_\mathrm{WT|\Delta1764,A1664G}}$', 'triple_del':'$\mathrm{W_\mathrm{\Delta1764|WT,A1664G}}$', 'triple_syn':'$\mathrm{W_\mathrm{A1664G|WT,\Delta1764}}$', 'triple_del_triple_syn_dratio':r'$\frac{\mathrm{W_\mathrm{\Delta1764|WT,A1664G}}}{\mathrm{W_\mathrm{A1664G|WT,\Delta1764}}}$', 'syn_with_syn_syn_with_del_dratio':r'$\frac{\mathrm{W_\mathrm{A1664G|A1664G}}}{\mathrm{W_\mathrm{A1664G|\Delta1764}}}$'}
        axes2[i].set_title(fix_symbols(titles[column]), fontsize=14)
    
    #posterior_patch2 = mpatches.Patch(color='green', label='Posterior')
    #axes2[1].legend(handles=[posterior_patch2], bbox_to_anchor=(0.5,-0.7), loc='upper center', borderaxespad=0.)
    axes2[0].set_xticks([1.0, 1.25, 1.5, 1.75])
    axes2[1].set_xticks([0, 1, 2, 3])
    axes2[2].set_xticks([0, 1, 2, 3])
    axes2[3].set_xticks([0, 1, 2, 3])
    fig2.text(0.93, 0.14, '>3', fontsize=9.5)
    fig2.text(0.92, 0.24, '/')
    fig2.text(0.69, 0.14, '>3', fontsize=9.5)
    fig2.text(0.68, 0.24, '/')
    fig2.text(0.45, 0.14, '>3', fontsize=9.5)
    fig2.text(0.44, 0.24, '/')
    fig2.text(0.06, 1.2, 'Which cheater has\na greater advantage?', fontsize=10, ma='center', style='italic')
    fig2.text(0.39, 1.2, 'Which cheater can grow\nwith which cheater?', fontsize=10, ma='center', style='italic')
    fig2.text(0.73, 1.2, 'Does the synonymous cheater\nexploit the deletion cheater?', fontsize=10, ma='center', style='italic')
    #fig.savefig(outpath, dpi=800, bbox_inches='tight')
    fig2.savefig(outpath.replace('.png', '.relationships.png'), dpi=800, bbox_inches='tight')
    return

def fix_symbols(string):
    return string.replace('\Delta1764', r'\bigtriangleup').replace('A1664G', r'\bigcirc')

#def prior_posterior_triple_only(data, outpath):
#    fig, axes = plt.subplots(nrows=1, ncols=3)
#    fig.set_size_inches(6,1.5)
#    fig.tight_layout()
#    axes = axes.flatten()
#    fig.text(0.5, -0.02, 'Payoff', ha='center', fontsize = 14)
#    fig.text(-0.05, 0.5, 'Frequency', va='center', rotation='vertical', fontsize = 14)
#    for i, column in enumerate(['triple_wt', 'triple_del', 'triple_syn']):
#        if column in []:
#            axes[i].hist(np.linspace(0,1,len(data)), alpha=0.2, color='#267BB8', bins=np.linspace(0, 4, 21))
#        elif column in ['triple_wt', 'triple_del', 'triple_syn']:
#            axes[i].hist(np.linspace(0,4,len(data)), alpha=0.2, color='#267BB8', bins=np.linspace(0, 4, 21))
#        data[column].plot.hist(ax=axes[i], color='#267BB8', bins=np.linspace(0, 4, 21))
#        axes[i].set_xlabel(None)
#        axes[i].set_ylabel(None)
#        axes[i].grid(which='major', alpha=0.2, linestyle='-')
#        axes[i].set_xticks([0,1,2,3,4])
#        axes[i].set_xlim(-0.1,4.1)
#        titles= {'triple_wt':'$\mathrm{W_\mathrm{WT|\Delta1764,A1664G}}$', 'triple_del':'$\mathrm{W_\mathrm{\Delta1764|WT,A1664G}}$', 'triple_syn':'$\mathrm{W_\mathrm{A1664G|WT,\Delta1764}}$'}
#        axes[i].set_title(fix_symbols(titles[column]), fontsize=14)
#    prior_patch = mpatches.Patch(color='#267BB8', label='Prior', alpha=0.2)
#    posterior_patch = mpatches.Patch(color='#267BB8', label='Posterior')
#    fig.legend(handles=[posterior_patch, prior_patch], bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
#    fig.savefig(outpath, dpi=800, bbox_inches='tight')
#    return

prior_posterior_graph4(get_best_parameters('X:/volume2/noam/cheater_model/abc_smc/new_lim_only_wt_with_del/37B_abc.txt', 10000), 'X:/volume2/noam/cheater_model/abc_smc/new_lim_only_wt_with_del/37B_abc.txt.prior_posterior.png')
prior_posterior_graph4(get_best_parameters('X:/volume2/noam/cheater_model/abc_smc/new_lim_only_wt_with_del/37A_abc.txt', 10000), 'X:/volume2/noam/cheater_model/abc_smc/new_lim_only_wt_with_del/37A_abc.txt.prior_posterior.png')
    

data = get_best_parameters('X:/volume2/noam/cheater_model/abc_smc/new_lim_only_wt_with_del/37B_abc.txt', 10000)
data['del_with_wt_syn_with_wt_ratio'] = data.del_with_wt_p / data.syn_with_wt_p
data['del_with_syn_syn_with_del_dratio'] = (data.del_with_syn_p) / (data.syn_with_del_p)
#data.loc[data['del_with_syn_syn_with_del_dratio'] > 3, 'del_with_syn_syn_with_del_dratio'] = 4.5
data['triple_del_triple_syn_dratio'] = (data.triple_del) / (data.triple_syn)
#data.loc[data['triple_del_triple_syn_dratio'] > 3, 'triple_del_triple_syn_dratio'] = 4.5
data['syn_with_syn_syn_with_del_dratio'] = (data.syn_with_syn_p / data.syn_with_del_p)
data.loc[data['syn_with_syn_syn_with_del_dratio'] > 3, 'syn_with_syn_syn_with_del_dratio'] = 4.5




## single cheater prior posterior graph
def prior_posterior_graph_single_cheater(astroabc_output, n, outpath):
    data = pd.read_csv(astroabc_output, '\t', index_col=False)
    data = data.tail(n)
    data['wt_with_wt_p'] = 1
    data['del_with_del_p'] = 0
    data = data.rename(columns={
     'param#0 ':'wt_with_del_p',
     ' param#1 ':'del_with_wt_p',})
    fig, axes = plt.subplots(nrows=2, ncols=2)
    fig.set_size_inches(4,2.75)
    fig.tight_layout()
    axes = axes.flatten()
    
    temp_fig, temp_ax = plt.subplots(nrows=1, ncols=1)
    
    #fig.subplots_adjust(hspace=0.5)
    for i, column in enumerate(['wt_with_wt_p', 'wt_with_del_p','del_with_wt_p', 'del_with_del_p']):
        if column:
            if column in ['wt_with_del_p']:
                a,b,c = temp_ax.hist(data[column], bins=np.linspace(0, 4, 21))
                axes[i].hist(np.linspace(0,1,max(a)/2), alpha=0.4, color='#C5C6C6', bins=1)
            elif column in ['del_with_wt_p']:
                a,b,c = temp_ax.hist(data[column], bins=np.linspace(0, 4, 21))
                axes[i].hist(np.linspace(0,4,max(a)/2), alpha=0.4, color='#C5C6C6', bins=1)
            if column not in ['del_with_del_p', 'wt_with_wt_p']:
                axes[i].set_xticks([0,1,2,3,4])
                axes[i].set_xlim(-0.1,4.1)
                data[column].plot.hist(ax=axes[i], color='#267BB8', bins=np.linspace(0, 4, 21))
                axes[i].set_xlabel(None)
                axes[i].set_ylabel(None)
                axes[i].set_yticks([])
                #axes[i].axvline(data[column].median(), color='black', linestyle='dashed', linewidth=1)
            else:    
                axes[i].set_visible(False)
            if 'with' in column:
                title = 'W$_\mathrm{' + '|'.join(column.replace('_p', '').split('_with_')).replace('wt', 'WT').replace('del', '\Delta1764').replace('syn', 'A1664G') + '}$'
                axes[i].set_title(fix_symbols(title), fontsize=12)
            axes[i].grid(which='major', alpha=0.3, linestyle='-') 
    prior_patch = mpatches.Patch(color='#C5C6C6', label='Prior', alpha=0.4)
    posterior_patch = mpatches.Patch(color='#267BB8', label='Posterior')
    axes[2].legend(handles=[posterior_patch, prior_patch], bbox_to_anchor=(1.1,-0.7), loc='upper center', borderaxespad=0., ncol=2)
    fig.text(0.26, 0.8, fix_symbols('W$_\mathrm{WT|WT}$\nfixed at 1'), ha='center', va='center', fontsize=12)
    fig.text(0.76, 0.3, fix_symbols('W$_\mathrm{\Delta1764|\Delta1764}$\nfixed at 0'), ha='center', va='center', fontsize=12)
    text1 = fig.text(0.5, -0.12, 'Fitness\n', ha='center', fontsize = 14)
    text2 = fig.text(0, 0.5, 'Frequency', va='center', rotation='vertical', fontsize = 14)
    fig.savefig(outpath, dpi=800, bbox_inches='tight')
    return

prior_posterior_graph_single_cheater('X:/volume2/noam/cheater_model/abc_smc/37B_single_cheater.txt', 10000, 'X:/volume2/noam/cheater_model/abc_smc/37B_single_cheater.txt.prior_posterior.png')
prior_posterior_graph_single_cheater('X:/volume2/noam/cheater_model/abc_smc/37A_single_cheater.txt', 10000, 'X:/volume2/noam/cheater_model/abc_smc/37A_single_cheater.txt.prior_posterior.png')
