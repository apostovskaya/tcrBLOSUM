import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import matplotlib.ticker as mtick
import seaborn as sns

main_dir = os.getcwd()
folder_in = f'{main_dir}/results/tcrBLOSUMmtx/tcrBLOSUMmtx_2024-02-06'

# Composition in percent for the complete database, UniProt Jan 2024
uniprot = {'A': 9.01, 'Q': 3.80, 'L': 9.84, 'S': 6.84, 'R': 5.84, 'E': 6.24,
           'K': 4.94, 'T': 5.55, 'N': 3.80, 'G': 7.26, 'M': 2.33, 'W': 1.30,
           'D': 5.47, 'H': 2.22, 'F': 3.88, 'Y': 2.87, 'C': 1.29, 'I': 5.52,
           'P': 5.00, 'V': 6.86}
uniprot_df = pd.DataFrame.from_dict(uniprot, orient='index',
                                    columns=['Frequency_UniProt'])
uniprot_df['Frequency_UniProt'] = uniprot_df['Frequency_UniProt'] / 100

alpha = [('C', 0.0005), ('H', 0.004), ('W', 0.0045), ('M', 0.0194),
         ('P', 0.0207), ('F', 0.0208), ('E', 0.0225), ('Q', 0.0319),
         ('R', 0.0372), ('D', 0.0409), ('Y', 0.0479), ('I', 0.0491),
         ('K', 0.0528), ('V', 0.0602), ('T', 0.0702), ('N', 0.0764),
         ('S', 0.0816), ('L', 0.0884), ('A', 0.1269), ('G', 0.1439)]
alpha_df = pd.DataFrame(alpha, columns=['AA', 'Frequency_CDR3a']).set_index(
    'AA')

beta = [('C', 0.0004), ('M', 0.0042), ('K', 0.0077), ('W', 0.0081),
        ('H', 0.0137), ('I', 0.0142), ('V', 0.0239), ('N', 0.0348),
        ('R', 0.0366), ('P', 0.0372), ('D', 0.0384), ('F', 0.0445),
        ('A', 0.0485), ('L', 0.0541), ('E', 0.0761), ('Y', 0.0778),
        ('T', 0.0793), ('Q', 0.0849), ('G', 0.1219), ('S', 0.1937)]
beta_df = pd.DataFrame(beta, columns=['AA', 'Frequency_CDR3b']).set_index('AA')

cd4 = [('C', 0.0016), ('W', 0.0054), ('M', 0.0073), ('H', 0.0123),
       ('I', 0.0218), ('K', 0.0243), ('P', 0.0296), ('V', 0.0297),
       ('F', 0.0327), ('D', 0.0333), ('E', 0.0456), ('N', 0.047),
       ('R', 0.0552), ('Q', 0.063), ('L', 0.0652), ('Y', 0.0655),
       ('A', 0.0771), ('T', 0.0878), ('G', 0.131), ('S', 0.1644)]
cd4_df = pd.DataFrame(cd4, columns=['AA', 'Frequency_CDR3_CD4']).set_index(
    'AA')

cd8 = [('C', 0.0004), ('W', 0.0071), ('M', 0.0094), ('H', 0.0105),
       ('K', 0.0227), ('I', 0.0268), ('P', 0.0317), ('V', 0.0362),
       ('F', 0.0364), ('D', 0.0378), ('R', 0.0379), ('N', 0.0491),
       ('E', 0.0578), ('L', 0.0651), ('Q', 0.0665), ('Y', 0.0672),
       ('T', 0.076), ('A', 0.0763), ('G', 0.1319), ('S', 0.1529)]
cd8_df = pd.DataFrame(cd8, columns=['AA', 'Frequency_CDR3_CD8']).set_index(
    'AA')

cdr3s = [('C', 0.0004), ('W', 0.0073), ('M', 0.0076), ('H', 0.0115),
         ('K', 0.0176), ('I', 0.0219), ('V', 0.0318), ('P', 0.0336),
         ('R', 0.0367), ('D', 0.039), ('F', 0.0393), ('N', 0.0439),
         ('L', 0.0616), ('E', 0.0644), ('A', 0.0657), ('Y', 0.0713),
         ('Q', 0.0733), ('T', 0.0773), ('G', 0.1268), ('S', 0.1691)]
cdr3s_df = pd.DataFrame(cdr3s, columns=['AA', 'Frequency_CDR3s']).set_index(
    'AA')

aa_groups = {'A': ['aliphatic', 'grey'], 'R': ['basic', 'blue'],
             'N': ['amide', 'white'], 'D': ['acidic', 'red'],
             'C': ['sulfur', 'yellow'], 'Q': ['amide', 'white'],
             'E': ['acidic', 'red'], 'G': ['aliphatic', 'grey'],
             'H': ['basic', 'blue'], 'I': ['aliphatic', 'grey'],
             'L': ['aliphatic', 'grey'], 'K': ['basic', 'blue'],
             'M': ['sulfur', 'yellow'], 'F': ['aromatic', 'black'],
             'P': ['aliphatic', 'grey'], 'S': ['small_hydroxy', 'green'],
             'T': ['small_hydroxy', 'green'], 'W': ['aromatic', 'black'],
             'Y': ['aromatic', 'black'], 'V': ['aliphatic', 'grey']}

aa_summary = pd.DataFrame.from_dict(aa_groups, orient='index',
                                    columns=['Property', 'Color'])

aa_summary.loc[aa_summary['Color'] == 'black', 'Color'] = 'dimgrey'
aa_summary.loc[aa_summary['Color'] == 'grey', 'Color'] = 'silver'

df_list = [uniprot_df, alpha_df, beta_df, cd4_df, cd8_df]

for df in df_list:
    aa_summary = aa_summary.merge(df, left_index=True, right_index=True)
aa_summary.sort_values('Frequency_UniProt', ascending=False, inplace=True)
aa_summary['AminoAcid'] = aa_summary.index

data_sources = ['CDR3a', 'CDR3b', 'CDR3_CD4', 'CDR3_CD8', 'CDR3s']
# PLOTTING amino acid distribution in the dataset
source = 'CDR3s'
for source in data_sources:
    fig, ax = plt.subplots(figsize=(10, 6))
    # Set the width of the bars
    bar_width = 0.35
    # Set the position of the bars on the x-axis
    r1 = np.arange(len(aa_summary['AminoAcid']))
    r2 = [x + bar_width for x in r1]

    # Create bars with different colors
    plt.bar(r1, aa_summary['Frequency_UniProt'], bar_width, hatch='//',
            color=list(aa_summary['Color']), edgecolor='black')
    plt.bar(r2, aa_summary[f'Frequency_{source}'], bar_width, hatch='\\\\',
            color=list(aa_summary['Color']), edgecolor='black')
    plt.xticks(r1 + bar_width / 2, aa_summary['AminoAcid'])
    # increase bar width
    # change_width(ax, 0.75)

    # create legend manually
    df = aa_summary[['Property', 'Color']].drop_duplicates()
    df.sort_values(by='Property', inplace=True)
    # map names to colors
    cmap = dict(zip(list(df['Property']), list(df['Color'])))
    # create the rectangles for the legend
    # patches = [Patch(color=v, label=k) for k, v in cmap.items()]
    # # Iterate through the handles to add edge color to the legend markers
    # for ha in patches:
    #     ha.set_edgecolor("black")
    #
    # # add the legend
    # plt.legend(title='Amino acid property',
    #            labels=list(df['Property']),
    #            handles=patches, frameon=False,
    #            bbox_to_anchor=(0.05, 0.5), loc='center left')  # , bbox_to_anchor=(1.04, 0.5),
    # # loc='center left', borderaxespad=0, fontsize=15, frameon=False)

    # Create legend manually for bars
    legend_labels_bars = ['Data source', 'UniProt', source]
    legend_handles_bars = [Patch(facecolor='none'),
                           Patch(facecolor='white', edgecolor='black',
                                 hatch='//'),
                           Patch(facecolor='white', edgecolor='black',
                                 hatch='\\\\')]

    legend_labels_aa = list(df['Property'])
    legend_handles_aa = [Patch(color=cmap[prop], label=prop, edgecolor='black')
                         for prop in legend_labels_aa]
    for ha in legend_handles_aa:
        ha.set_edgecolor("black")

    # Combine handles and labels from both sets of legends
    legend_labels = legend_labels_aa + legend_labels_bars
    legend_handles = legend_handles_aa + legend_handles_bars

    # Add the combined legend outside the plot
    plt.legend(legend_handles, legend_labels, title='Amino acid property',
               frameon=False, bbox_to_anchor=(1.05, 1), loc='upper left')

    # plt.legend(legend_handles_bars, legend_labels_bars, title='Data source',
    #            frameon=False,
    #            bbox_to_anchor=(1.05, 0.5), loc='center left')

    plt.xlabel('amino acid', size=15)
    plt.ylabel('Frequency in the data', size=15)
    plt.ylim(0, 0.2)
    # plt.title(
    #     f'The frequency of each amino acid \nin the CDR3 sequences in the {prefix}',
    #     size=15)
    plt.tight_layout()
    # plt.show()
    plt.savefig(f'{folder_in}/Plots/aaDistrib_{source}.jpg', dpi=1200)
    plt.close()

aa_summary['CDR3/UniProt'] = aa_summary['Frequency_CDR3s'] / aa_summary[
    'Frequency_UniProt']
aa_summary.to_csv(f'{folder_in}/AA_frequency_summary.tsv', sep='\t')
# aa_summary.plot(x='AminoAcid', y='CDR3/UniProt', kind="bar", rot=0)

# plot up/down ratio of frequencies in CDR3s vs UniProt
baseline = 1
# Set colors based on whether values are above or below the baseline
colors = ['tab:orange' if x > baseline else 'tab:blue' for x in
          aa_summary['CDR3/UniProt']]
plt.bar(aa_summary['AminoAcid'],
        [x - baseline for x in aa_summary['CDR3/UniProt']], color=colors)
plt.gca().yaxis.set_major_formatter(
    mtick.FuncFormatter(lambda x, _: x + baseline))
plt.xlabel('amino acid', size=15)
plt.ylabel('Frequency in the data', size=15)
plt.tight_layout()
# plt.show()
plt.savefig(f'{folder_in}/Plots/AA_frequency_updown.jpg', dpi=1200)
plt.close()

# up down but color by property
# reorder AAs according to BLOSUM62
order = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F',
         'P', 'S', 'T', 'W', 'Y', 'V']
aa_summary['AA_sorter'] = pd.Categorical(aa_summary['AminoAcid'],
                                  ordered=True, categories=order)
aa_summary = aa_summary.sort_values('AA_sorter')

fig, ax = plt.subplots(figsize=(10, 6))
plt.bar(aa_summary['AminoAcid'],
        [x - baseline for x in aa_summary['CDR3s/UniProt']],
        color=list(aa_summary['Color']), edgecolor='black')
plt.gca().yaxis.set_major_formatter(
    mtick.FuncFormatter(lambda x, _: x + baseline))
df = aa_summary[['Property', 'Color']].drop_duplicates()
df.sort_values(by='Property', inplace=True)
# map names to colors
cmap = dict(zip(list(df['Property']), list(df['Color'])))
# create the rectangles for the legend
patches = [Patch(color=v, label=k) for k, v in cmap.items()]
# Iterate through the handles to add edge color to the legend markers
for ha in patches:
    ha.set_edgecolor("black")
# add the legend
plt.legend(title='Amino acid property',
           labels=list(df['Property']),
           handles=patches, frameon=False,
           bbox_to_anchor=(1.05, 1), loc='upper left')
plt.xlabel('amino acid', size=15)
plt.ylabel('AA frequency in CDR3s relative to UniProt', size=15)
plt.tight_layout()
# plt.show()
plt.savefig(f'{folder_in}/Plots/AA_frequency_updown_property_ordered.jpg',
            dpi=1200)
plt.close()

# per group
aa_summary.loc[aa_summary['Color'] == 'dimgrey', 'Color'] = 'black'
aa_summary.loc[aa_summary['Color'] == 'silver', 'Color'] = 'grey'
for source in data_sources:
    aa_summary[f'{source}/UniProt'] = aa_summary[f'Frequency_{source}'] \
                                      / aa_summary['Frequency_UniProt']
    aa_summary.sort_values(by=[f'Frequency_{source}'], inplace=True,
                           ascending=False)

    fig, ax = plt.subplots(figsize=(10, 6))
    plt.bar(aa_summary['AminoAcid'],
            [x - baseline for x in aa_summary[f'{source}/UniProt']],
            color=list(aa_summary['Color']), edgecolor='black')
    plt.gca().yaxis.set_major_formatter(
        mtick.FuncFormatter(lambda x, _: x + baseline))
    df = aa_summary[['Property', 'Color']].drop_duplicates()
    df.sort_values(by='Property', inplace=True)
    # map names to colors
    cmap = dict(zip(list(df['Property']), list(df['Color'])))
    # create the rectangles for the legend
    patches = [Patch(color=v, label=k) for k, v in cmap.items()]
    # Iterate through the handles to add edge color to the legend markers
    for ha in patches:
        ha.set_edgecolor("black")
    # add the legend
    plt.legend(title='Amino acid property',
               labels=list(df['Property']),
               handles=patches, frameon=False,
               bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.xlabel('amino acid', size=15)
    plt.ylabel(f'AA frequency in {source} relative to UniProt', size=15)
    plt.tight_layout()
    # plt.show()
    plt.savefig(f'{folder_in}/Plots/AA_frequency_updown_property_{source}.jpg',
                dpi=1200)
    plt.close()

# all together, stacked bar plot
fig, ax = plt.subplots(figsize=(10, 6))
ax = aa_summary.plot(x='AminoAcid', y=['Frequency_CDR3a', 'Frequency_CDR3b',
                                       'Frequency_CDR3_CD4',
                                       'Frequency_CDR3_CD8'],
                     kind="bar", rot=0, stacked=True, width=0.5,
                     color=['chocolate', 'teal', 'darkorange',
                            'lightseagreen'])
plt.xlabel('amino acid', size=15)
plt.ylabel('Frequency in the data', size=15)
plt.legend(frameon=False)
plt.tight_layout()
# plt.show()
plt.savefig(f'{folder_in}/Plots/AA_frequency_stacked.jpg', dpi=1200)
plt.close()

# check symmetry of tcrBLOSUMs - only C is asymmetric
path_in = '//results/tcrBlosumMtx/tcrBlosumMtx_2023-11-11'
for filename in os.listdir(path_in):
    if filename.startswith('tcrBlosum_'):
        print(filename)
        data_in = pd.read_csv(f'{path_in}/{filename}', sep='\t', index_col=0)
        # data_in.values
        # np.array_equal(data_in.values, np.transpose(data_in.values))
        # data_in[data_in.values != np.transpose(data_in.values)]
        data2 = pd.DataFrame(data_in.values != np.transpose(data_in.values),
                             columns=data_in.columns, index=data_in.index)
        print(f'Rows: {data2[data2.any()].index}')
        print([np.where(row)[0].tolist() for _, row in data2.iterrows()])
