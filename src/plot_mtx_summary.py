import os
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns


def plot_score_sums(mtx, save_to, prefix, colname):
    # plot hist of scores
    mtx.hist(column=colname, color='tab:green')
    plt.xlabel(f'Similarity score value ({colname})')
    plt.ylabel(f'Frequency, {colname}')
    plt.title(f'Distribution of similarity score {colname}')
    plt.savefig(f'{save_to}/{prefix}_score_hist_{colname}.png')
    plt.close()

    # plot distribution of the sum to look at the edge AAs
    mtx_sort = mtx.sort_values(by=colname)
    ax = plt.gca()
    mtx_sort[[colname]].plot(kind='bar', legend=False, ax=ax, color='tab:green', rot=0)
    plt.xlabel('amino acid')
    plt.ylabel(f'Similarity scores, {colname}')
    plt.title(f'Distribution of similarity score {colname} per amino acid')
    # plt.show()
    plt.savefig(f'{save_to}/{prefix}_score_distr_{colname}.png')
    plt.close()


def plot_clustermap(matrix, matrix_handle, path_out, vmin=-10):
    sns.clustermap(matrix, col_cluster=True, cmap='seismic', annot=True,
                   robust=True,
                   vmin=vmin, vmax=10, center=0, linewidths=0.5,
                   linecolor='white')
    # plt.show()
    plt.savefig(f'{path_out}/clustermap_{matrix_handle}.jpg', dpi=1200)
    plt.close()


main_dir = os.getcwd()
folder = f'{main_dir}/results'
folder_tcr = f'{folder}/tcrBLOSUMmtx/tcrBLOSUMmtx_2024-02-06'
z_mtx = pd.read_csv(f'{folder}/zScoresBlosum_ints.tsv', sep='\t', index_col=0)
z_mtx['sums'] = z_mtx.sum(axis=1)
plot_score_sums(z_mtx, f'{folder}', 'zScores', 'sums')

# folder = '/Users/apost/Documents/CloudMail/PhD_2020/SimilarityMTXs'
tanimoto = pd.read_csv(f'{folder}/otherMtxs/aa_10xTanimotoSimilarity.tsv',
                       sep='\t', index_col=0)
tanimoto['sums'] = tanimoto.sum(axis=1)
plot_score_sums(tanimoto, f'{folder}/otherMtxs/Tanimoto', 'Tanimoto', 'sums')

for descriptor in ['PhysChemSim', 'TopoSim']:
    descr_mtx = pd.read_csv(f'{folder}/otherMtxs/{descriptor}Mtx.tsv',
                            sep='\t', index_col=0)
    if descriptor == 'TopoSim':
        plot_clustermap(descr_mtx, descriptor, f'{folder_tcr}/Plots', vmin=1)
    else:
        plot_clustermap(descr_mtx, descriptor, f'{folder_tcr}/Plots')


# the same for blosum62 for comparison
blosum62 = pd.read_csv(f'{folder_tcr}/tcrBLOSUM_BLOSUM62.tsv',
                       sep='\t', index_col=0)
blosum62['favorable'] = blosum62[blosum62 > 0].count(axis=1)
blosum62['unfavorable'] = blosum62[blosum62 < 0].count(axis=1)
blosum62['neutral'] = blosum62[blosum62 == 0].count(axis=1)
# useless plots
# plot_score_sums(blosum62, f'{folder}/otherMtxs/blosum62', 'blosum62', 'unfavorable')

# folder = '/Users/apost/Documents/CloudMail/PhD_2020/SimilarityMTXs'
chain_mtx = [('alpha', 'tcrBLOSUMa'), ('beta', 'tcrBLOSUMb')]

for chain, mtx in chain_mtx:
    tcr_blosum = pd.read_csv(f'{folder_tcr}/tcrBLOSUM_all_{chain}.tsv',
                             sep='\t', index_col=0)

    # heatmap with dendrogram on top
    sns.clustermap(tcr_blosum, col_cluster=True, cmap='seismic', annot=True,
                   robust=True,
                   vmin=-10, vmax=10, center=0, linewidths=0.5,
                   linecolor='white')
    # plt.show()
    plt.savefig(f'{folder_tcr}/Plots/clustermap_{mtx}.jpg', dpi=1200)
    plt.close()

    tcr_blosum['favorable'] = tcr_blosum[tcr_blosum > 0].count(axis=1)
    tcr_blosum['unfavorable'] = tcr_blosum[tcr_blosum < 0].count(axis=1)
    tcr_blosum['neutral'] = tcr_blosum[tcr_blosum == 0].count(axis=1)
    # useless plots
    # plot_score_sums(tcr_blosum,
    #                 f'{folder}/tcrBlosumMtx/tcrBlosumMtx_2024-02-06/Plots',
    #                 'tcrBLOSUMb', 'neutral')

    # compare the number of (un)favorable substitutions,
    # plot blosum62 and tcrBLOSUM together
    for subs_type in ['favorable', 'unfavorable', 'neutral']:
        blosum62_sort = blosum62.sort_values(by=subs_type)
        tcr_blosum_sort = tcr_blosum.loc[blosum62_sort.index.to_list(), :]
        ax = plt.gca()
        ax.set_ylim([0, 19])
        blosum62_sort[[subs_type]].plot(kind='bar', ax=ax, color='tab:blue',
                                        legend=False, position=1,
                                        width=0.2, rot=0)
        tcr_blosum_sort[[subs_type]].plot(kind='bar', ax=ax, color='tab:orange',
                                          legend=False, position=0,
                                          width=0.2, rot=0)
        plt.hlines(y=9.5, xmin=0, xmax=20, linewidth=0.4, color='black',
                   linestyles='dotted')
        plt.xlabel('amino acid')
        plt.ylabel(f'The number of {subs_type} substitutions')
        # plt.title(f'Distribution of the numbers of {subs_type} substitutions')
        # Get artists and labels for legend and chose which ones to display
        handles, labels = ax.get_legend_handles_labels()
        plt.legend(handles=handles[:2], labels=['BLOSUM62', mtx])
        plt.tight_layout()
        # plt.show()
        plt.savefig(f'{folder_tcr}/Plots/N{subs_type}_blosum62-vs-{mtx}.jpg',
                    dpi=1200)
        plt.close()
