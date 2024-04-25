import os

import pandas as pd
from scipy.spatial.distance import pdist, squareform
# https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.distance.pdist.html

from src.calculate_tcrBLOSUMs import plot_heatmap
from src.Matrix_evaluation.constants import AA_ALPHABET

# main directory
main_dir = os.getcwd()
folder_in = f'{main_dir}/data'
folder_out = f'{main_dir}/results/otherMTXs'

# z-scores from
# Peptide Quantitative Structure-Activity Relationships, Multivariate Approach
z_table = pd.read_csv(f'{folder_in}/z-scores.tsv', sep='\t', index_col=0)

# order of AA matching BLOSUM62
aa_list = [aa for aa in AA_ALPHABET]
assert all(aa in aa_list for aa in z_table.index.to_list()), \
    'Amino acids in z-scores file do not match amino acids in BLOSUM62 matrix.'
z_table = z_table.reindex(aa_list)

# a "flat" array that consists only of the upper triangle
# of the distance matrix, not including the diagonal of 0s.
z_blosum = pdist(z_table[['z1', 'z2', 'z3']].values, metric='euclidean')
# square matrix nd-array
z_blosum_square = squareform(z_blosum)
# back to DF
z_blosum_square_df = pd.DataFrame(z_blosum_square,
                                  columns=z_table.index.to_list(),
                                  index=z_table.index.to_list())
z_blosum_square_df_rounded = z_blosum_square_df.astype(int)  # or .round(decimals=1)

print('\nPhysChemSim matrix based on z-scored has been computed.')

# save calculated z-blosum matrix in tsv
# z_blosum_square_df.to_csv(f'{folder_out}/PhysChem_distMtx_raw.tsv', sep='\t')
z_blosum_square_df_rounded.to_csv(f'{folder_out}/PhysChemDistMtx.tsv',
                                  sep='\t')

# reverse values of distance mtx for similarity
z_sim_mtx = z_blosum_square_df_rounded * -1
# rescale 0 to mean self-substitution score in Blosum62 which is +6
z_sim_mtx[z_sim_mtx == 0] = 6
z_sim_mtx.to_csv(f'{folder_out}/PhysChemSimMtx.tsv', sep='\t')

# plot heatmap of similarity mtx
plot_heatmap(z_sim_mtx, f'{folder_out}/Heatmap_PhysChemSimMtx.png',
             'PhysChemSim', cent=-0.5)

# plotting heatmap of zScore dist mtx
# plot_heatmap(z_blosum_square_df_rounded, f'{folder}/Heatmap_zScoresSimMtx.png', 'zScores',
#              col_sheme='seismic_r', cent=5, min_v=0, max_v=10)
# print(f'Heatmap_zScoresSimMtx.png plot has been saved in the {folder}.')

# BLOSUM62
blosum62 = pd.read_csv(f'{folder_out}/blosum62_20aa.tsv', sep='\t', index_col=0)
# plot heatmap of blosum62
plot_heatmap(blosum62, f'{folder_out}/Heatmap_blosum62.png', 'BLOSUM62')

# VDJam
# vdjam = pd.read_csv(f'{folder_out}/vdjam.tsv', sep='\t')
