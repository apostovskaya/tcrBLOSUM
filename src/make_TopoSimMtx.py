# might be needed for rdkit
# import sys
# sys.path.append('/usr/local/lib/python3.7/site-packages/')

import os

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

from rdkit.Chem import AllChem
from rdkit.Chem import PandasTools
from rdkit import DataStructs
# from rdkit import Chem
# from rdkit.Chem import Draw
# from rdkit.Chem.Draw import IPythonConsole
# from rdkit.Chem import Descriptors

from src.calculate_tcrBLOSUMs import plot_heatmap
from src.Matrix_evaluation.constants import AA_ALPHABET

# main directory
main_dir = os.getcwd()
folder_in = f'{main_dir}/data'
folder_out = f'{main_dir}/results/otherMTXs'

# a file with aa info from this source:
# https://binfalse.de/software/parseable-biodata/amino-acids/,
# except smiles, they were substituted from Wikipedia
# since RDkit couldn't convert initial smiles to molecule -> NONE
aa_raw = pd.read_csv(f'{folder_in}/aminoacids_raw_forTopoSim.tsv',
                     sep='\t', index_col=False,
                     names=['Name', '3-Letter-Code', '1-Letter-Code',
                            'Molecular_formula', 'SMILES', 'polar/nonpolar',
                            'neutral/pos/neg', 'Hydropathy_index', 'pI',
                            'pK1(-COOH)', 'pK2(-NH3)'])

# add a column with a molecule encoding,
# can be used for drawing or to create a fingerprint
PandasTools.AddMoleculeColumnToFrame(aa_raw, 'SMILES', 'Molecule_picture')
# img = Draw.MolsToGridImage(aa_raw['Molecule_picture'], molsPerRow=5)

# create Morgan fingerprints (fp) (ECFP, extended connectivity fp) for all AAs
# variable to store fps for each aa
aa_fingerprints = {}
# for each aa - dic of 3 fp (fp of the default type, non-readable;
# converted to numpy array; indexes of non-zero bits)

for i in range(len(aa_raw.index)):
    aa = aa_raw.iloc[i, aa_raw.columns.get_loc("1-Letter-Code")]
    aa_mol = aa_raw.iloc[i, aa_raw.columns.get_loc("Molecule_picture")]
    bi = {}
    # 1024 bits, radius=2 (2 bonds)
    aa_fp = AllChem.GetMorganFingerprintAsBitVect(aa_mol, 2,
                                                  nBits=1024, bitInfo=bi)
    aa_fp_arr = np.zeros((1,))
    # get np array of fp instead of bits (= list of 0s and 1s)
    DataStructs.ConvertToNumpyArray(aa_fp, aa_fp_arr)
    # indexes of non-0 elements in the np array same as np.nonzero(aa_fp_arr)
    aa_fp_nz = set(aa_fp.GetOnBits())
    prints = [(aa, x, bi) for x in aa_fp_nz]
    # drawing for IPython
    # Draw.DrawMorganBits(prints, molsPerRow=4,
    #                     legends=[str(x) for x in aa_fp_nz])
    aa_fingerprints[aa] = {'fp': aa_fp, 'fp_arr': aa_fp_arr,
                           'fp_arr_nz': aa_fp_nz}

# ~~~~~~~~~~ TEST BLOCK ~~~~~~~~~~
# hand-made Tanimoto
common = aa_fingerprints['A']['fp_arr_nz'] & aa_fingerprints['G']['fp_arr_nz']
combined = aa_fingerprints['A']['fp_arr_nz'] | aa_fingerprints['G']['fp_arr_nz']
print(len(common)/len(combined))  # 0.29411764705882354

# in-build Tanimoto
# 0.29411764705882354
DataStructs.FingerprintSimilarity(aa_fingerprints['A']['fp'],
                                  aa_fingerprints['G']['fp'])
# the same
DataStructs.TanimotoSimilarity(aa_fingerprints['A']['fp'],
                               aa_fingerprints['G']['fp'])
# ~~~~~~~~~~ END of TEST BLOCK ~~~~~~~~~~

aa_list = [aa for aa in AA_ALPHABET]

aa_fp_sim = {}
for aa_1 in aa_list:
    # calculate fp similarity (TanimotoSimilarity is the default,
    # same as DataStructs.TanimotoSimilarity)
    aa1_fp_sim = [DataStructs.FingerprintSimilarity(
        aa_fingerprints[aa_1]['fp'],
        aa_fingerprints[aa_2]['fp']) for aa_2 in aa_list]
    aa1_fp_sim_scaled = list(map(lambda x: round(10*x), aa1_fp_sim))
    aa_fp_sim[aa_1] = aa1_fp_sim_scaled
    print(aa1_fp_sim_scaled)

aa_ScaledTanimoto_mtx = pd.DataFrame.from_dict(aa_fp_sim,
                                               orient='index', columns=aa_list)
aa_ScaledTanimoto_mtx.to_csv(f'{folder_out}/TopoSimMtx.tsv', sep='\t')

# PLOT HEATMAP
tanimoto = pd.read_csv(f'{folder_out}/TopoSimMtx.tsv', sep='\t', index_col=0)
# Centering the colormap to 0 by passing the center parameter as 0.
# Displaying the cell values
plot_heatmap(tanimoto, f'{folder_out}/Heatmap_TopoSim.png', 'TopoSim',
             col_sheme='seismic', cent=5, min_v=1, max_v=10)
