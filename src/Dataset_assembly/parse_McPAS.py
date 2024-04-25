import sys
import os
import pandas as pd

# To parse McPAS data to VDJdb format with selected columns needed
# for TCR-BLOSUM calculations
# numbers in the comments are for this file: McPAS_2020-06-18
# to run from terminal:
# python /Users/apost/Documents/CloudMail/PhD_2020/tcrBlosum/src/McPAS_parsing.py 'data' 'McPAS_2022-09-10.csv'

if len(sys.argv) > 2:
    print("~ Script: " + sys.argv[0])
    print("~ Folder: " + sys.argv[1])
    print("~ Filename: " + sys.argv[2])
else:
    print(" No file has been provided. File is required. ")

# input/output folder path
folder = os.path.abspath(sys.argv[1])  # "mydir/myfile.txt"
filename = sys.argv[2]
# folder = '/Users/apost/Documents/CloudMail/PhD_2020/SimilarityMTXs/'
# folder = '/Users/apost/Documents/CloudMail/PhD_2020/tcrBlosum/data'
# filename = 'McPAS_2021_10_04'
# filename = 'McPAS_2022-09-10'

# read in data -> pandas DF
mcpas = pd.read_csv(f'{folder}/raw/{filename}', sep=',',
                    encoding='latin1', low_memory=False)  # 21689
filename = filename.split('.')[0]
# DtypeWarning: Columns (13,17,24,26,28) have mixed types.
# Specify dtype option on import or set low_memory=False.
print('Initial data (species, counts):\n',
      mcpas['Species'].value_counts().to_dict())
# Human 17983;  Mouse 3662
# {'Human': 36219, 'Mouse': 3722}

# keep only human data (remove mouse)
mcpas = mcpas.loc[mcpas['Species'] == 'Human']
print('Only human data is kept (species, counts):\n',
      mcpas['Species'].value_counts().to_dict())  # Human 17983

# remove incomplete epitope entries
mcpas = mcpas.dropna(subset=['Epitope.peptide'], how='any')
# len(mcpas.index)  # 13070 -> 13701
# mcpas[mcpas['MHC'] == 'H-2Kb']
# Out[26]:
#      CDR3.alpha.aa CDR3.beta.aa Species  ... Mouse.strain   PubMed.ID Remarks
# 7440  CAVRNNNARLMF  CASSVVNEAFF   Human  ...      C57/BL6  28636592.0     NaN
# remove rows with provided Mouse.strain
mcpas = mcpas[mcpas['Mouse.strain'].isna()]
# 13566

# keep only useful columns
mcpas = mcpas[['CDR3.alpha.aa', 'CDR3.beta.aa', 'Epitope.peptide',
               'MHC', 'TRAV', 'TRAJ', 'TRBV', 'TRBJ',
               'Pathology', 'T.Cell.Type']]

# Introduce MHC_class
mcpas['MHC'] = mcpas['MHC'].str.upper()
mcpas.loc[mcpas['MHC'].str.contains(pat="^HLA-[ABC]", na=False),
          'MHC_class'] = 'MHCI'
mcpas.loc[mcpas['MHC'].str.contains(pat="^HLA-D|^D", na=False),
          'MHC_class'] = 'MHCII'

# indicate T cell type where MHC_class is known
# be aware, some MHC-II epitopes have double CD4, CD8 cell assignment for some reason
# mcpas.loc[mcpas['MHC_class'] == 'MHCII', 'T.Cell.Type'].unique()
# Out[11]: array([nan, 'CD4', 'CD4,CD8'], dtype=object)
# mcpas.loc[mcpas['T.Cell.Type'] == 'CD4,CD8', 'MHC_class'].unique()
# array(['MHCII'], dtype=object)
# 558 rows, all specific to 1 HIV epitope 'RFYKTLRAEQASQ', will consider them CD4
mcpas.loc[mcpas['MHC_class'] == 'MHCI', 'T.Cell.Type'] = 'CD8'
mcpas.loc[mcpas['MHC_class'] == 'MHCII', 'T.Cell.Type'] = 'CD4'

# vice versa, assign MHC_class based on T cell type
mcpas.loc[mcpas['T.Cell.Type'] == 'CD8', 'MHC_class'] = 'MHCI'
mcpas.loc[mcpas['T.Cell.Type'] == 'CD4', 'MHC_class'] = 'MHCII'

# split into alpha and beta df
mcpas_alpha = mcpas.dropna(subset=['CDR3.alpha.aa', 'Epitope.peptide'],
                           how='any')
alpha_cdrs = mcpas_alpha.shape[0]  # 3732
mcpas_beta = mcpas.dropna(subset=['CDR3.beta.aa', 'Epitope.peptide'],
                          how='any')
beta_cdrs = mcpas_beta.shape[0]  # 11625
entries_before = alpha_cdrs + beta_cdrs

# removing incomplete CDR3 or Epitope entries
mcpas_alpha = mcpas_alpha[
    ['CDR3.alpha.aa', 'Epitope.peptide', 'MHC_class', 'Pathology',
     'T.Cell.Type', 'TRAV', 'TRAJ']]
mcpas_beta = mcpas_beta[
    ['CDR3.beta.aa', 'Epitope.peptide', 'MHC_class', 'Pathology',
     'T.Cell.Type', 'TRBV', 'TRBJ']]

# reformat columns
#   alpha
mcpas_alpha.rename(
    columns={'CDR3.alpha.aa': 'CDR3', 'Epitope.peptide': 'Epitope',
             'TRAV': 'V', 'TRAJ': 'J', 'T.Cell.Type': 'T_cell'}, inplace=True,
    errors="raise")
mcpas_alpha['Gene'] = 'TRA'
mcpas_alpha = mcpas_alpha[
    ['CDR3', 'Gene', 'V', 'J', 'Epitope', 'MHC_class', 'Pathology', 'T_cell']]
#   beta
mcpas_beta.rename(
    columns={'CDR3.beta.aa': 'CDR3', 'Epitope.peptide': 'Epitope',
             'TRBV': 'V', 'TRBJ': 'J', 'T.Cell.Type': 'T_cell'}, inplace=True,
    errors="raise")
mcpas_beta['Gene'] = 'TRB'
mcpas_beta = mcpas_beta[
    ['CDR3', 'Gene', 'V', 'J', 'Epitope', 'MHC_class', 'Pathology', 'T_cell']]
print('Columns were renamed and reordered to VDJdb format.')

# all CDR3s and Epitopes to capital letters
#   alpha
mcpas_alpha['CDR3'] = mcpas_alpha['CDR3'].str.upper()
mcpas_alpha['Epitope'] = mcpas_alpha['Epitope'].str.upper()
#   beta
mcpas_beta['CDR3'] = mcpas_beta['CDR3'].str.upper()
mcpas_beta['Epitope'] = mcpas_beta['Epitope'].str.upper()

# remove non-canonical symbols (O, *, (, ), #, ?)
# from the beginning of CDR3s
mcpas_alpha['CDR3'] = mcpas_alpha['CDR3'].str.replace(
    pat="^O|^\*|^\(|^\)|^\#|^\?", repl='')
mcpas_beta['CDR3'] = mcpas_beta['CDR3'].str.replace(
    pat="^O|^\*|^\(|^\)|^\#^\?", repl='')
# from the end of CDR3s
mcpas_alpha['CDR3'] = mcpas_alpha['CDR3'].str.replace(
    pat="O$|\*$|\($|\)$|\#$|\?$", repl='')
mcpas_beta['CDR3'] = mcpas_beta['CDR3'].str.replace(
    pat="O$|\*$|\($|\)$|\#$|\?$", repl='')

# remove rows with non-canonical CDR3 amino acid notation
mcpas_alpha = mcpas_alpha.loc[
    ~mcpas_alpha['CDR3'].str.contains("O|\*|\(|\)|\#|\?")]  # 3728
alpha_cdrs_cut = alpha_cdrs - mcpas_alpha.shape[0]
mcpas_beta = mcpas_beta.loc[
    ~mcpas_beta['CDR3'].str.contains("O|\*|\(|\)|\#|\?")]  # 11526
beta_cdrs_cut = beta_cdrs - mcpas_beta.shape[0]
print(alpha_cdrs_cut, 'alpha and', beta_cdrs_cut,
      'beta CDR3s with non-classical symbols (O, *, (, ), #, ?) were removed.')

# drop duplicates
mcpas_alpha.drop_duplicates(subset=['CDR3', 'V', 'J', 'Epitope'],
                            inplace=True)  # 2765
# len(mcpas_alpha.drop_duplicates(subset=['CDR3', 'Epitope']).index)  # 2658
mcpas_beta.drop_duplicates(subset=['CDR3', 'V', 'J', 'Epitope'],
                           inplace=True)  # 10354
# len(mcpas_beta.drop_duplicates(subset=['CDR3', 'Epitope']).index)  # 10092

# merge alpha and beta into parsed mcpas
mcpas_parsed = pd.concat([mcpas_alpha, mcpas_beta], ignore_index=True)
entries_cut = entries_before - len(mcpas_parsed.index)
print(mcpas_parsed.shape[0], 'entries remained after processing,', entries_cut,
      'entries were filtered out.')
# 13119, 2238
del mcpas_alpha, mcpas_beta

# add Source column
mcpas_parsed['Source'] = 'McPAS'

# save parsed file in tsv
mcpas_parsed.to_csv(path_or_buf=f'{folder}/parsed/parsed_{filename}.tsv',
                    sep='\t', index=False)
print(f'The following file has been created in the {folder}/parsed:\n',
      f'parsed_{filename}.tsv')
