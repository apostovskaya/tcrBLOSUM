import sys
import os
import pandas as pd

# To parse VDJdb data to format with selected columns needed for TCR-BLOSUM calculations
# numbers in the comments are for this file: VDJdb_Human_TRA_TRB_2020-05-26
# to run: python /Users/apost/Documents/CloudMail/PhD_2020/tcrBlosum/src/Dataset_assembly/VDJdb_parsing.py 'data' 'VDJdb_ab_2023-06-13.tsv'

if len(sys.argv) > 2:
    print("~ Script: " + sys.argv[0])
    print("~ Folder: " + sys.argv[1])
    print("~ Filename: " + sys.argv[2])
else:
    print(" No file has been provided. File is required. ")

# input/output folder path
folder = os.path.abspath(sys.argv[1])  # "mydir/myfile.txt"
filename = sys.argv[2]
# folder = 'Z:/PhD_2020/SimilarityMTXs/'
# filename = 'VDJdb_Human_TRA_TRB_2021-10-04'

# read in alpha/beta cdr3 data -> pandas DF
VDJdb = pd.read_csv(f'{folder}/raw/{filename}', sep='\t')
if len(VDJdb['Species'].value_counts().to_dict()) > 1:  # to keep only human data
    VDJdb = VDJdb.loc[[VDJdb['Species'] == 'HomoSapiens']]
print('Initial data (species, counts):\n', VDJdb['Species'].value_counts().to_dict(), '\n')  # 'HomoSapiens': 69519

# filter out 10x suspicious data (cdrs specific to 1 epitope)
VDJdb = VDJdb.loc[~(VDJdb["Epitope"] == "KLGGALQAK") | ~(VDJdb["Reference"] == \
                 "https://www.10xgenomics.com/resources/application-notes/a-new-way-of-exploring-immunity-linking"
                 "-highly-multiplexed-antigen-recognition-to-immune-repertoire-and-phenotype/#")]  # 41575

# remove incomplete entries
VDJdb = VDJdb.dropna(subset=['CDR3', 'Gene', 'V', 'J', 'Epitope'], how='any')  # 41325

# keep only columns from the VDJdb that I might need
VDJdb = VDJdb[['CDR3', 'Gene', 'V', 'J', 'Epitope', 'MHC class', 'Epitope species']]

# rename columns
VDJdb.rename(columns={'MHC class': 'MHC_class', 'Epitope species': 'Pathology'}, inplace=True, errors="raise")

# indicate T cell type where MHC_class is known
VDJdb.loc[VDJdb['MHC_class'] == 'MHCI', 'T_cell'] = 'CD8'
VDJdb.loc[VDJdb['MHC_class'] == 'MHCII', 'T_cell'] = 'CD4'

# remove duplicated rows
VDJdb = VDJdb.drop_duplicates(subset=['CDR3', 'Gene', 'V', 'J', 'Epitope'])  # 41325

# add Source column
VDJdb['Source'] = 'VDJdb'

print('Filtered data (number of CDR3 entries):\n', len(VDJdb.index), '\n')  # 32853
print('Filtered data (number of unique CDR3-Epitope entries):\n',
      len(VDJdb.drop_duplicates(subset=['CDR3', 'Epitope']).index), '\n')  # 32412

# save parsed file in tsv
VDJdb.to_csv(path_or_buf=f'{folder}/parsed/parsed_{filename}', sep='\t', index=False)
print(f'The following file has been created in the {folder}/parsed:\n',
      f'parsed_{filename}')
