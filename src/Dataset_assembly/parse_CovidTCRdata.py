import pandas as pd
import sys
import os

# To parse Pieter Meysman (PMe) TCR data to VDJdb format with selected columns needed for TCR-BLOSUM calculations
# small dataset from Sofie and immuneCode Covid-data, alterations are in comments
# to run:
#


if len(sys.argv) > 3:
    print("~ Script: " + sys.argv[0])
    print("~ Folder: " + sys.argv[1])
    print("~ Filename1: " + sys.argv[2])
    print("~ Filename2: " + sys.argv[3])
else:
    print(" No file has been provided. File is required. ")

# input/output folder path
folder = os.path.abspath(sys.argv[1])  # "mydir/myfile.txt"
filenames = [sys.argv[2], sys.argv[3]]
# folder = 'Z:/PhD_2020/SimilarityMTXs/'
# filename = 'otherTCRdata.csv'  # PMe data from Sofie csv, used to be PMe_TCRdata_parsed-by-Sofie
# filename = 'Covid_iCode.txt'  # PMe Covid data from immuneCode txt

for filename in filenames:
    if filename.endswith('csv'):
        # read in data -> pandas DF, manual alternatives for 2 datasets
        tcr_data = pd.read_csv(f'{folder}/raw/{filename}', sep=',')
        # 296 PMe data from Sofie
        # keep only useful columns
        # not applicable for # PMe Covid data from immuneCode
        # (HLA_peptide is called Epitope)
        tcr_data = tcr_data[
            ['TRBV_gene', 'CDR3_beta', 'TRBJ_gene', 'HLA_peptide',
             'Epitope species']]

        # rename & reorder columns
        tcr_data = tcr_data.rename(
            columns={'CDR3_beta': 'CDR3', 'HLA_peptide': 'Epitope',
                     'TRBV_gene': 'V', 'TRBJ_gene': 'J',
                     'Epitope species': 'Pathology'}, errors="raise")
    else:
        # PMe Covid data from immuneCode; 27828 CDR3 sequences.
        tcr_data = pd.read_csv(f'{folder}/raw/{filename}', sep='\t')
        tcr_data = tcr_data.rename(
            columns={'CDR3_beta': 'CDR3', 'TRBV_gene': 'V', 'TRBJ_gene': 'J'},
            errors="raise")
        tcr_data['Pathology'] = 'SARS-CoV-2'

    print('\nInitial data:', len(tcr_data.index), 'CDR3 sequences.')

    # everything else is applicable to both
    filename = filename.split('.')[0]
    tcr_data['Gene'] = 'TRB'
    tcr_data['MHC_class'] = 'nan'
    tcr_data['T_cell'] = 'nan'
    tcr_data = tcr_data[['CDR3', 'Gene', 'V', 'J', 'Epitope', 'MHC_class', 'Pathology', 'T_cell']]
    print('Columns were renamed and reordered to VDJdb format.')

    # remove incomplete cdr3 or epitope entries
    tcr_data.dropna(subset=['CDR3', 'Epitope'], inplace=True)  # 296, no nan

    # all CDR3s and Epitopes to capital letters
    tcr_data['CDR3'] = tcr_data['CDR3'].str.upper()
    tcr_data['Epitope'] = tcr_data['Epitope'].str.upper()

    # remove non-canonical symbols (O, *, (, ), #, ?)
    # from the beginning of CDR3s
    tcr_data['CDR3'] = tcr_data['CDR3'].str.replace(pat="^O|^\*|^\(|^\)|^\#|^\?", repl='')
    # from the end of CDR3s
    tcr_data['CDR3'] = tcr_data['CDR3'].str.replace(pat="O$|\*$|\($|\)$|\#$|\?$", repl='')

    # remove rows with non-canonical CDR3 amino acid notation
    tcr_data = tcr_data.loc[~tcr_data['CDR3'].str.contains("O|\*|\(|\)|\#|\?")]  # 0 removed

    # remove duplicates
    tcr_data.drop_duplicates(subset=['CDR3', 'V', 'J', 'Epitope'], inplace=True)  # 296
    # tcr_data.drop_duplicates(subset=['CDR3', 'Epitope'], inplace=True)  # 289

    print('Remaining data after filtering:', len(tcr_data.index), 'CDR3 sequences.')

    # add Source column
    tcr_data['Source'] = 'ADREM'
    # tcr_data['Source'] = 'immuneCode'

    # save parsed file in tsv
    tcr_data.to_csv(path_or_buf=f'{folder}/parsed/parsed_{filename}.tsv', sep='\t', index=False)
    print(f'The following file has been created in the {folder}/parsed:\n', f'parsed_{filename}.tsv')
