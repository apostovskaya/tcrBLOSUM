import sys
import os
import pandas as pd

if len(sys.argv) > 3:
    print("~ Executed script: " + sys.argv[0])
    print("~ Input folder: " + sys.argv[1])
else:
    print(" No file has been provided. Files are required. ")

# input/output folder paths
folder_in = os.path.abspath(sys.argv[1])
# main_dir = os.path.dirname(os.path.realpath(__file__))
folder_out_all = f'{folder_in}/all_data'
folder_out_all_covid = f'{folder_in}/all_covid'
folder_out_all_no_covid = f'{folder_in}/all_no_covid'

for folder in [folder_out_all, folder_out_all_covid, folder_out_all_no_covid]:
    if not os.path.exists(folder):
        os.makedirs(folder)
        print(f'{folder} folder is created')

# merge all the data (.tsv) into one dataframe
all_filenames = sys.argv[2:]
tcr_data_all = pd.concat([pd.read_csv(f'{folder_in}/{filename}', sep='\t') for filename in all_filenames])

# create file with and without Covid to compare how addition/deletion of epitopes/pathogen affects tcrBlosum
tcr_data_covid = tcr_data_all.loc[tcr_data_all['Pathology'] == 'SARS-CoV-2']
tcr_data_no_covid = tcr_data_all.loc[~(tcr_data_all['Pathology'] == 'SARS-CoV-2')]

# alpha and beta
tcr_data_covid_a = tcr_data_covid.loc[tcr_data_covid['Gene'] == 'TRA']
tcr_data_covid_b = tcr_data_covid.loc[tcr_data_covid['Gene'] == 'TRB']
tcr_data_no_covid_a = tcr_data_no_covid.loc[tcr_data_no_covid['Gene'] == 'TRA']
tcr_data_no_covid_b = tcr_data_no_covid.loc[tcr_data_no_covid['Gene'] == 'TRB']

# CD4 and CD8
tcr_data_covid_cd4 = tcr_data_covid.loc[tcr_data_covid['T_cell'] == 'CD4']
tcr_data_covid_cd8 = tcr_data_covid.loc[tcr_data_covid['T_cell'] == 'CD8']
tcr_data_no_covid_cd4 = tcr_data_no_covid.loc[tcr_data_no_covid['T_cell'] == 'CD4']
tcr_data_no_covid_cd8 = tcr_data_no_covid.loc[tcr_data_no_covid['T_cell'] == 'CD8']

# all data together
# alpha and beta
tcr_data_a = tcr_data_all.loc[tcr_data_all['Gene'] == 'TRA']
tcr_data_b = tcr_data_all.loc[tcr_data_all['Gene'] == 'TRB']

# CD4 and CD8
tcr_data_cd4 = tcr_data_all.loc[tcr_data_all['T_cell'] == 'CD4']
tcr_data_cd8 = tcr_data_all.loc[tcr_data_all['T_cell'] == 'CD8']

# save in tsv
tcr_data_all.to_csv(path_or_buf=f'{folder_out_all}/all_tcrData.tsv', sep='\t', index=False)
tcr_data_a.to_csv(path_or_buf=f'{folder_out_all}/all_alpha_tcrData.tsv', sep='\t', index=False)
tcr_data_b.to_csv(path_or_buf=f'{folder_out_all}/all_beta_tcrData.tsv', sep='\t', index=False)
tcr_data_cd4.to_csv(path_or_buf=f'{folder_out_all}/all_CD4_tcrData.tsv', sep='\t', index=False)
tcr_data_cd8.to_csv(path_or_buf=f'{folder_out_all}/all_CD8_tcrData.tsv', sep='\t', index=False)

print(f'The following files have been created in the {folder_out_all}:\n',
      'all_tcrData.tsv\n',
      'all_alpha_tcrData.tsv\n', 'all_beta_tcrData.tsv\n',
      'all_CD4_tcrData.tsv\n', 'all_CD8_tcrData.tsv')

tcr_data_covid.to_csv(path_or_buf=f'{folder_out_all_covid}/covid_tcrData.tsv', sep='\t', index=False)
tcr_data_covid_a.to_csv(path_or_buf=f'{folder_out_all_covid}/covid_alpha_tcrData.tsv', sep='\t', index=False)
tcr_data_covid_b.to_csv(path_or_buf=f'{folder_out_all_covid}/covid_beta_tcrData.tsv', sep='\t', index=False)
tcr_data_covid_cd4.to_csv(path_or_buf=f'{folder_out_all_covid}/covid_CD4_tcrData.tsv', sep='\t', index=False)
tcr_data_covid_cd8.to_csv(path_or_buf=f'{folder_out_all_covid}/covid_CD8_tcrData.tsv', sep='\t', index=False)
print(f'The following files have been created in the {folder_out_all_covid}:\n',
      'covid_tcrData.tsv\n',
      'covid_alpha_tcrData.tsv\n', 'covid_beta_tcrData.tsv\n',
      'covid_CD4_tcrData.tsv\n', 'covid_CD8_tcrData.tsv')

tcr_data_no_covid.to_csv(path_or_buf=f'{folder_out_all_no_covid}/no_covid_tcrData.tsv', sep='\t', index=False)
tcr_data_no_covid_a.to_csv(path_or_buf=f'{folder_out_all_no_covid}/no_covid_alpha_tcrData.tsv', sep='\t', index=False)
tcr_data_no_covid_b.to_csv(path_or_buf=f'{folder_out_all_no_covid}/no_covid_beta_tcrData.tsv', sep='\t', index=False)
tcr_data_no_covid_cd4.to_csv(path_or_buf=f'{folder_out_all_no_covid}/no_covid_CD4_tcrData.tsv', sep='\t', index=False)
tcr_data_no_covid_cd8.to_csv(path_or_buf=f'{folder_out_all_no_covid}/no_covid_CD8_tcrData.tsv', sep='\t', index=False)
print(f'The following files have been created in the {folder_out_all_no_covid}:'
      f'\nno_covid_tcrData.tsv\n',
      'no_covid_alpha_tcrData.tsv\n', 'no_covid_beta_tcrData.tsv\n',
      'no_covid_CD4_tcrData.tsv\n', 'no_covid_CD8_tcrData.tsv')
