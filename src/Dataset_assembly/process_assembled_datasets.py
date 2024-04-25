import sys
import os
import pandas as pd
import copy
from datetime import date


# To process all TCR datasets for TCR-BLOSUM calculations
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# TO-DO LIST
# add comments?
# Need to delete CDR3s with other non-canonical symbols (currently, those are filtered out in individual
# preprocessing scripts: O, #, (, ), *, ?) - change RegEx?
# V, J columns are currently useless, not unified or even absent -> change them
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Right now I disregard VDJ-notations and remove all duplicated sequences, as well as kind of duplicates that occur
# only after the trimming of the end amino acids. That removes existing biological bias, so why do I do that?
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


def filter_data_for_tcr_blosum(data_raw):
    """
    To, first, remove CDR3 sequences which don't start with C and end with F or W;
    then trim C/F-or-W amino acids from the start/end;
    remove duplicated CDR3 sequences specific to the same epitope
    and originating from the same chain (alpha/beta).
    :param data_raw: df, original DataFrame
    :return: df, processed DataFrame
    """
    assert all(elem in data_raw.columns.to_list() for elem in
               ['CDR3', 'Gene', 'V', 'J', 'Epitope']), \
        "`filter_data_for_tcr_blosum`: " \
        "The input data doesn't contain required columns, please, check the " \
        "spelling and/or case. The required columns are: " \
        "'CDR3', 'Gene', 'V', 'J', 'Epitope'."

    # all CDR3s and Epitopes to capital letters
    data_raw['CDR3'] = data_raw['CDR3'].str.upper()
    data_raw['Epitope'] = data_raw['Epitope'].str.upper()
    data_raw['Gene'] = data_raw['Gene'].str.upper()

    # CHANGE V/J gene notation to IMGT without zeroes and *, /
    # remove white spaces
    data_raw['V'] = data_raw['V'].apply(str).str.replace(' ', '')
    data_raw['J'] = data_raw['J'].apply(str).str.replace(' ', '')

    # replace .
    data_raw['V'] = data_raw['V'].apply(str).str.replace('.', '')
    data_raw['J'] = data_raw['J'].apply(str).str.replace('.', '')

    # remove everything after *
    data_raw['V'] = data_raw['V'].apply(str).str.replace('\*\d+', '')
    data_raw['J'] = data_raw['J'].apply(str).str.replace('\*\d+', '')

    # remove everything after /
    data_raw['V'] = data_raw['V'].apply(str).str.replace('\/\S+', '')
    data_raw['J'] = data_raw['J'].apply(str).str.replace('\/\S+', '')

    # convert to TRBV12-3 (no zeroes) as per IMGT
    data_raw['V'] = data_raw['V'].apply(str).str.replace('V0', 'V')
    data_raw['V'] = data_raw['V'].apply(str).str.replace('-0', '-')
    data_raw['J'] = data_raw['J'].apply(str).str.replace('V0', 'V')
    data_raw['J'] = data_raw['J'].apply(str).str.replace('-0', '-')

    # FILTERING
    # remove incomplete entries
    data_raw = data_raw.dropna(subset=['CDR3', 'Gene', 'V', 'J', 'Epitope'],
                               how='any')

    # remove duplicated rows
    data = data_raw.drop_duplicates(
        subset=['CDR3', 'Gene', 'V', 'J', 'Epitope'])  # 46373

    entries_before = data.shape[0]
    print(
        f'\nThe number of non-duplicated entries before processing: '
        f'{entries_before}',
        '\nThe number of unique CDR3 entries before processing:',
        len(data['CDR3'].drop_duplicates()))  # 46266 # 37821

    # check which/how many CDRs don't start/end with C/F-or-W
    n_CFW = 0  # 37553
    n_nonCFW = 0  # 268
    nonCFW_cdrs = []
    for cdr in data['CDR3'].drop_duplicates():
        if cdr[0] == 'C' and (cdr[-1] == 'F' or cdr[-1] == 'W'):
            n_CFW += 1
        else:
            n_nonCFW += 1
            nonCFW_cdrs.append(cdr)
    nonCFW_cdrs.sort()
    print(
        f'including {n_CFW} classical C-F/W CDR3s and {n_nonCFW} non C-F/W CDR3s\n',
        'Examples of non-classical CDRs which will be removed:\n',
        nonCFW_cdrs[:15])

    # remove rows with non-classical CDRs (don't start with C and end with either F or W)
    data_filtered = data.loc[~data['CDR3'].isin(nonCFW_cdrs)].copy()

    # after:
    entries_after = data_filtered.shape[0]  # 45997
    print(
        f'The number of CDR3 entries after filtering out non-classical '
        f'(non C-F/W) CDR3s: {entries_after};\n',
        (entries_before - entries_after), 'entries removed.')  # 268
    print(
        'The number of unique CDR3 entries after filtering out non-classical '
        '(non C-F/W) CDR3s:',
        len(data_filtered['CDR3'].drop_duplicates()))  # 37552

    # trim AAs: -2…-1 for beta, -1….-1 for alpha
    # (99% of CDR3s have the same 1 or 2 AAs in these positions)
    # do -1 trimming on both ends for all CDR3s
    data_filtered['preCDR3'] = data_filtered['CDR3'].str[:1]
    data_filtered['postCDR3'] = data_filtered['CDR3'].str[-1:]
    data_filtered['CDR3'] = data_filtered['CDR3'].str[1:-1]

    # for beta CDR3s, trim 1 more AA from the start of CDR3s
    data_filtered.loc[data_filtered['Gene'] == 'TRB', 'preCDR3'] = \
        data_filtered.loc[data_filtered['Gene'] == 'TRB', 'preCDR3'].str[:] \
        + data_filtered.loc[data_filtered['Gene'] == 'TRB', 'CDR3'].str[0]
    data_filtered.loc[data_filtered['Gene'] == 'TRB', 'CDR3'] = \
        data_filtered.loc[data_filtered['Gene'] == 'TRB', 'CDR3'].str[1:]


    # remove the 1st C & the last F/W
    # (-> after that some CDRs become identical but maybe better to keep
    # them because they have different origin?)
    # data['CDR3'] = data['CDR3'].apply(lambda x: x[1:-1] if (
    #             x[0] == 'C' and (x[-1] == 'F' or x[-1] == 'W')) else x)

    # remove duplicated rows
    # data_filtered = copy.deepcopy(
    #     data.drop_duplicates(subset=['CDR3', 'Gene', 'V', 'J', 'Epitope']))
    # subset=['CDR3', 'Gene', 'Epitope']

    # after:
    # entries_after_2 = len(data_filtered['CDR3'])  # 40722
    # print(
    #     'The number of CDR3 entries after removing duplicates in CDR3 and Epitope columns:',
    #     entries_after_2, ';',
    #     (entries_after - entries_after_2), 'entries removed.')  # 5275
    print('The number of unique CDR3 entries (without duplicates '
          'in CDR3 and Epitope columns):',
          len(data_filtered['CDR3'].drop_duplicates()))  # 37524
    print('The number of unique Epitopes:',
          len(data_filtered['Epitope'].drop_duplicates()))  # 450
    return data_filtered


def create_blocks(data_filtered):
    """
    To add CDR3 length column ('CDR3_length') and divide CDR3 sequences
    into blocks, where one block is a group (>= 2) of CDR3s of the same chain
    (alpha and beta separately) and length, specific to the same epitope.
    Blocks are indicated in 'Block_N' (the number of the block) column
    of the returned dataframe.
    :param data_filtered: DataFrame with CDR3 sequences after filtering
    :return: processed DataFrame
    """
    assert all(elem in data_filtered.columns.to_list() for elem in
               ['CDR3', 'Gene', 'Epitope']), \
        "`filter_data_for_tcr_blosum`: The input data doesn't contain " \
        "required columns, please, check the spelling and/or case. " \
        "The required columns are: 'CDR3', 'Gene', 'Epitope'."

    # PROCESSING
    # add length column - added +4 but didn't rerun
    # (to have original length and not trimmed length as it is now)
    data_filtered['CDR3_length'] = data_filtered['CDR3'].apply(
        lambda x: len(x) + 4)

    # sort by values
    data_filtered.sort_values(by=['Epitope', 'Gene', 'CDR3_length', 'CDR3'],
                              inplace=True)
    data_filtered.reset_index(inplace=True, drop=True)

    # distribution of CDR3 lengths per Ep, alpha and beta chains separately
    cdr3a_lengths_perEp = {}
    cdr3b_lengths_perEp = {}
    for ep in data_filtered['Epitope'].drop_duplicates().to_list():
        # get the number of CDR3s of the same length of the same Ep for each CDR3 length
        cdr3a_lengths_perEp[ep] = data_filtered.loc[
            (data_filtered['Epitope'] == ep) &
            (data_filtered['Gene'] == 'TRA'),
            'CDR3_length'].value_counts().to_dict()
        cdr3b_lengths_perEp[ep] = data_filtered.loc[
            (data_filtered['Epitope'] == ep) &
            (data_filtered['Gene'] == 'TRB'),
            'CDR3_length'].value_counts().to_dict()
    cdr3_lengths_perEp = {'TRA': cdr3a_lengths_perEp,
                          'TRB': cdr3b_lengths_perEp}
    # structure is {'TRA': {epitope: {cdr_length: number of cdrs}}, ...}

    # CREATE BLOCKS
    # reassemble data with indicated separate alpha and beta blocks
    # for each epitope; block counts start with 1
    data_tmp = pd.DataFrame()
    block_num = 1
    for ep in data_filtered['Epitope'].drop_duplicates():  # 450 epitopes
        for gene in cdr3_lengths_perEp.keys():  # TRA or TRB
            for cdr_len in cdr3_lengths_perEp[gene][ep].keys():
                # set of lengths of alpha/beta CDR3s of one epitope
                if cdr3_lengths_perEp[gene][ep][cdr_len] >= 2:
                    data_tmp = pd.concat(
                        [data_tmp,
                         data_filtered.loc[(data_filtered['Epitope'] == ep) &
                                           (data_filtered['Gene'] == gene) &
                                           (data_filtered[
                                                'CDR3_length'] == cdr_len)]],
                        ignore_index=True, sort=True)
                    # C:\Users\User\AppData\Local\conda\conda\envs\HIV\lib\site-packages\pandas\core\frame.py:7138:
                    # FutureWarning: Sorting because non-concatenation axis is not aligned. A future version
                    # of pandas will change to not sort by default. To accept the future behavior, pass 'sort=False'.
                    # To retain the current behavior and silence the warning, pass 'sort=True'. sort=sort,
                    data_tmp.loc[(data_tmp['Epitope'] == ep) & (
                            data_tmp['Gene'] == gene) &
                                 (data_tmp['CDR3_length'] == cdr_len),
                                 'Block_N'] = block_num
                    block_num += 1
    data_blosum = copy.deepcopy(data_tmp)

    del data_tmp, ep, cdr_len, gene

    # ~~~~~~~~~~~~~~~~~ check block
    # before:
    # len(data_filtered.index)  # 40816; 37621 unique cdrs; 450 unique Eps
    # after:
    # len(data_blosum.index)  # 39564; 36550 unique cdrs; 263 unique Eps
    # to check if numbers add up
    # n_singles = 0
    # n_blocks = 0
    # for k1 in cdr3_lengths_perEp.keys():
    #     for k2 in cdr3_lengths_perEp[k1].keys():
    #         for k3 in cdr3_lengths_perEp[k1][k2].keys():
    #             if cdr3_lengths_perEp[k1][k2][k3] == 1:
    #                 n_singles += 1
    #             else:
    #                 n_blocks += 1
    # print(n_singles, n_blocks)  # 1252; 1593
    # they do: len(data_blosum.index) + n_singles = len(data_filtered.index)

    print('\nBlocks of CDR3 data were assembled. One block is a group '
          '(at least 2) of CDR3s of the same chain and length,'
          ' specific to the same epitope.\nThere are',
          int(max(data_blosum['Block_N'])), 'blocks in total, covering',
          data_blosum.shape[0], 'CDR3 sequences specific to',
          len(data_blosum.Epitope.drop_duplicates()),
          'unique epitopes;', len(data_blosum['CDR3'].drop_duplicates()),
          'CDR3s are unique.\nDuplicated CDR3 sequences, which are specific '
          'to different Eps or originate from different genes, '
          'are not removed since they were probably generated independently, '
          'thus, reflect sth about amino acid (pairs) frequencies.')
    # There are 1592 blocks in total, covering 39470 CDR3 sequences
    # specific to 263 unique epitopes; 36453 CDR3s are unique.

    # final sorting
    data_blosum.sort_values(by=['Epitope', 'Gene', 'CDR3_length', 'CDR3'],
                            inplace=True)
    data_blosum.reset_index(drop=True, inplace=True)
    return data_blosum


if __name__ == '__main__':
    main_dir = os.getcwd()
    today = date.today().strftime("%Y-%m-%d")

    folder_in = os.path.abspath(sys.argv[1])
    folder_out = f'{main_dir}/data/parsed/processed_for_tcr_blosum_calculations'
    log_folder = f'{main_dir}/logs/LogFiles_{today}/processing_logs'

    if len(sys.argv) != 2:
        print("Wrong input parameters. Expected: script_path input_folder")
    # the name of the module is:', __name__

    # input/output folder path
    # folder = os.path.abspath(sys.argv[1])
    # filename = sys.argv[2]
    # prefix = sys.argv[3]
    
    # Check if the folder exists, if not, create
    if not os.path.exists(folder_out):
        os.makedirs(folder_out)
        print(f'{folder_out} directory is created')
    if not os.path.exists(log_folder):
        os.makedirs(log_folder)
        print(f'{log_folder} directory is created')

    for filename in os.listdir(folder_in):
        if filename.endswith('.tsv'):
            prefix = '_'.join(filename.split('_')[:-1])
            # https://stackoverflow.com/questions/7152762/how-to-redirect-print-output-to-a-file
            orig_stdout = sys.stdout
            f = open(f'{log_folder}/{prefix}.txt', 'w+')
            sys.stdout = f

            print('\nThe date of the execution is:', today)
            print(f"~ Executed script: {sys.argv[0]}")
            print(f"~ Input folder: {folder_in}")
            print(f"~ Input file: {filename}")
            print(f"~ Log file: {log_folder}/{prefix}.txt")
            print(f"~ Output folder: {folder_out}")

            tcr_data_raw = pd.read_csv(f'{folder_in}/{filename}',
                                       sep='\t', low_memory=False)
    # sys:1: DtypeWarning: Columns (5,7) have mixed types.
    # Specify dtype option on import or set low_memory=False.

            # filter TCR data: remove CDR3 sequences which don't start with C and
            # end with F or W; trim 2 amino acids from the start/end;
            # remove duplicated sequences specific to the same epitope
            # and originating from the same chain (a/b) & VJ-gene.
            tcr_data_filtered = filter_data_for_tcr_blosum(tcr_data_raw)

            # process filtered TCR data: add CDR3 length column ('CDR3_length')
            # and divide CDR3 sequences into blocks, where one block is a group (>= 2)
            # of CDR3s of the same chain (alpha and beta separately) and length,
            # specific to the same epitope.
            # Blocks are indicated in 'Block_N' column of the returned dataframe
            tcr_data_blosum = create_blocks(tcr_data_filtered)

            # reorder columns
            tcr_data_blosum = tcr_data_blosum[
                ['Block_N', 'CDR3', 'CDR3_length', 'Gene', 'Epitope',
                 'MHC_class', 'T_cell', 'Pathology', 'V', 'J', 'Source']]

            # save assembled file in tsv
            tcr_data_blosum.to_csv(
                path_or_buf=f'{folder_out}/processed_TCRdata_{prefix}.tsv',
                sep='\t', index=False)
            print(f'\nThe following file has been created in the '
                  f'{folder_out}:\n', f'processed_TCRdata_{prefix}.tsv')

    sys.stdout = orig_stdout
    f.close()

    print(f'The following file has been created in the '
          f'{folder_out}:\n', f'processed_TCRdata_{prefix}.tsv')
    print(f'Log file {prefix}.txt is saved in {log_folder}')

# else:
#     print(f'Module name is: "{__name__}". '
#           f'Import of the "tcrBLOSUM_data_processing" file without execution.')
