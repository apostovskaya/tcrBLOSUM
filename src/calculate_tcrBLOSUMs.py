# TCR-BLOSUM
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# The idea is to take known, trusted pairwise/multiple CDR3 alignments similar to what we expect our next
# alignment to look like, and count the frequency at which each amino acid pair occurs,
# i.e. the frequencies of aa-pairs in closely related seqs (blocks)
# one block is a group (at least 2) of CDR3s of the same gene and length, specific to the same epitope

# The BLOSUM score for a particular substitution is a log-odds score that provides a measure of the biological
# probability of a substitution relative to the chance probability of the substitution:
# log(Prob(a1a2)/[Prob(a1) * Prob(a2)])
# Prob(a1a2) - the frequency of the substitution in functionally similar CDR3s;
# Prob(a1) and Prob(a2) are the frequencies of amino acids a1 and a2 in the database.
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# TO-DO LIST
# decide whether to include CDR3s of similar length (delta=1, if yes, how to align); currently not included
# think about asserts inside functions, checks for no CFW cdrs, for the file format
# add asserts that lists of provided amino acids match (individ freqs vs paira vs provided list)
# for pair calcs: what if there are weird aas in the DF but not in the aa_list? vice versa?
# cal_score: if aas in freq_1aa don't match freq_2aas, then what?
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
import os
import sys
from datetime import date
from math import log
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import seaborn as sns
from matplotlib.patches import Patch
from src.Matrix_evaluation.constants import AA_20


def count_aa_freq(data, all_aa=None):
    """
    calculates frequency of each amino acid (aa) in a list based on its occurrence in CDR3s of the provided data,
    excluding first C and the last F/W of each CDR3
    :param all_aa: a list/set/tuple, one letter amino acid code, the default is a tuple of standard 20 amino acids
    :param data: a Series, all CDR3 sequences with trimmed 1st C and the last F or W
    :return: two dictionaries aa_n and aa_fr where keys are 1-letter amino acid codes,
    values - amino acid counts and frequencies respectively
    """
    if all_aa is None:
        all_aa = AA_20
        # amino acid groups as in UniProt
        aa_groups = {'A': ['aliphatic', 'grey'], 'R': ['basic', 'blue'], 'N': ['amide', 'white'],
                     'D': ['acidic', 'red'],
                     'C': ['sulfur', 'yellow'], 'Q': ['amide', 'white'], 'E': ['acidic', 'red'],
                     'G': ['aliphatic', 'grey'],
                     'H': ['basic', 'blue'], 'I': ['aliphatic', 'grey'], 'L': ['aliphatic', 'grey'],
                     'K': ['basic', 'blue'], 'M': ['sulfur', 'yellow'], 'F': ['aromatic', 'black'],
                     'P': ['aliphatic', 'grey'], 'S': ['small_hydroxy', 'green'], 'T': ['small_hydroxy', 'green'],
                     'W': ['aromatic', 'black'], 'Y': ['aromatic', 'black'], 'V': ['aliphatic', 'grey']}
    assert len(all_aa) > 0, '\n`count_aa_freq`: the list of amino acids is empty. ' \
                            'You can use the default value containing all 20 AAs.'
    assert 'CDR3' in data.columns.to_list(), \
        "\n`count_aa_freq` requires 'CDR3' column, please, check the spelling and/or case in the input data."

    # aa_n = {}  # total occurrence for each individual aa in all cdr3s
    # aa_fr = {}  # frequency for each individual aa in all cdr3s
    # aa_n_sum = 0  # total occurrence of all aa-s in all cdr3s

    aa_n = {aa: data['CDR3'].str.count(aa).sum() for aa in all_aa}
    aa_n_sum = data['CDR3'].str.len().sum()
    aa_fr = {aa: round(aa_n[aa] / aa_n_sum, 4) for aa in aa_n.keys()}

    # for aa in all_aa:
    #     # n = 0
    #     # for i in range(len(data['CDR3'])):  # iteration to get occurrence of i-aa in each cdr3
    #     #     n += data['CDR3'][i].count(aa)  # total occurrence of i-aa in all cdr3s
    #     n = data['CDR3'].str.count(aa).sum()
    #     aa_n[aa] = n
    #     aa_n_sum += n
    # for aa in all_aa:
    #     aa_fr[aa] = round(aa_n[aa] / aa_n_sum, 4)

    if aa_groups:
        for aa in all_aa:
            aa_groups[aa].append(aa_fr[aa])
            aa_groups[aa].append(aa_n[aa])
        aa_df = pd.DataFrame.from_dict(aa_groups, orient='index', columns=['Property', 'Color', 'Frequency', 'Counts'])
        aa_df.sort_values(by=['Frequency', 'Property'], inplace=True)
    else:
        aa_dict = {aa: [aa_fr[aa], aa_n[aa]] for aa in aa_n.keys()}
        aa_df = pd.DataFrame.from_dict(aa_dict, orient='index', columns=['Frequency', 'Counts'])
    aa_df.reset_index(inplace=True)
    aa_df.rename(columns={'index': 'AminoAcid'}, inplace=True)

    print('\nFrequencies of individual amino acids have been calculated.')

    return aa_n, aa_fr, aa_df


def count_aa_pairs_freqs(data_with_blocks, all_aa=AA_20):
    """
    Traverses all pairwise combinations of two CDR3 sequences (strings) of the same length in every cdr3 block and
    counts amino acid pairs at the positions with the same indexes
    :param data_with_blocks: DataFrame with columns: 'CDR3' (w/o the 1st C & the last F/W), 'CDR3_length', TCR 'Gene',
    'Epitope'
    :param all_aa: a list/set/tuple, one letter amino acid code, the default is a tuple of all 20 aa
    :return: dict of dicts, {aminoAcid-1: {aminoAcid-2: freq-12, aminoAcid-3: freq-13, ..., aminoAcid-20: freq-120},
    ..., aminoAcid-20: {aminoAcid-1: freq-120, aminoAcid-2: freq-220, ..., aminoAcid-19: freq-1920}),
    freqs A-V = freqs V-A
    """
    # if all_aa is None:
    #     all_aa = ('A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V')
    assert len(all_aa) > 0, \
        '\nin `count_aa_pairs_freqs`: the list of amino acids is empty. ' \
        'You can use the default value containing all 20 AAs.'

    assert all(elem in data_with_blocks.columns.to_list() for elem in ['CDR3', 'Block_N']), \
        "\nin `count_aa_pairs_freqs`: The input data doesn't contain required " \
        "columns, please, check the spelling and/or case. " \
        "The following columns are required: 'CDR3', 'Block_N'."

    # dictionary to store frequencies of all AA pairs
    aa_pair_cnts = {aa: {aa: 0 for aa in all_aa} for aa in all_aa}
    # {aa1: {aa1: N11, aa2: N12, ...}, aa2: {aa1: N21, aa2: N22, ...}, ...}
    n_pairs = 0

    # count aa-pairs in each cdr3 block
    for block in data_with_blocks['Block_N'].drop_duplicates():
        # separate a block - all cdr3s of the same length recognizing the same Ep
        cdr3_epi_block = data_with_blocks.loc[data_with_blocks['Block_N'] == block, 'CDR3'].to_list()
        # should only contain those with >=2 (at least 2 CDR3s per block), but check just in case
        assert len(cdr3_epi_block) >= 2, \
            '\nin `count_aa_pairs_freqs`: Not enough CDR3s of the same length' \
            ' per epitope: block has < 2 sequences, `count_aa_pairs_freqs` ' \
            'requires at least 2 sequences of the same length.'
        # CDR3s from the same block should be of the same length, just in case check
        cdr3_block_lengths = list(map(lambda x: len(x), cdr3_epi_block))
        assert len(set(cdr3_block_lengths)) == 1, \
            '\nin `count_aa_pairs_freqs`: The block of CDR3s contains ' \
            'sequences of different length or is empty. ' \
            '`count_aa_pairs_freqs` requires at least 2 ' \
            'sequences of the same length.'

        for i in range(len(cdr3_epi_block) - 1):
            for j in range(i + 1, len(cdr3_epi_block)):
                # confirm that string (CDR3) lengths match
                if len(cdr3_epi_block[i]) != len(cdr3_epi_block[j]):
                    raise Exception("StringLengthError: "
                                    "`count_aa_pairs_freqs` requires two input"
                                    " sequences to be of equal length.")
                cdr3_pair_ij = zip(cdr3_epi_block[i], cdr3_epi_block[j])  # iterator (~tuple) with pairs of aa-s
                for aa1, aa2 in cdr3_pair_ij:
                    # list with alphabetically sorted aa pairs, so it's D-N pair for both, D-N and N-D
                    # aa_pair_i = sorted([aa1, aa2])
                    if aa1 == aa2:  # D-D pair
                        try:
                            aa_pair_cnts[aa1][aa2] += 1
                            n_pairs += 1
                        except KeyError:
                            print(f'KeyError: unsupported amino acid symbol '
                                  f'encountered: {aa1} or {aa2}')
                    else:
                        try:
                            aa_pair_cnts[aa1][aa2] += 1  # D-N
                            aa_pair_cnts[aa2][aa1] += 1  # N-D
                            n_pairs += 1
                        except KeyError:
                            print(f'KeyError: unsupported amino acid symbol '
                                  f'encountered: {aa1} or {aa2}')
    # convert occurrences into probabilities (frequencies)
    # (N-of-one-aa-pair/ total-N-of-all-aa-pairs)
    # aa_pair_cnts is {aa1: {aa1: N11, aa2: N12, ...},
    # aa2: {aa1: N21, aa2: N22, ...}, ...},
    # aa_pair_freqs is the same structure but all Nij-s are normalized
    # to become Frequencies of ij-pair
    aa_pair_freqs = {k1: {k2: round(v2 / n_pairs, 4) for k2, v2 in v1.items()}
                     for k1, v1 in aa_pair_cnts.items()}
    print('\nFrequencies of amino acid pairs have been calculated.')
    print(f'\nFrequencies of pairs with C are:\n{aa_pair_freqs["C"]}')

    return aa_pair_freqs


# def calc_aa_pairs_exp_freqs(pairs_freqs):
#
#     return


def calc_score(freq_1aa, freq_2aa, all_aa=AA_20, lmbd=1, log_base=2):
    """
    Computes BLOSUM-style matrix filled with substitution scores for all possible amino acid pairs
    based on the provided frequencies of individual amino acids and frequencies of amino acid pairs
    :param freq_1aa: dictionary (aa: float), where the floats are the probabilities that we expect to observe
    an aa on average in any CDR3 sequence (background frequency)
    :param freq_2aa: the target frequencies: the probability that we expect to observe residues a1 and a2 aligned in
    alignment of CDR3s with the same specificity, for all possible aa-pairs
    :param all_aa: a list/set/tuple, one letter amino acid code, the default is a tuple of all 20 aa in the same order
    as in the original BLOSUM62 matrix
    :param lmbd: int?, a scaling factor, default is set to 10 (BLOSUM62 used 0.347?)
    :param log_base: int, 2 or 10, or None (for natural log, base e); default is log2 as in BLOSUM62 (right?)
    :return: df, a calculated table (symmetrical matrix) of scores for all pairs of amino acids
    """
    # if all_aa is None:
    #     # blosum62 order without the last three B, Z, X amino acids
    #     all_aa = ('A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V')
    assert len(all_aa) > 0, \
        '\n`calc_score`: the list of amino acids is empty. ' \
        'You can use the default value containing all 20 AAs.'
    assert len(freq_1aa) == len(all_aa), \
        '\n`calc_score`: the list of amino acid frequencies is incomplete.'
    assert len(freq_2aa) == len(all_aa), \
        '\n`calc_score`: the list of amino acid pair frequencies is incomplete.'
    assert all(el1 >= 0 for el1 in freq_1aa.values()), \
        f'\n`calc_score`: provided frequencies of the following individual ' \
        f'amino acids are negative: ' \
        f'{dict((k, v) for k, v in freq_1aa.items() if v < 0)}. ' \
        f'Frequencies should be greater or equal to zero.'
    # assert all(el2 >= 0 for el2 in freq_2aa.values().values()), f'`calc_score`: the following provided frequencies ' \
    #                                                             f'of amino acid pairs are negative. ' \
    #                                                             f'Frequencies should be greater or equal to zero.'
    # d = {'a': {'a': 1, 'b': -1}, 'b': {'a': -2, 'b': 3}}
    # assert log_base, 'Invalid logarithm base is provided. Supported options are: 2, 10, None (for base e).'  # WHY?

    # substitute 0 frequencies with the smallest non-zero frequency / 100
    # to make the score calculation possible by avoiding 0 division error
    min_freq = min(i for i in freq_1aa.values() if i != 0) / 100
    for aa in all_aa:
        if freq_1aa[aa] == 0:
            # print(f"\nFrequency of the {aa} amino acid is {freq_1aa[aa]}. "
            #       f"To avoid zero division error and "
            #       f"proceed with the score calculations, "
            #       f"the frequency will be substituted with 1e-09
            #       as a 0 approximation (just a very small number).")
            # freq_1aa[aa] = 10 ** (-9)
            print(f"\nFrequency of the {aa} amino acid is {freq_1aa[aa]}. "
                  f"To avoid zero division error and proceed with "
                  f"the score calculations, the frequency will be substituted "
                  f"with {min_freq} - the smallest non-zero frequency of any "
                  f"individual amino acid divided by 100 as an approximation of 0.")
            freq_1aa[aa] = min_freq

    # find the smallest frequency of AA pair to define substitute for 0 freq
    min_obs_probs = []
    for key, val in freq_2aa.items():
        min_obs_probs.append(min(i for i in val.values() if i != 0))
    zero_prob = min(min_obs_probs) / 100

    aa_scores = {}
    for aa1 in all_aa:  # I need sth like try-except here for KeyError when not all aas are present in the data/list
        scores_aa1 = []
        for aa2 in all_aa:
            obs_prob_a1a2 = freq_2aa[aa1][aa2]
            assert obs_prob_a1a2 == freq_2aa[aa2][aa1], \
                '\n`in calc_score`: Provided frequencies of amino acid pairs' \
                ' are asymmetric.'

            if obs_prob_a1a2 == 0:
                # substitute 0 frequencies with a tiny number to avoid log0
                # print(f"\nFrequency of the {aa1}{aa2} pair is {obs_prob_a1a2}."
                #       f" To avoid taking log of 0 and  proceed with the score "
                #       f"calculations, the frequency will be substituted with "
                #       f"1e-09 as an approximation of 0 "
                #       f"(just a very small number).")
                # obs_prob_a1a2 = 10 ** (-9)
                # zero_prob = min(i for i in freq_2aa[aa1].values() if i != 0)/10
                print(f"\nFrequency of the {aa1}{aa2} pair is {obs_prob_a1a2}."
                      f" To avoid taking log of 0 and proceed with the score "
                      f"calculations, the frequency will be substituted with "
                      f"{zero_prob} - \nthe minimum non-zero frequency of "
                      f"{aa1}-pair with any amino acid "
                      f"divided by 100 as an approximation of 0.")
                obs_prob_a1a2 = zero_prob

            # compute expected probability of amino acid pair
            # based on frequencies of individual amino acids
            if aa1 == aa2:
                exp_prob_a1a2 = freq_1aa[aa1] * freq_1aa[aa2]
            else:
                # AB or BA -> double chance of occurring
                exp_prob_a1a2 = 2 * freq_1aa[aa1] * freq_1aa[aa2]

            # score calculation for different log bases
            if log_base is None:
                scores_aa1.append(int(lmbd * log(obs_prob_a1a2 / exp_prob_a1a2)))  # ln
            else:
                scores_aa1.append(int(lmbd * log(obs_prob_a1a2 / exp_prob_a1a2, log_base)))
        aa_scores[aa1] = scores_aa1

    aa_scores_df = pd.DataFrame.from_dict(aa_scores, orient='index', columns=all_aa)

    # check whether some scores are not symmetric:
    if not np.array_equal(aa_scores_df.values,
                          np.transpose(aa_scores_df.values)):
        data2 = pd.DataFrame(aa_scores_df.values != np.transpose(aa_scores_df.values),
                             columns=aa_scores_df.columns, index=aa_scores_df.index)
        print(f'Rows with asymmetric values: {data2[data2.any()].index}')
        print('For each column, print indexes of rows with asymmetric values:')
        print(f'Columns: {data2.columns.to_list()}')
        print('Rows (indexes) per column:')
        print([np.where(row)[0].tolist() for _, row in data2.iterrows()])

    sys.stdout = orig_stdout

    # resulting matrix should be symmetric
    assert np.array_equal(aa_scores_df.values,
                          np.transpose(aa_scores_df.values)), \
        '\nin `calc_score`: Constructed matrix is not symmetric! ' \
        'Check provided AA frequencies and repeat the calculation.'

    print('tcrBLOSUM matrix has been computed successfully.')

    return aa_scores_df


def change_width(ax, new_value):
    for patch in ax.patches:
        current_width = patch.get_width()
        diff = new_value - current_width
        # change the bar width
        patch.set_width(new_value)
        # recenter the bar
        patch.set_x(patch.get_x() + diff / 2)


def plot_heatmap(mtx, path, mtx_name, col_sheme='seismic',
                 cent=0, min_v=-10, max_v=10):
    # PLOTTING HEATMAP of tcrBLOSUM
    # Centering the colormap to 0 by passing the center parameter as 0.
    # Displaying the cell values with annot
    sns.heatmap(mtx, cmap=col_sheme, center=cent, annot=True, robust=True,
                vmin=min_v, vmax=max_v, linewidths=0.5, linecolor='white')
    plt.title(f"{mtx_name}", size=18)
    plt.tight_layout()
    plt.savefig(path, dpi=1200)
    plt.close()


if __name__ == '__main__':
    today = date.today().strftime("%Y-%m-%d")
    assert len(sys.argv) == 2, print("Some of the required input is missing. "
                                     "Expected format is: "
                                     "script_name input_file_name")

    # print(f'The name of the module is: {__name__}')

    # input/output folder path
    main_dir = os.getcwd()
    # folder = os.path.abspath(sys.argv[1])  # "mydir/myfile.txt"  'SimilarityMTXs'
    # filename = sys.argv[1]

    # prefix = filename.split('raw_')[1].split('_for')[0]

    folder_in = f'{main_dir}/data/{sys.argv[1]}'
    folder_out = f'{main_dir}/results'
    # folder_in = f'./data/parsed'

    # filename = filename.split('.')[0]

    # Check if the folder exists, if not, create
    if not os.path.exists(f'{main_dir}/logs/tcrBLOSUMmtx_{today}'):
        os.makedirs(f'{main_dir}/logs/tcrBLOSUMmtx_{today}')
        print(f'tcrBLOSUMmtx_{today} folder is created in logs directory')
    if not os.path.exists(f'{folder_out}/tcrBLOSUMmtx/tcrBLOSUMmtx_{today}/Plots'):
        os.makedirs(f'{folder_out}/tcrBLOSUMmtx/tcrBLOSUMmtx_{today}/Plots')
        print(f'Plots folder is created in the {folder_out}/tcrBLOSUMmtx/tcrBLOSUMmtx_{today}')
    if not os.path.exists(f'{folder_out}/tcrBLOSUMmtx/tcrBLOSUMmtx_{today}/aaSum'):
        os.makedirs(f'{folder_out}/tcrBLOSUMmtx/tcrBLOSUMmtx_{today}/aaSum')
        print(f'aaSum folder is created in the {folder_out}/tcrBLOSUMmtx/tcrBLOSUMmtx_{today}')

    for filename in os.listdir(folder_in):
        if filename.endswith('.tsv'):
            prefix = filename.split('TCRdata_')[1].split('.')[0]
            tcr_data = pd.read_csv(f'{folder_in}/{filename}', sep='\t')
            # tcr_data = pd.read_csv(f'{folder_in}/'
            #                        f'processed_for_tcr_blosum_calculations/'
            #                        f'processed_TCRdata_all_beta.tsv', sep='\t')

            orig_stdout = sys.stdout
            print("~ Prefix is: " + prefix)
            f = open(f'{main_dir}/logs/tcrBLOSUMmtx_{today}/LogFile_tcrBLOSUM_{prefix}.txt', 'w+')
            sys.stdout = f

            print('\nThe date of the execution is:', today)

            print("\n~ Executed script: " + sys.argv[0])
            print(f"~ Output folder: {folder_out}")
            print(f"~ Input file: {folder_in}/{filename}")
            print("~ Prefix is: " + prefix)

            # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            # COMPUTING
            # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            # calculate frequencies for each individual aa in a provided "list" by their occurrence in CDR3s of the tcr_data,
            # excluding first C and the last F/W of each CDR3
            aa_cnts, aa_freqs, aa_summary = count_aa_freq(tcr_data)
            sys.stdout = orig_stdout
            print(
                'Frequencies of individual amino acids have been calculated.')
            sys.stdout = f
            # .items gives a list of key-value tuples => 2nd el ([1]) is value -> sort by value
            aa_cnts_sorted = sorted(aa_cnts.items(), key=lambda kvtuple: kvtuple[1])  # logged in a file
            print(f'\nThe number of times each amino acid occurs in the CDR3 sequences in the dataset: \n{aa_cnts_sorted}')
            aa_freqs_sorted = sorted(aa_freqs.items(), key=lambda kvtuple: kvtuple[1])  # logged in a file
            print(f'\nThe frequency of each amino acid in the CDR3 sequences in the dataset: \n{aa_freqs_sorted}')

            # count frequencies of the amino acid pairs at the positions with the same indexes
            # by traversing all pairwise combinations of two CDR3 sequences (strings) of the same length in every cdr3 block
            aa_pairs_freqs = count_aa_pairs_freqs(tcr_data)
            aa_pairs_summary = pd.DataFrame.from_dict(aa_pairs_freqs)
            sys.stdout = orig_stdout
            print('Frequencies of amino acid pairs have been calculated.')
            sys.stdout = f

            # create tcr-blosum matrix for all aa-pairs
            tcr_blosum = calc_score(aa_freqs, aa_pairs_freqs, lmbd=1, log_base=2)
            sys.stdout = f
            print('\ntcrBLOSUM matrix has been computed successfully.')

            # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            # PLOTTING
            # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            # PLOTTING amino acid distribution in the dataset
            # Create bars with different colors
            ax = plt.bar(aa_summary.sort_values('Frequency', ascending=False)['AminoAcid'],
                         aa_summary.sort_values('Frequency', ascending=False)['Frequency'],
                         color=list(aa_summary.sort_values('Frequency', ascending=False)['Color']),
                         edgecolor='black')
            # increase bar width
            change_width(ax, 0.75)
            # create legend manually
            df = aa_summary[['Property', 'Color']].drop_duplicates()
            df.sort_values(by='Property', inplace=True)
            # map names to colors
            cmap = dict(zip(list(df['Property']),
                            list(df['Color'])))
            # create the rectangles for the legend
            patches = [Patch(color=v, label=k) for k, v in cmap.items()]
            # Iterate through the handles to add edge color to the legend markers
            for ha in patches:
                ha.set_edgecolor("black")

            # add the legend
            plt.legend(title='Amino acid property',
                       labels=list(df['Property']),
                       handles=patches)  # , bbox_to_anchor=(1.04, 0.5),
            # loc='center left', borderaxespad=0, fontsize=15, frameon=False)

            plt.xlabel('amino acid', size=13)
            plt.ylabel('Frequency in the data', size=13)
            plt.title(f'The frequency of each amino acid \nin the CDR3 sequences in the {prefix}', size=15)
            # plt.legend(bbox_to_anchor=(1, 1), loc=2)
            plt.tight_layout()
            plt.savefig(f'{folder_out}/tcrBLOSUMmtx/tcrBLOSUMmtx_{today}/Plots/aaDistrib_{prefix}.png')
            plt.close()

            # PLOTTING HEATMAP of tcrBLOSUM
            plot_heatmap(tcr_blosum, f'{folder_out}/tcrBLOSUMmtx/tcrBLOSUMmtx_{today}/Plots/Heatmap_{prefix}.png', prefix)
            print(f'\nHeatmap_{prefix}.png plot has been saved in the {folder_out}/tcrBLOSUMmtx/tcrBLOSUMmtx_{today}/Plots.')
            # sns.heatmap(tcr_blosum, cmap='seismic', center=0, annot=True, robust=True, vmin=-18, vmax=18,
            #             linewidths=0.5, linecolor='yellow')
            # plt.savefig(f'{folder}/tcrBlosumMtx/tcrBlosumMtx_{today}/Plots/Heatmap_{prefix}.png')
            # plt.close()

            # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            # SAVING
            # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            # SAVING summary data for individual and paired amino acids
            aa_summary.to_csv(path_or_buf=f'{folder_out}/tcrBLOSUMmtx/tcrBLOSUMmtx_{today}/aaSum/aaFreqs_{prefix}.tsv', sep='\t',
                              index=False)
            aa_pairs_summary.to_csv(path_or_buf=f'{folder_out}/tcrBLOSUMmtx/tcrBLOSUMmtx_{today}/aaSum/aaPairsFreqs_{prefix}.tsv',
                                    sep='\t', index=False)

            # SAVING calculated tcr-blosum matrix in tsv
            tcr_blosum.to_csv(path_or_buf=f'{folder_out}/tcrBLOSUMmtx/tcrBLOSUMmtx_{today}/tcrBLOSUM_{prefix}.tsv', sep='\t')
            print(f'\nThe following files have been created in the {folder_out}/tcrBLOSUMmtx/tcrBLOSUMmtx_{today}/aaSum:\n'
                  f'    aaFreqs_{prefix}.tsv, \n    aaPairsFreqs_{prefix}.tsv')
            print(f'\ntcrBLOSUM_{prefix}.tsv is saved in the {folder_out}/tcrBLOSUMmtx/tcrBLOSUMmtx_{today};\n'
                  f'\nLogFile_tcrBLOSUM_{prefix}.txt is saved in the {main_dir}/logs/tcrBLOSUMmtx_{today}')

            sys.stdout = orig_stdout
            f.close()

else:
    print('Module name is: "' + __name__ + '". Import of the "tcrBLOSUM" functions without execution.')

# https://stackoverflow.com/questions/37958360/comparing-characters-in-strings/37958526

# import sys
# if not sys.version_info.major == 3 and sys.version_info.minor >= 7:
#     print("Python 3.7 or higher is required.")
#     print("You are using Python {}.{}.".format(sys.version_info.major, sys.version_info.minor))
#     sys.exit(1)
