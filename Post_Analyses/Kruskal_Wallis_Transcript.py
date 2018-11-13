'''
Description: See Title
Author: Alexander Brown
Date: 1/15/2018
'''

import csv
import os
from os import listdir
from os.path import isfile, join
from scipy import stats

import numpy as np
from itertools import combinations
from statsmodels.stats.multitest import multipletests
from statsmodels.stats.libqsturng import psturng
import warnings

def import_BFL():
    import sys
    my_dir = "I:\Function_Libraries"
    sys.path.insert(1, my_dir)
    import Brown_Function_Library as BFL
    global BFL
    print("BFL Imported" + '\n')
    
def get_file_list(directory):
    file_list = [directory + '/' + f for f in listdir(directory) if isfile(join(directory, f))]
    return file_list

def get_combined_list(file_list, normalize):
    combined_list = []
    for file in file_list:
        reader = csv.reader(open(file, "r"), delimiter = '\t')
        header = next(reader)
        observed_count_list = []
        expected_count_list = []
        y_list = []
        for line in reader:
            o_count = int(line[1])
            e_count = float(line[2])
            try:
                enrichment = o_count/e_count
            except ZeroDivisionError:
                enrichment = 0
            observed_count_list.append(o_count)
            expected_count_list.append(e_count)
            y_list.append(enrichment)

        del y_list[0]
        del y_list[-1]
        del observed_count_list[0]
        del observed_count_list[-1]
        del expected_count_list[0]
        del expected_count_list[-1]

        if normalize == True:
            enrichment_normalization = sum(expected_count_list)/sum(observed_count_list)
            y_list = [y*enrichment_normalization for y in y_list]
        combined_list.extend(y_list)
    return combined_list

def make_range_list(min_value, max_value, bin_total):
    range_list = []
    value_difference = max_value - min_value
    bin_interval = value_difference/bin_total
    bin_start = min_value
    bin_end = bin_start + bin_interval
    bin_num = 1
    while bin_num < bin_total +1:
        bin_list = [bin_start, bin_end]
        range_list.append(bin_list)
        bin_start = bin_end
        bin_end = bin_start + bin_interval
        bin_num += 1
    return range_list

def find_nth_occurrence_rindex(search_string, find_string, nth):
    start = 0
    end = len(search_string)
    while nth > 0:
        index = search_string.rindex(find_string, start, end)
        end = index - 1
        nth -= 1
    return index

def make_tally_dict(file_ID_list, range_list):
    tally_dict = {}

    for bin_list in range_list:
        tally_key = '-'.join([str(s) for s in bin_list])
        tally_dict[tally_key] = {}
        for file_ID in file_ID_list:
            tally_dict[tally_key][file_ID] = 0
    return tally_dict

def bisect_search_version_between_values(query, search_list):
    list_len = len(search_list)
    low_bound = 0
    high_bound = list_len-1
    guess = round((high_bound + low_bound)/2)
    exist = 1

    start_value = search_list[guess][0]
    end_value = search_list[guess][1]

    pg_ticker = 0
    while not start_value < query < end_value:
        previous_guess = guess

        if end_value < query:
            low_bound = guess
            guess = round((high_bound + low_bound)/2)
        if start_value > query:
            high_bound = guess
            guess = round((high_bound + low_bound)/2)
        if guess == previous_guess:
            if pg_ticker == 1:
                exist = 0
                break
            if guess == list_len -2:
                guess = guess + 1
                pg_ticker = 1
            else:
                exist = 0
                break

        start_value = search_list[guess][0]
        end_value = search_list[guess][1]
    return guess, exist

def tally_banana(file, file_ID, range_list, tally_dict, normalize):
    x_list = []
    y_list = []
    observed_count_list = []
    expected_count_list = []

    reader = csv.reader(open(file, "r"), delimiter = '\t')
    header = next(reader)
    nuc_count = header[-1]
    for line in reader:
        bp = int(line[0]) -75
        o_count = int(line[1])
        e_count = float(line[2])
        try:
            enrichment = o_count/e_count
        except ZeroDivisionError:
            enrichment = 0
        x_list.append(bp)
        observed_count_list.append(o_count)
        expected_count_list.append(e_count)
        y_list.append(enrichment)

    del x_list[0]
    del x_list[-1]
    del y_list[0]
    del y_list[-1]
    del observed_count_list[0]
    del observed_count_list[-1]
    del expected_count_list[0]
    del expected_count_list[-1]

    if normalize == True:
        enrichment_normalization = sum(expected_count_list)/sum(observed_count_list)
        y_list = [y*enrichment_normalization for y in y_list]

    for y_value in y_list:
        guess, exist = bisect_search_version_between_values(y_value, range_list)
        if exist == 0:
            continue
        bin_list = range_list[guess]
        tally_key = '-'.join([str(s) for s in bin_list])
        tally_dict[tally_key][file_ID] += 1

def kw_dunn(groups, to_compare, alpha, method):
    """ Modified From: https://gist.github.com/alimuldal/fbb19b73fa25423f02e8

    Kruskal-Wallis 1-way ANOVA with Dunn's multiple comparison test

    Arguments:
    ---------------
    groups: sequence
        arrays corresponding to k mutually independent samples from
        continuous populations

    to_compare: sequence
        tuples specifying the indices of pairs of groups to compare, e.g.
        [(0, 1), (0, 2)] would compare group 0 with 1 & 2. by default, all
        possible pairwise comparisons between groups are performed.

    alpha: float
        family-wise error rate used for correcting for multiple comparisons
        (see statsmodels.stats.multitest.multipletests for details)

    method: string
        method used to adjust p-values to account for multiple corrections (see
        statsmodels.stats.multitest.multipletests for options)

    Returns:
    ---------------
    H: float
        Kruskal-Wallis H-statistic

    p_omnibus: float
        p-value corresponding to the global null hypothesis that the medians of
        the groups are all equal

    Z_pairs: float array
        Z-scores computed for the absolute difference in mean ranks for each
        pairwise comparison

    p_corrected: float array
        corrected p-values for each pairwise comparison, corresponding to the
        null hypothesis that the pair of groups has equal medians. note that
        these are only meaningful if the global null hypothesis is rejected.

    reject: bool array
        True for pairs where the null hypothesis can be rejected for the given
        alpha

    Reference:
    ---------------
    Gibbons, J. D., & Chakraborti, S. (2011). Nonparametric Statistical
    Inference (5th ed., pp. 353-357). Boca Raton, FL: Chapman & Hall.

    """

    # omnibus test (K-W ANOVA)
    # -------------------------------------------------------------------------

    groups = [np.array(gg) for gg in groups]

    k = len(groups)

    n = np.array([len(gg) for gg in groups])
    if np.any(n < 5):
        warnings.warn("Sample sizes < 5 are not recommended (K-W test assumes "
                      "a chi square distribution)")

    allgroups = np.concatenate(groups)
    N = len(allgroups)
    ranked = stats.rankdata(allgroups)

    # correction factor for ties
    T = stats.tiecorrect(ranked)
    if T == 0:
        raise ValueError('All numbers are identical in kruskal')

    # sum of ranks for each group
    j = np.insert(np.cumsum(n), 0, 0)
    R = np.empty(k, dtype=np.float)
    for ii in range(k):
        R[ii] = ranked[j[ii]:j[ii + 1]].sum()

    # the Kruskal-Wallis H-statistic
    H = (12. / (N * (N + 1.))) * ((R ** 2.) / n).sum() - 3 * (N + 1)

    # apply correction factor for ties
    H /= T

    df_omnibus = k - 1
    p_omnibus = stats.chisqprob(H, df_omnibus)

    # multiple comparisons
    # -------------------------------------------------------------------------

    # by default we compare every possible pair of groups
    if to_compare is None:
        to_compare = tuple(combinations(range(k), 2))

    ncomp = len(to_compare)

    Z_pairs = np.empty(ncomp, dtype=np.float)
    p_uncorrected = np.empty(ncomp, dtype=np.float)
    Rmean = R / n

    for pp, (ii, jj) in enumerate(to_compare):
        # standardized score
        ts3_ts = list(np.unique(allgroups, return_counts=True)[1])
        E_ts3_ts = sum([x**3 - x for x in ts3_ts if x>1])

        if sum([x>1 for x in ts3_ts]) > 0:
            warnings.warn("We see ties.")

            yi = np.abs(Rmean[ii] - Rmean[jj])
            theta10 = (N * (N + 1)) / 12
            theta11 =  E_ts3_ts / ( 12* (N - 1) )
            theta2 = (1 / n[ii] + 1 / n[jj])
            theta = np.sqrt( (theta10 - theta11) * theta2 )
            Zij = yi / theta
        else:
            Zij = (np.abs(Rmean[ii] - Rmean[jj]) /
                   np.sqrt((1. / 12.) * N * (N + 1) * (1. / n[ii] + 1. / n[jj])))

        Z_pairs[pp] = Zij

    # corresponding p-values obtained from upper quantiles of the standard
    # normal distribution
    p_uncorrected = stats.norm.sf(Z_pairs) * 2.

    # correction for multiple comparisons
    reject, p_corrected, alphac_sidak, alphac_bonf = multipletests(
        p_uncorrected, alpha=alpha, method=method
    )

    return H, p_omnibus, Z_pairs, p_corrected, reject

def prepare_out_new_line(new_line):
    new_line = '\t'.join([str(s) for s in new_line])
    return new_line

def saddle_up():
    histone_tally_dir = r"I:\Melanoma_2018_Nuc_Project\Transcript_Nucleosome_Tallies"
    histone_tally_files = get_file_list(histone_tally_dir)

    normalize = False

    combined_list = get_combined_list(histone_tally_files, normalize)

    ylim_min = round(min(combined_list)-.4)
    ylim_max = round(max(combined_list)+.6)
    if normalize == True:
        ylim_min = round(min(combined_list)-.04, 2)
        ylim_max = round(max(combined_list)+.06, 2)
    bin_total = 400

    file_ID_list = BFL.unique_ID_list_from_string_list(histone_tally_files, "_", "_")
    print(file_ID_list)

    range_list = make_range_list(ylim_min, ylim_max, bin_total)
    tally_dict = make_tally_dict(file_ID_list, range_list)

    i = 0
    while i < len(histone_tally_files):
        file = histone_tally_files[i]
        file_ID = file_ID_list[i]
        tally_banana(file, file_ID, range_list, tally_dict, normalize)
        i += 1

    args = []
    for file_ID in file_ID_list:
        array_list = []
        for tally_key in tally_dict:
            tally = tally_dict[tally_key][file_ID]
            array_list.append(tally)
        args.append(array_list)

    H, p_omnibus, Z_pairs, p_corrected, reject = kw_dunn(args, to_compare=None, alpha=0.05, method='b')
    name_combinations = tuple(combinations(file_ID_list, 2))

    outfile = open("KW_Transcript_Mutations.txt", "w")
    header = prepare_out_new_line(["", "Dunn Test", "Reject Null", "KW p-value:", p_omnibus])
    outfile.write(header + '\n')
    i = 0
    while i < len(name_combinations):
        comparison_name = ' vs '.join(name_combinations[i])
        new_line = prepare_out_new_line([comparison_name, p_corrected[i], reject[i]])
        outfile.write(new_line + '\n')
        i+=1

    outfile = open("KW_Transcript_Binned_Mutations.txt", "w")
    header = prepare_out_new_line([""] + file_ID_list)
    outfile.write(header + '\n')
    for tally_key in tally_dict:
        new_line = [tally_key]
        for file_ID in file_ID_list:
            tally = tally_dict[tally_key][file_ID]
            new_line.append(tally)
        new_line = prepare_out_new_line(new_line)
        outfile.write(new_line + '\n')

import_BFL()
saddle_up()
