'''
Description: See Title
Author: Alexander Brown
Date: 2/28/2018
'''

import csv
import scipy.stats as ss
import numpy as np

from itertools import combinations

def import_BFL():
    import sys
    my_dir = "I:\Function_Libraries"
    sys.path.insert(1, my_dir)
    import Brown_Function_Library as BFL
    global BFL
    print("BFL Imported")

def work_horse():
    histone_tally_dir = r"I:\Melanoma_2018_Nuc_Project\Transcript_Nucleosome_Tallies"
    histone_tally_file_paths = BFL.get_file_path_list(histone_tally_dir)

    ylim_min = 1
    ylim_max = 147
    bin_total = 16
    range_list = BFL.make_range_list(ylim_min, ylim_max, bin_total)

    file_ID_list = BFL.unique_ID_list_from_string_list(histone_tally_file_paths, "_", "_")

    list_of_o_count_binned_lists = []
    list_of_e_count_binned_lists = []
    list_of_enrichment_binned_lists = []
    for file_path in  histone_tally_file_paths:
        reader = csv.reader(open(file_path, "r"), delimiter='\t')
        header = next(reader)
        enrichment_binned_list = [0 for pair in range_list]
        o_count_list = [0 for pair in range_list]
        e_count_list = [0 for pair in range_list]
        for line in reader:
            bp = int(line[0])-1
            guess, exist = BFL.bisect_search_V_between_including_values(bp, range_list)
            if exist == 0:
                continue
            enrichment = float(line[3])
            o_count = int(line[1])
            e_count = float(line[2])
            enrichment_binned_list[guess] += enrichment
            o_count_list[guess] += o_count
            e_count_list[guess] += e_count
#        normalization_factor = sum(e_count_list)/sum(o_count_list)
 #       enrichment_binned_list = [y*normalization_factor for y in enrichment_binned_list]
        list_of_o_count_binned_lists.append(o_count_list)
        list_of_e_count_binned_lists.append(e_count_list)
        list_of_enrichment_binned_lists.append(enrichment_binned_list)

    index_combinations = tuple(combinations(range(len(list_of_enrichment_binned_lists)), 2))

    outfile = open("Chi2_Transcript_Mutations.txt", "w")
    header = BFL.prepare_outfile_new_line(["", "p observed", "p expected"])
    outfile.write(header + '\n')
    for index_pair in index_combinations:
        file_ID_pair = file_ID_list[index_pair[0]] + " vs " + file_ID_list[index_pair[1]]
        o_data_set_1 = list_of_o_count_binned_lists[index_pair[0]]
        o_data_set_2 = list_of_o_count_binned_lists[index_pair[1]]
        obs = np.array([o_data_set_1, o_data_set_2])
        chi2, o_p, dof, expected = ss.chi2_contingency(obs)

        e_data_set_1 = list_of_e_count_binned_lists[index_pair[0]]
        e_data_set_2 = list_of_e_count_binned_lists[index_pair[1]]
        exp = np.array([e_data_set_1, e_data_set_2])
        chi2, e_p, dof, expected = ss.chi2_contingency(exp)

        new_line = BFL.prepare_outfile_new_line([file_ID_pair, o_p, e_p])
        outfile.write(new_line + '\n')

import_BFL()
work_horse()
