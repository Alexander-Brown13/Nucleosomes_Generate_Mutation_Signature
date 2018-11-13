'''
Description: See Title
Author: Alexander Brown
Date: 2/9/2018
'''

import sys
import csv
import os
from os import listdir
from os.path import isfile, join

def get_file_list(directory):
    file_list = [directory + '/' + f for f in listdir(directory) if isfile(join(directory, f))]
    return file_list

def generate_bin_dict(num_of_intervals):
    bin_dict = {}
    x = 1
    while x <= num_of_intervals:
        bin_dict[x] = 0
        x += 1

    return bin_dict

def find_nth_occurrence_rindex(search_string, find_string, nth):
    start = 0
    end = len(search_string)
    while nth > 0:
        index = search_string.rindex(find_string, start, end)
        end = index - 1
        nth -= 1
    return index

def generate_chr_dict():
    directory = r"I:\Melanoma\Cat_Cap_Chrs"
    chr_files = [directory + '/' + f for f in listdir(directory)]
    chr_dict = {'chr1_':None, 'chr2_':None, 'chr3_':None, 'chr4_':None, 'chr5_':None, 'chr6_':None, 'chr7_':None, 'chr8_':None, 'chr9_':None, 'chr10_':None, 
		'chr11_':None, 'chr12_':None, 'chr13_':None, 'chr14_':None, 'chr15_':None, 'chr16_':None, 'chr17_':None, 'chr18_':None, 'chr19_':None, 'chr20_':None, 
                'chr21_':None, 'chr22_':None, 'chrX_':None, 'chrY_':None, 'chrM_':None}

    maxInt = sys.maxsize
    decrement = True
    while decrement:
        decrement = False
        try:
            csv.field_size_limit(maxInt)
        except OverflowError:
            maxInt = int(maxInt/10)
            decrement = True

    for chr_file in chr_files:
        directory_index = find_nth_occurrence_rindex(chr_file, "/", 1)
        underscore_index = find_nth_occurrence_rindex(chr_file, "_", 2)
        chr_file_chr = chr_file[directory_index+1:underscore_index+1]
        for chr_key in chr_dict:
            if chr_key == chr_file_chr:
                chr_fh = csv.reader(open(chr_file, "r"), delimiter = '\t')
                first_line = next(chr_fh)
                second_line = next(chr_fh)
                chr_dict[chr_key] = ''.join(second_line)
    return chr_dict

def generate_dyad_dict(dyad_file):
    chr_dict = {'chr1_':None, 'chr2_':None, 'chr3_':None, 'chr4_':None, 'chr5_':None, 'chr6_':None, 'chr7_':None, 'chr8_':None, 'chr9_':None, 'chr10_':None, 
		'chr11_':None, 'chr12_':None, 'chr13_':None, 'chr14_':None, 'chr15_':None, 'chr16_':None, 'chr17_':None, 'chr18_':None, 'chr19_':None, 'chr20_':None, 
                'chr21_':None, 'chr22_':None, 'chrX_':None, 'chrY_':None, 'chrM_':None}
    for chr_key in chr_dict:
        chr_dict[chr_key] = []

    chr_to_chr_dict = {'chr1':'chr1_', 'chr2':'chr2_', 'chr3':'chr3_', 'chr4':'chr4_', 'chr5':'chr5_', 'chr6':'chr6_', 'chr7':'chr7_', 'chr8':'chr8_', 'chr9':'chr9_', 'chr10':'chr10_', 
		'chr11':'chr11_', 'chr12':'chr12_', 'chr13':'chr13_', 'chr14':'chr14_', 'chr15':'chr15_', 'chr16':'chr16_', 'chr17':'chr17_', 'chr18':'chr18_', 'chr19':'chr19_', 'chr20':'chr20_', 
                'chr21':'chr21_', 'chr22':'chr22_', 'chrX':'chrX_', 'chrY':'chrY_', 'chrMT':'chrM_'}

    dyad_reader = csv.reader(open(dyad_file, "r"), delimiter = '\t')
    for dyad_line in dyad_reader:
        chromosome = chr_to_chr_dict[dyad_line[0]]
        start_pos = int(dyad_line[1]) - 500 #-73 for nucleosome, -1 to tally dipyrimidine correctly
        end_pos = int(dyad_line[2]) + 500 #+73 fpr nucleosome, +1 to tally dipyrimidine correctly
        crd_list = [start_pos, end_pos]
        chr_dict[chromosome].append(crd_list)

    for chr_key in chr_dict:
        chr_dict[chr_key].sort()
    return chr_dict

def bisect_search(a_list, query):
    low_bound = 0
    high_bound = len(a_list)-1
    guess = round(((high_bound - low_bound)/2) + low_bound)
    exist = 1

    left_flank_start = a_list[guess][0]
    right_flank_end = a_list[guess][1]-1

    while not left_flank_start <= query <= right_flank_end:
        previous_guess = guess

        if right_flank_end < query:
            low_bound = guess
            guess = round(((high_bound - low_bound)/2) + low_bound)
        if left_flank_start > query:
            high_bound = guess
            guess = round(((high_bound - low_bound)/2) + low_bound)
        if guess == previous_guess:
            exist = 0
            break

        left_flank_start = a_list[guess][0]
        right_flank_end = a_list[guess][1]-1
    return guess, exist

def make_value_indices_list(guess, a_list, query):
    left_flank_start = a_list[guess][0]
    right_flank_end = a_list[guess][1] -1#-1 because of end

    try:
        while left_flank_start <= query <= right_flank_end and guess > -1:
            guess -= 1
            left_flank_start = a_list[guess][0]
            right_flank_end = a_list[guess][1] -1#-1 because of end

    except IndexError:
        pass
    guess += 1
    low_bound = guess

    left_flank_start = a_list[guess][0]
    right_flank_end = a_list[guess][1] -1#-1 because of end
    try:
        while left_flank_start <= query <= right_flank_end:
            guess += 1
            left_flank_start = a_list[guess][0]
            right_flank_end = a_list[guess][1] -1#-1 because of end

    except IndexError:
        pass
    guess -= 1
    high_bound = guess
    index_list = []
    while low_bound <= high_bound:
        index_list.append(low_bound)
        low_bound += 1
    return index_list

def work_horse(mut_file, dyad_dict, mut_tally_dict, dinuc_tally_dict):
    chr_to_chr_dict = {'1':'chr1_', '2':'chr2_', '3':'chr3_', '4':'chr4_', '5':'chr5_', '6':'chr6_', '7':'chr7_', '8':'chr8_', '9':'chr9_', '10':'chr10_', 
		'11':'chr11_', '12':'chr12_', '13':'chr13_', '14':'chr14_', '15':'chr15_', '16':'chr16_', '17':'chr17_', '18':'chr18_', '19':'chr19_', '20':'chr20_', 
                '21':'chr21_', '22':'chr22_', 'X':'chrX_', 'Y':'chrY_', 'MT':'chrM_'}

    dinuc_list = ['TT', 'TC', 'CT', 'CC', 'AA', 'AG', 'GA', 'GG']
    T_C_list = ["TT", "TC", "CT", "CC"]
    A_G_list = ["AA", "AG", "GA", "GG"]

    mut_reader = csv.reader(open(mut_file, "r"), delimiter = '\t')
    for line in mut_reader:
        chromosome = line[0] #so that it matches the key for the dyad_dict
        dyad_list = dyad_dict[chromosome]
        if len(dyad_list) == 0:
            continue
        directionality = line[3]
        pos_1 = int(line[1]) - 1
        pos_2 = int(line[2]) - 1
        pos_list = [pos_1, pos_2]

        for pos in pos_list:
            guess, exist = bisect_search(dyad_list, pos)
            if exist == 0:
                continue

            index_list = make_value_indices_list(guess, dyad_list, pos)
            for index in index_list:                
                dyad_ranges = dyad_list[index]
                start_pos = dyad_ranges[0]
                end_pos = dyad_ranges[1]

                bp_index = pos - start_pos + 1
                mut_tally_dict[bp_index] += 1

def stable_away(chr_dict, dyad_dict, mut_tally_dict, dinuc_tally_dict, number_of_intervals, outfile):
    outfile.write("Bin" + '\t' + "Lesion_Tally" + '\n')
    key = 1
    while key <= number_of_intervals:
        new_line = '\t'.join([str(s) for s in [key, mut_tally_dict[key]]])
        outfile.write(new_line + '\n')
        key += 1

def different_horse(mut_files, dyad_dict):
    example_file = mut_files[0]
    directory_index_1 = find_nth_occurrence_rindex(example_file, "\\", 1)
    directory_index_2 = find_nth_occurrence_rindex(example_file, "/", 1)
    outfile_name = example_file[directory_index_1+1 : directory_index_2] + "_Strong_Dyad_24hr_XR_Lesions_Tallies.txt"
#    outfile = open(outfile_name, "w")

    number_of_intervals = 1001
    mut_tally_dict = generate_bin_dict(number_of_intervals)
    dinuc_tally_dict = generate_bin_dict(number_of_intervals)
    for mut_file in mut_files:
        work_horse(mut_file, dyad_dict, mut_tally_dict, dinuc_tally_dict)

    return mut_tally_dict

def get_80J_tally_dict(file):
    number_of_intervals = 149
    tally_dict = generate_bin_dict(number_of_intervals)
    reader = csv.reader(open(file, "r"), delimiter='\t')
    header = next(reader)
    for line in reader:
        bp_pos = int(line[0])-1
        count_80J = int(line[2])
        tally_dict[bp_pos] += count_80J
    return tally_dict

def retrieve_nuc_count(dyad_dict):
    count = 0
    for chromosome in dyad_dict:
        for dyad in dyad_dict[chromosome]:
            count += 1
    return count

def saddle_up():
    dyad_file = r"I:\Melanoma\Nucleosomes\nucs_allchr_hg19_t10_notBLIST_NoOverlaps.txt"
    dyad_dict = generate_dyad_dict(dyad_file)

    files_1hr_XR = get_file_list(r"I:\XR-Seq\Mapped_SAMs\BEDs\1hr")
    files_naked_DMG_Seq = get_file_list(r"I:\HS-Damage-Seq\GM12878_CPD_20J_nakedDNA\BEDs")

    tally_dict_1hr_XR = different_horse(files_1hr_XR, dyad_dict)
    tally_dict_DMG_Seq = different_horse(files_naked_DMG_Seq, dyad_dict)

    outfile = open("1hr_XR_by_Damage_Seq_Strong_Dyad_Lesions_And_Dinuc_Tallies_1001_Window.txt", "w")
    nuc_count = retrieve_nuc_count(dyad_dict)
    header = '\t'.join(["Bin", "1hr_XR_Lesion_Tally", "Damage_Lesion_Tally", "1hr_XR_by_Damage_Tally", "Nuc_Count:", str(nuc_count)])
    outfile.write(header + '\n')
    len_of_dict = len(tally_dict_1hr_XR)
    i=1
    while i < len_of_dict+1:
        count_1hr_XR = tally_dict_1hr_XR[i]
        count_DMG = tally_dict_DMG_Seq[i]
        try:
            normalized_value = count_1hr_XR/count_DMG
        except ZeroDivisionError:
            normalized_value = 0
        new_line = '\t'.join([str(s) for s in [i, count_1hr_XR, count_DMG, normalized_value]])
        outfile.write(new_line + '\n')
        i += 1

saddle_up()
