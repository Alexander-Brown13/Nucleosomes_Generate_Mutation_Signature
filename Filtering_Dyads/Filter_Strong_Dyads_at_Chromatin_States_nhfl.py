'''
Description: See Title
Author: Alexander Brown
Date: 1/30/2018
'''

import gzip
import sys
import csv
import os
from os import listdir
from os.path import isfile, join

def get_file_list(directory):
    file_list = [directory + '/' + f for f in listdir(directory) if isfile(join(directory, f))]
    return file_list

def find_nth_occurrence_rindex(search_string, find_string, nth):
    start = 0
    end = len(search_string)
    while nth > 0:
        index = search_string.rindex(find_string, start, end)
        end = index - 1
        nth -= 1
    return index

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
        start_pos = int(dyad_line[1])# - 74 #-73 for nucleosome, -1 to tally dipyrimidine correctly
        end_pos = int(dyad_line[2])# + 74 #+73 fpr nucleosome, +1 to tally dipyrimidine correctly
        crd_list = [start_pos, end_pos]
        chr_dict[chromosome].append(crd_list)

    for chr_key in chr_dict:
        chr_dict[chr_key].sort()
    return chr_dict

def generate_chromatin_dict(file):
    chr_dict = {'chr1_':None, 'chr2_':None, 'chr3_':None, 'chr4_':None, 'chr5_':None, 'chr6_':None, 'chr7_':None, 'chr8_':None, 'chr9_':None, 'chr10_':None, 
		'chr11_':None, 'chr12_':None, 'chr13_':None, 'chr14_':None, 'chr15_':None, 'chr16_':None, 'chr17_':None, 'chr18_':None, 'chr19_':None, 'chr20_':None, 
                'chr21_':None, 'chr22_':None, 'chrX_':None, 'chrY_':None, 'chrM_':None}
    for chr_key in chr_dict:
        chr_dict[chr_key] = []

    chr_to_chr_dict = {'chr1':'chr1_', 'chr2':'chr2_', 'chr3':'chr3_', 'chr4':'chr4_', 'chr5':'chr5_', 'chr6':'chr6_', 'chr7':'chr7_', 'chr8':'chr8_', 'chr9':'chr9_', 'chr10':'chr10_', 
		'chr11':'chr11_', 'chr12':'chr12_', 'chr13':'chr13_', 'chr14':'chr14_', 'chr15':'chr15_', 'chr16':'chr16_', 'chr17':'chr17_', 'chr18':'chr18_', 'chr19':'chr19_', 'chr20':'chr20_', 
                'chr21':'chr21_', 'chr22':'chr22_', 'chrX':'chrX_', 'chrY':'chrY_', 'chrM':'chrM_'}

    reader = csv.reader(open(file, "r"), delimiter = '\t')
    for line in reader:
        chromosome = chr_to_chr_dict[line[0]]
        start_pos = int(line[1]) -1 #base 0
        end_pos = int(line[2]) #base 1
        chromatin_state = line[3]
        if "/" in chromatin_state:
            index = find_nth_occurrence_rindex(chromatin_state, "/", 1)
            chromatin_state = chromatin_state[:index]
        crd_list = [start_pos, end_pos, chromatin_state]
        chr_dict[chromosome].append(crd_list)

    for chr_key in chr_dict:
        chr_dict[chr_key].sort()
    return chr_dict

def get_chromatin_outfile_dict(file, outfile_dir):
    outfile_dict = {}
    reader = csv.reader(open(file, "r"), delimiter='\t')
    for line in reader:
        chromatin_state = line[3]
        if "/" in chromatin_state:
            index = find_nth_occurrence_rindex(chromatin_state, "/", 1)
            chromatin_state = chromatin_state[:index]
        if chromatin_state in outfile_dict:
            continue
        outfile = open(outfile_dir + "/" + chromatin_state + "_Nucleosomes.txt", "w")
        outfile_dict[chromatin_state] = outfile
    return outfile_dict

def bisect_search(a_list, query):
    low_bound = 0
    high_bound = len(a_list)-1
    guess = round((high_bound + low_bound)/2)
    exist = 1

    left_flank_start = a_list[guess][0]
    right_flank_end = a_list[guess][1]-1

    while not left_flank_start <= query <= right_flank_end:
        previous_guess = guess

        if right_flank_end < query:
            low_bound = guess
            guess = round((high_bound + low_bound)/2)
        if left_flank_start > query:
            high_bound = guess
            guess = round((high_bound + low_bound)/2)
        if guess == previous_guess:
            exist = 0
            break

        left_flank_start = a_list[guess][0]
        right_flank_end = a_list[guess][1]-1
    return guess, exist

def work_horse(dyad_dict, chromatin_state_dict, outfile_dict):
    chr_to_chr_dict = {'1':'chr1_', '2':'chr2_', '3':'chr3_', '4':'chr4_', '5':'chr5_', '6':'chr6_', '7':'chr7_', '8':'chr8_', '9':'chr9_', '10':'chr10_', 
		'11':'chr11_', '12':'chr12_', '13':'chr13_', '14':'chr14_', '15':'chr15_', '16':'chr16_', '17':'chr17_', '18':'chr18_', '19':'chr19_', '20':'chr20_', 
                '21':'chr21_', '22':'chr22_', 'X':'chrX_', 'Y':'chrY_', 'MT':'chrM_'}

    dinuc_list = ['TT', 'TC', 'CT', 'CC', 'AA', 'AG', 'GA', 'GG']
    T_C_list = ["TT", "TC", "CT", "CC"]
    A_G_list = ["AA", "AG", "GA", "GG"]

    for chromosome in dyad_dict:
        for coordinate_pair in dyad_dict[chromosome]:
            start_pos = coordinate_pair[0]
            end_pos = coordinate_pair[1]
            dyad_pos = round((start_pos+end_pos)/2)

            chromatin_state_list = chromatin_state_dict[chromosome]
            if len(chromatin_state_list) == 0:
                continue
            guess, exist = bisect_search(chromatin_state_list, dyad_pos)
            if exist == 0:
                continue

            chromatin_state = chromatin_state_list[guess][2]
            outfile = outfile_dict[chromatin_state]

            new_line = '\t'.join([str(s) for s in [chromosome[:-1], start_pos, end_pos]])
            outfile.write(new_line + '\n')

def saddle_up():
    dyad_file = r"I:\Melanoma\Nucleosomes\nucs_allchr_hg19_t10_notBLIST_NoOverlaps.txt"
    dyad_dict = generate_dyad_dict(dyad_file)

    chromatin_state_file = r"I:\Melanoma_2018_Nuc_Project\Chromatin_States\wgEncodeBroadHmmNhlfHMM.bed"
    outfile_dir = r"I:\Melanoma_2018_Nuc_Project\Chromatin_States\Chromatin_State_Sorted_Nucleosomes"
    outfile_dict = get_chromatin_outfile_dict(chromatin_state_file, outfile_dir)

    chromatin_state_dict = generate_chromatin_dict(chromatin_state_file)

    work_horse(dyad_dict, chromatin_state_dict, outfile_dict)

saddle_up()
