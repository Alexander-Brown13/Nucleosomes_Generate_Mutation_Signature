'''
Description: See Title
Author: Alexander Brown
Date: 1/27/2018
'''

import sys
import csv
import os
from os import listdir
from os.path import isfile, join

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

def make_trinuc_list():
    bp_list = ["A", "C", "G", "T"]
    ref_mut_list = ["C", "T"]
    dinuc_list = ["CC", "CT", "TC", "TT"]
    trinuc_list = []

    for ref_bp in ref_mut_list:
        for bp in bp_list:
            for bp_2 in bp_list:
                tri_nuc = bp + ref_bp + bp_2
                trinuc_list.append(tri_nuc)
    return trinuc_list

def make_trinuc_dict():
    bp_list = ["A", "C", "G", "T"]
    ref_mut_list = ["C", "T"]
    dinuc_list = ["CC", "CT", "TC", "TT"]
    trinuc_dict = {}

    for ref_bp in ref_mut_list:
        for bp in bp_list:
            for bp_2 in bp_list:
                tri_nuc = bp + ref_bp + bp_2
                trinuc_dict[tri_nuc] = 0
    return trinuc_dict

def reverse_complement_seq(seq):
    rev_comp_seq = seq.replace('A', 'X').replace('T', 'A').replace('X', 'T').replace('C', 'X').replace('G', 'C').replace('X', 'G')
    rev_comp_seq = rev_comp_seq[::-1]
    return rev_comp_seq

def complement_seq(seq):
    comp_seq = seq.replace('A', 'X').replace('T', 'A').replace('X', 'T').replace('C', 'X').replace('G', 'C').replace('X', 'G')
    return comp_seq

def work_horse(chr_dict, trinuc_tally_dict, total_dict):
    chr_to_chr_dict = {'1':'chr1_', '2':'chr2_', '3':'chr3_', '4':'chr4_', '5':'chr5_', '6':'chr6_', '7':'chr7_', '8':'chr8_', '9':'chr9_', '10':'chr10_', 
		'11':'chr11_', '12':'chr12_', '13':'chr13_', '14':'chr14_', '15':'chr15_', '16':'chr16_', '17':'chr17_', '18':'chr18_', '19':'chr19_', '20':'chr20_', 
                '21':'chr21_', '22':'chr22_', 'X':'chrX_', 'Y':'chrY_', 'MT':'chrM_'}

    mutation_list = ['TA', 'TC', 'TG', 'CA', 'CG', 'CT']

    total_length = 0
    for chromosome in chr_dict:
        if chromosome == "chrM_":
            continue
        chr_seq = chr_dict[chromosome]
        i = 0
        j = i + 3
        chr_len = len(chr_seq)
        while j < chr_len + 1:
            trinuc = chr_seq[i:j]
            try:
                trinuc_tally_dict[trinuc] += 1
            except KeyError:
                trinuc = reverse_complement_seq(trinuc) #make sure this is nested
                try:
                    trinuc_tally_dict[trinuc] += 1
                except KeyError:
                    pass
            i += 1
            j += 1
            total_dict["Total"] += 1
            total_length += 1
        print(chromosome)
        print(chr_len)
        print(trinuc_tally_dict)
        print(total_dict["Total"])
        print(total_length)

def different_horse(chr_dict):
    trinuc_tally_dict = {}
    total_dict = {"Total":0}
    trinuc_tally_dict = make_trinuc_dict()

    outfile = open("All_Trinuc_Tallies_Genome_Wide.txt", "w")
    work_horse(chr_dict, trinuc_tally_dict, total_dict)

    for trinuc in trinuc_tally_dict:
        count = trinuc_tally_dict[trinuc]
        new_line = '\t'.join([str(s) for s in [trinuc, count]])
        outfile.write(new_line + '\n')
    print(total_dict)

def saddle_up():
    trinuc_list = make_trinuc_list()
    chr_dict = generate_chr_dict()

    different_horse(chr_dict)

saddle_up()
