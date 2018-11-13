'''
Description: See Title
Author: Alexander Brown
Date: 1/18/2018
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

def get_genome_tally_dict():
    genome_tally_dict = {}
    genome_trinuc_file = r"I:\Melanoma_2018_Nuc_Project\All_Trinuc_Tallies_Genome_Wide.txt"
    reader = csv.reader(open(genome_trinuc_file, "r"), delimiter = '\t')
    for line in reader:
        trinuc = line[0]
        count = int(line[1])
        genome_tally_dict[trinuc] = count
    return genome_tally_dict

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

def work_horse(mut_file, chr_dict, trinuc_tally_dict):
    chr_to_chr_dict = {'1':'chr1_', '2':'chr2_', '3':'chr3_', '4':'chr4_', '5':'chr5_', '6':'chr6_', '7':'chr7_', '8':'chr8_', '9':'chr9_', '10':'chr10_', 
		'11':'chr11_', '12':'chr12_', '13':'chr13_', '14':'chr14_', '15':'chr15_', '16':'chr16_', '17':'chr17_', '18':'chr18_', '19':'chr19_', '20':'chr20_', 
                '21':'chr21_', '22':'chr22_', 'X':'chrX_', 'Y':'chrY_', 'MT':'chrM_'}

    mutation_list = ['TA', 'TC', 'TG', 'CA', 'CG', 'CT']

    mut_reader = csv.reader(open(mut_file, "r"), delimiter = '\t')
    line_count = 0
    header = next(mut_reader)
    for mut_line in mut_reader:
        line_count += 1

        chromosome = chr_to_chr_dict[mut_line[4]]
        if chromosome == 'chrM_':
            continue
        event_name = mut_line[9]
        reference_bp = mut_line[10]
        variant_bp = mut_line[12]
        nucleotide_list = ["A", "T", "C", "G"]
        if reference_bp not in nucleotide_list or variant_bp not in nucleotide_list:
            continue
        if reference_bp == variant_bp:
            continue
        mutation = reference_bp + variant_bp
        
        start_pos = int(mut_line[5]) - 1
        end_pos = start_pos + 1
        tri_nuc_start = start_pos - 1
        tri_nuc_end = end_pos + 1
        tri_nuc = chr_dict[chromosome][tri_nuc_start:tri_nuc_end]

        if mutation not in trinuc_tally_dict:
            mutation = complement_seq(mutation)
            tri_nuc = reverse_complement_seq(tri_nuc)
        try:
            trinuc_tally_dict[mutation][tri_nuc] += 1
        except KeyError:
            continue
    print(line_count)

def different_horse(mut_files, chr_dict):
    genome_tally_dict = get_genome_tally_dict()
    trinuc_tally_dict = {}
    mutation_list = ['TA', 'TC', 'TG', 'CA', 'CG', 'CT']
    for mutation in mutation_list:
        trinuc_tally_dict[mutation] = make_trinuc_dict()

    outfile = open("Acral_Cutaneous_Uveal_All_Trinuc_Tallies.txt", "w")
    for mut_file in mut_files:
        work_horse(mut_file, chr_dict, trinuc_tally_dict)

    for mutation in trinuc_tally_dict:
        outfile.write(mutation + '\n')
        for key in trinuc_tally_dict[mutation]:
            if mutation[0] != key[1]:
                continue
            mutation_count = trinuc_tally_dict[mutation][key]
            genome_trinuc_count = genome_tally_dict[key]
            mutation_fraction = mutation_count/genome_trinuc_count
            new_line = '\t'.join([str(s) for s in [key, mutation_count, mutation_fraction]])
            outfile.write(new_line + '\n')

def saddle_up():
    trinuc_list = make_trinuc_list()
    chr_dict = generate_chr_dict()

    aT_mut_file = r"I:\Melanoma_2018_Nuc_Project\Acral_Cutaneous_Mutations\Acral_Dinuc_T_Mutations_H.txt"
    aC_mut_file = r"I:\Melanoma_2018_Nuc_Project\Acral_Cutaneous_Mutations\Acral_Dinuc_C_Mutations_H.txt"
#    acral_mut_files = [aT_mut_file, aC_mut_file]
    a_mut_file = r"I:\Melanoma_2018_Nuc_Project\Acral_Cutaneous_Mutations\Acral_simple_somatic_mutation.open.MELA-AU_sort_dedup3.txt"
    acral_mut_files = [a_mut_file]

    cT_mut_file = r"I:\Melanoma_2018_Nuc_Project\Acral_Cutaneous_Mutations\Cutaneous_Dinuc_T_Mutations_H.txt"
    cC_mut_file = r"I:\Melanoma_2018_Nuc_Project\Acral_Cutaneous_Mutations\Cutaneous_Dinuc_C_Mutations_H.txt"
#    cutaneous_mut_files = [cT_mut_file, cC_mut_file]
    c_mut_file = r"I:\Melanoma_2018_Nuc_Project\Acral_Cutaneous_Mutations\Cutaneous_simple_somatic_mutation.open.MELA-AU_sort_dedup3.txt"
    cutaneous_mut_files = [c_mut_file, a_mut_file]

    original_mut_file = r"I:\Melanoma\Melanoma_Mutations\simple_somatic_mutation.open.MELA-AU_sort_dedup3.txt"

#    different_horse(acral_mut_files, chr_dict)
    different_horse([original_mut_file], chr_dict)

saddle_up()
