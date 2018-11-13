'''
Description: See Title
Author: Alexander Brown
Date: 1/10/2018
'''

import sys
import csv

def find_nth_occurrence_rindex(search_string, find_string, nth):
    start = 0
    end = len(search_string)
    while nth > 0:
        index = search_string.rindex(find_string, start, end)
        end = index - 1
        nth -= 1
    return index

def get_cutaneous_acral_ID_dict(cutaneous_acral_ID_file, donor_ID_to_ICGC_ID_file):
    cutaneous_acral_ID_dict = {}
    reader = csv.reader(open(cutaneous_acral_ID_file, "r"), delimiter = '\t')
    header = next(reader)
    for line in reader:
        donor_ID = line[0]
        subtype = line[6]
        donor_ID_list = donor_ID.split("_")
        donor_ID = '-'.join(donor_ID_list)
        cutaneous_acral_ID_dict[donor_ID] = subtype

    cutaneous_acral_ID_dict_2 = {}
    reader = csv.reader(open(donor_ID_to_ICGC_ID_file, "r"), delimiter = '\t')
    header = next(reader)
    for line in reader:
        donor_ID = line[6]
        ICGC_ID = line[5]
        cutaneous_acral_ID_dict_2[ICGC_ID] = cutaneous_acral_ID_dict[donor_ID]
    return cutaneous_acral_ID_dict_2

def work_horse(file, outfile_dict, cutaneous_acral_ID_dict):
    reader = csv.reader(open(file, "r"), delimiter = '\t')
    for line in reader:
#        chromosome = chr_to_chr_dict[mut_line[8]]
 #       reference_bp = mut_line[14]
  #      variant_bp = mut_line[16]
   #     mutation = reference_bp + variant_bp
        
    #    start_pos = int(mut_line[9]) - 1
        ICGC_ID = line[1]
        try:
            subtype = cutaneous_acral_ID_dict[ICGC_ID]
            if subtype == "Mucosal":
                continue
        except KeyError:
            continue
        outfile = outfile_dict[subtype]
        new_line = '\t'.join(line)
        outfile.write(new_line + '\n')

def saddle_up():
    cutaneous_acral_ID_file = r"I:\Melanoma_2018_Nuc_Project\Acral_Cutaneous_Mutations\Supp_Table_1_Hayward_et_al.txt"
    donor_ID_to_ICGC_ID_file = r"I:\Melanoma_2018_Nuc_Project\Acral_Cutaneous_Mutations\sample.MELA-AU.tsv"
    T_mut_file = r"I:\Melanoma\Melanoma_Mutations\Dinuc_T_Mutations_H.txt"
    C_mut_file = r"I:\Melanoma\Melanoma_Mutations\Dinuc_C_Mutations_H.txt"
#    mut_files = [T_mut_file, C_mut_file]
    original_mut_file = r"I:\Melanoma\Melanoma_Mutations\simple_somatic_mutation.open.MELA-AU_sort_dedup3.txt"
    mut_files = [original_mut_file]

    cutaneous_acral_ID_dict = get_cutaneous_acral_ID_dict(cutaneous_acral_ID_file, donor_ID_to_ICGC_ID_file)
    for f in mut_files:
        directory_index = find_nth_occurrence_rindex(f, "\\", 1)
        cutaneous_outfile = open("Cutaneous_" + f[directory_index+1:], "w")
        acral_outfile = open("Acral_" + f[directory_index+1:], "w")
        outfile_dict = {"Cutaneous" : cutaneous_outfile, "Acral" : acral_outfile}
        work_horse(f, outfile_dict, cutaneous_acral_ID_dict)

saddle_up()
