'''
Description: See Title
Author: Alexander Brown
Date: 6/4/2018
'''

import multiprocessing

import sys
import csv
import os
from os import listdir
from os.path import isfile, join

import astropy.stats as astro

def get_file_list(directory):
    file_list = [directory + '/' + f for f in listdir(directory) if isfile(join(directory, f))]
    return file_list

def get_expected_mut_dict(total_mut_path):
    expected_mut_dict = {}
    reader = csv.reader(open(total_mut_path, "r"), delimiter='\t')
    header = next(reader)
    for line in reader:
        bp = int(line[0]) -75
        e_count = float(line[2])
        expected_mut_dict[bp] = e_count
    return expected_mut_dict

def generate_bin_dict(num_of_intervals):
    bin_dict = {}
    x = 1
    while x <= num_of_intervals:
        bin_dict[x] = 0
        x += 1
    return bin_dict

def work_horse(path_list, expected_mut_dict, outfile):
    periodicity_list = []
    for path in path_list:
        x_list = []
        y_list = []
        reader = csv.reader(open(path, "r"), delimiter='\t')
        header = next(reader)
        for line in reader:
            bp = int(line[0]) -75
            o_count = int(line[1])
            e_count = expected_mut_dict[bp]
            try:
                enrichment = o_count/e_count
            except ZeroDivisionError:
                enrichment = 0
            x_list.append(bp)
            y_list.append(enrichment)

        del x_list[0]
        del x_list[-1]
        del y_list[0]
        del y_list[-1]

###############LOMB SCARGLE###############################
        f, Pxx_den = astro.LombScargle(x_list, y_list, fit_mean=False, center_data=True, nterms=1).autopower()
        new_Pxx_den = [y for x,y in zip(f, Pxx_den) if x < 0.5] #0.01 < x < 0.5
        Pxx_den_max = max(new_Pxx_den)
        len_f_array = len(f)
        i = 0
        while i < len_f_array:
            if Pxx_den[i] != Pxx_den_max:
                i += 1
                continue
            f_of_Pxx_den_max = f[i]
            i += 1
        periodicity_from_Pxx_den_max = 1/f_of_Pxx_den_max
#########################################################
        periodicity_list.append(periodicity_from_Pxx_den_max)

    periodicity_list.sort()
    periodicity_dict = {}
    for periodicity in periodicity_list:
        try:
            periodicity_dict[periodicity] += 1
        except KeyError:
            periodicity_dict[periodicity] = 1

    median_index = round(len(periodicity_list)/2)
    median = periodicity_list[median_index]
    header = '\t'.join([str(s) for s in ["Median:", median]])
    outfile.write(header + '\n')

    for key in periodicity_dict:
        count = periodicity_dict[key]
        new_line = '\t'.join([str(s) for s in [key, count]])
        outfile.write(new_line + '\n')

def saddle_up():
    directory = r"I:\Melanoma_2018_Nuc_Project\1000_Subsets"
    path_list = get_file_list(directory)

    total_mut_path = r"I:\Melanoma_2018_Nuc_Project\Cutaneous_Strong_Dyad_Mutation_And_Dinuc_Tallies.txt"
    expected_mut_dict = get_expected_mut_dict(total_mut_path)

    outfile = open("1000_Subsets_Periodicities.txt", "w")

    work_horse(path_list, expected_mut_dict, outfile)

saddle_up()
