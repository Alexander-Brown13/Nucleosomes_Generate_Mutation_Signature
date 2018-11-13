'''
Description: See Title
Author: Alexander Brown
Date: 1/15/2018
'''

import csv
import os
from os import listdir
from os.path import isfile, join

import numpy as np
import statistics
import astropy.stats as astro
import matplotlib
import matplotlib.pyplot as plt
matplotlib.rcParams['pdf.fonttype'] = 42

def get_file_list(directory):
    file_list = [directory + '/' + f for f in listdir(directory) if isfile(join(directory, f))]
    return file_list

def get_enrichment_list(file_list):
    enrichment_list = []
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
        
#        enrichment_normalization = sum(expected_count_list)/sum(observed_count_list)
 #       y_list = [y*enrichment_normalization for y in y_list]
        enrichment_list.extend(y_list)
    return enrichment_list

def squared_num(number):
    squared_num = number*number
    return squared_num

def find_nth_occurrence_rindex(search_string, find_string, nth):
    start = 0
    end = len(search_string)
    while nth > 0:
        index = search_string.rindex(find_string, start, end)
        end = index - 1
        nth -= 1
    return index

def round_up_to_largest_base_10(number):
    str_number = str(number)
    len_number = len(str_number) + 1
    first_digit = int(str_number[0])
    i = 1
    str_largest_base_10 = "1"
    while i < len_number:
        str_largest_base_10 = str_largest_base_10 + "0"
        i += 1
    largest_base_10 = int(str_largest_base_10)
    return largest_base_10

def round_down_to_largest_base_10(number):
    str_number = str(number)
    len_number = len(str_number)
    first_digit = int(str_number[0])
    i = 1
    str_largest_base_10 = "1"
    while i < len_number:
        str_largest_base_10 = str_largest_base_10 + "0"
        i += 1
    largest_base_10 = int(str_largest_base_10)
    return largest_base_10

def work_horse(file, num_of_files, figure, subplot_index, ylim_dict):
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

#    enrichment_normalization = sum(expected_count_list)/sum(observed_count_list)
 #   y_list = [y*enrichment_normalization for y in y_list]
    y_list_mean = round(statistics.mean(y_list))

    best_fit_coefficients = np.polyfit(x_list, y_list, 2)
    best_fit_polynomial = np.poly1d(best_fit_coefficients)
    best_fit_y_list = []
    for x_value in x_list:
        best_fit_y_list.append(best_fit_polynomial(x_value))
    first_deriv = np.polyder(best_fit_polynomial)
    second_deriv = np.polyder(first_deriv)

    SSR = sum([squared_num(y - best_fit_polynomial(x)) for x,y in zip(x_list, y_list)])
    SST = sum([squared_num(y - y_list_mean) for y in y_list])
    R_2 = 1 - (SSR/SST)

###############LOMB SCARGLE###############################
    f, Pxx_den = astro.LombScargle(x_list, y_list, fit_mean=False, center_data=True, nterms=1).autopower()
    new_Pxx_den = [y for x,y in zip(f, Pxx_den) if x < 0.5]
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
    
    columns = 3
    rows = round((num_of_files/columns)+.5)
    plt.subplot(rows, columns, subplot_index)
    plt.plot(x_list, y_list, 'royalblue', x_list, best_fit_y_list, 'c--')

    axv_line_list = [-72.1, -61.8, -51.5, -41.2, -30.9, -20.6, -10.3, 0, 10.3, 20.6, 30.9, 41.2, 51.5, 61.8, 72.1]
    for x_pos in axv_line_list:
        plt.axvline(x=x_pos, color='gray', linestyle='--', lw=1)
    every_other_axv_line_list = []
    i = 0
    while i < len(axv_line_list):
        every_other_axv_line_list.append(axv_line_list[i])
        i += 2
    plt.xticks(visible=False)
    if subplot_index > num_of_files - columns:
        plt.xlabel("bp Position")
        plt.xticks(every_other_axv_line_list, visible=True)
    if (subplot_index-1) % 3 == 0:
        plt.ylabel("Enrichment (Local/Global)")

    mut_count_mean = round(statistics.mean(observed_count_list))
    x_list_min = min(x_list)
    x_list_max = max(x_list)
    y_min = ylim_dict["Min"]*.6
    y_max = ylim_dict["Max"]*1.1
    plt.ylim(y_min*1.1, y_max*1)
    bbox_color = 'red' if mut_count_mean < 200 else 'green'
    bbox_specs = {'facecolor':bbox_color, 'alpha':0.5, 'pad':2}
    y_min_shift = 1.1
    y_text_coordinate = y_min * y_min_shift
    plt.text(x_list_min, y_text_coordinate, "Nuc Count: " + '\n' + str(nuc_count),
             horizontalalignment='left', color='black', fontsize=12,
             bbox=bbox_specs)

    plt.text(x_list_max, y_text_coordinate, "Second Der: " + str(second_deriv),
             horizontalalignment='right', color='black', weight='bold', fontsize=10)

    plt.text(0, y_text_coordinate, "Period: " + '\n' + str(round(periodicity_from_Pxx_den_max, 4)),
             horizontalalignment='center', color='black', weight='bold', fontsize=10)

    ninth_underscore_index = find_nth_occurrence_rindex(file, "_", 9)
    eighth_underscore_index = find_nth_occurrence_rindex(file, "_", 8)
    subplot_title = file[ninth_underscore_index+1:eighth_underscore_index]
    plt.title(subplot_title, size="medium")

def saddle_up():
    histone_tally_dir = r"I:\Melanoma_2018_Nuc_Project\Transcript_Nucleosome_Tallies"
    histone_tally_files = get_file_list(histone_tally_dir)
    num_of_files = len(histone_tally_files)

    figure = plt.figure(figsize=(15,10))

    enrichment_list = get_enrichment_list(histone_tally_files)
    ylim_min = min(enrichment_list)
    ylim_max = max(enrichment_list)
    ylim_dict = {"Min" : ylim_min, "Max" : ylim_max}

    subplot_index_dict = {"Low":1, "Middle":2, "High":3}
    for file in histone_tally_files:
        ninth_underscore_index = find_nth_occurrence_rindex(file, "_", 9)
        eighth_underscore_index = find_nth_occurrence_rindex(file, "_", 8)
        subplot_index = subplot_index_dict[file[ninth_underscore_index+1:eighth_underscore_index]]

        work_horse(file, num_of_files, figure, subplot_index, ylim_dict)
    plt.tight_layout()
    plt.savefig("Transcript_Nuc_Tallies_Graphs.pdf", bbox_inches="tight", transparent=True)
#    plt.savefig("Transcript_Nuc_Tallies_Graphs_Normalized.png", bbox_inches="tight", transparent=True)
  #  plt.show()


saddle_up()
