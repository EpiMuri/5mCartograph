# 5mCartograph_smk_3_summarize_calculate_and_output.py -> summarizes data of subsets and calculates data of interest

import re
from datetime import datetime
import pandas as pd
import statistics
import csv
import sys
import os

# start time
startTime = datetime.now()
print('-> script 3 STARTS')

############################################################################################################################
############################################################################################################################

# inputs
input_number = len(sys.argv)-2
index_for_split = int(input_number/2+1)
fileCounts5mer_list = sys.argv[1:index_for_split]
fileCountsChr_list = sys.argv[index_for_split:len(sys.argv)-1]

# create location string for output 
file_out = sys.argv[-1].strip()

# file path log output (also includes results)
fileLOG_3_1 = file_out + '/5mCartograph_out3/5mCartograph_out3_log_and_results.txt'
fileLOG_3_2 = file_out + '/5mCartograph_out3/5mCartograph_out3_5mer_results.txt'
fileLOG_3_3 = file_out + '/5mCartograph_out3/5mCartograph_out3_chromosome_contig_results.txt'
print(fileLOG_3_1)
print(fileLOG_3_2)
print(fileLOG_3_3)

############################################################################################################################
############################################################################################################################

# create result folder
os.makedirs(file_out + '/5mCartograph_out3/', exist_ok=True)

############################################################################################################################
############################################################################################################################

# create variables

# collect all 5-mers and how often they are unmethylated and methylated
# [[5-mer, #_unmethylated, %_unmethylated, #_methylated, %_methylated], [...], ...]
all_5_mers = []

# list for all distributions per chromosome
chr_tigs_distribution = []

# list of all 2/3-mers by context
CG_context = 'CG'
CHG_context =  ['CAG', 'CCG', 'CTG']
CHH_context = ['CAA', 'CAC', 'CAT', 'CCA', 'CCC', 'CCT', 'CTA', 'CTC', 'CTT']

# lost of all 5-mers alphabetically
contexts = ['AACAA', 'AACAC', 'AACAG', 'AACAT', 'AACCA', 'AACCC', 'AACCG', 'AACCT', 'AACGA', 'AACGC', 'AACGG', 'AACGT', 'AACTA', 'AACTC', 'AACTG', 'AACTT', 
    'ACCAA', 'ACCAC', 'ACCAG', 'ACCAT', 'ACCCA', 'ACCCC', 'ACCCG', 'ACCCT', 'ACCGA', 'ACCGC', 'ACCGG', 'ACCGT', 'ACCTA', 'ACCTC', 'ACCTG', 'ACCTT', 
    'AGCAA', 'AGCAC', 'AGCAG', 'AGCAT', 'AGCCA', 'AGCCC', 'AGCCG', 'AGCCT', 'AGCGA', 'AGCGC', 'AGCGG', 'AGCGT', 'AGCTA', 'AGCTC', 'AGCTG', 'AGCTT', 
    'ATCAA', 'ATCAC', 'ATCAG', 'ATCAT', 'ATCCA', 'ATCCC', 'ATCCG', 'ATCCT', 'ATCGA', 'ATCGC', 'ATCGG', 'ATCGT', 'ATCTA', 'ATCTC', 'ATCTG', 'ATCTT',
    'CACAA', 'CACAC', 'CACAG', 'CACAT', 'CACCA', 'CACCC', 'CACCG', 'CACCT', 'CACGA', 'CACGC', 'CACGG', 'CACGT', 'CACTA', 'CACTC', 'CACTG', 'CACTT', 
    'CCCAA', 'CCCAC', 'CCCAG', 'CCCAT', 'CCCCA', 'CCCCC', 'CCCCG', 'CCCCT', 'CCCGA', 'CCCGC', 'CCCGG', 'CCCGT', 'CCCTA', 'CCCTC', 'CCCTG', 'CCCTT', 
    'CGCAA', 'CGCAC', 'CGCAG', 'CGCAT', 'CGCCA', 'CGCCC', 'CGCCG', 'CGCCT', 'CGCGA', 'CGCGC', 'CGCGG', 'CGCGT', 'CGCTA', 'CGCTC', 'CGCTG', 'CGCTT', 
    'CTCAA', 'CTCAC', 'CTCAG', 'CTCAT', 'CTCCA', 'CTCCC', 'CTCCG', 'CTCCT', 'CTCGA', 'CTCGC', 'CTCGG', 'CTCGT', 'CTCTA', 'CTCTC', 'CTCTG', 'CTCTT', 
    'GACAA', 'GACAC', 'GACAG', 'GACAT', 'GACCA', 'GACCC', 'GACCG', 'GACCT', 'GACGA', 'GACGC', 'GACGG', 'GACGT', 'GACTA', 'GACTC', 'GACTG', 'GACTT', 
    'GCCAA', 'GCCAC', 'GCCAG', 'GCCAT', 'GCCCA', 'GCCCC', 'GCCCG', 'GCCCT', 'GCCGA', 'GCCGC', 'GCCGG', 'GCCGT', 'GCCTA', 'GCCTC', 'GCCTG', 'GCCTT', 
    'GGCAA', 'GGCAC', 'GGCAG', 'GGCAT', 'GGCCA', 'GGCCC', 'GGCCG', 'GGCCT', 'GGCGA', 'GGCGC', 'GGCGG', 'GGCGT', 'GGCTA', 'GGCTC', 'GGCTG', 'GGCTT', 
    'GTCAA', 'GTCAC', 'GTCAG', 'GTCAT', 'GTCCA', 'GTCCC', 'GTCCG', 'GTCCT', 'GTCGA', 'GTCGC', 'GTCGG', 'GTCGT', 'GTCTA', 'GTCTC', 'GTCTG', 'GTCTT',
    'TACAA', 'TACAC', 'TACAG', 'TACAT', 'TACCA', 'TACCC', 'TACCG', 'TACCT', 'TACGA', 'TACGC', 'TACGG', 'TACGT', 'TACTA', 'TACTC', 'TACTG', 'TACTT', 
    'TCCAA', 'TCCAC', 'TCCAG', 'TCCAT', 'TCCCA', 'TCCCC', 'TCCCG', 'TCCCT', 'TCCGA', 'TCCGC', 'TCCGG', 'TCCGT', 'TCCTA', 'TCCTC', 'TCCTG', 'TCCTT', 
    'TGCAA', 'TGCAC', 'TGCAG', 'TGCAT', 'TGCCA', 'TGCCC', 'TGCCG', 'TGCCT', 'TGCGA', 'TGCGC', 'TGCGG', 'TGCGT', 'TGCTA', 'TGCTC', 'TGCTG', 'TGCTT', 
    'TTCAA', 'TTCAC', 'TTCAG', 'TTCAT', 'TTCCA', 'TTCCC', 'TTCCG', 'TTCCT', 'TTCGA', 'TTCGC', 'TTCGG', 'TTCGT', 'TTCTA', 'TTCTC', 'TTCTG', 'TTCTT']

# counter for all Cs
c_counter = 0

# 0 == 0 so that lines are collapsable
if 0 == 0:
    # (1) + (2) + (3)
    CG_unmet_CDS = 0
    CG_unmet_UTR = 0
    CG_unmet_intron = 0
    CG_unmet_non_genic = 0
    CG_unmet_upstream = 0
    CG_unmet_downstream = 0
    CG_rather_unmet_CDS = 0
    CG_rather_unmet_UTR = 0
    CG_rather_unmet_intron = 0
    CG_rather_unmet_non_genic = 0
    CG_rather_unmet_upstream = 0
    CG_rather_unmet_downstream = 0
    CG_50met_CDS = 0
    CG_50met_UTR = 0
    CG_50met_intron = 0
    CG_50met_non_genic = 0
    CG_50met_upstream = 0
    CG_50met_downstream = 0
    CG_rather_met_CDS = 0
    CG_rather_met_UTR = 0
    CG_rather_met_intron = 0
    CG_rather_met_non_genic = 0
    CG_rather_met_upstream = 0
    CG_rather_met_downstream = 0
    CG_met_CDS = 0
    CG_met_UTR = 0
    CG_met_intron = 0
    CG_met_non_genic = 0
    CG_met_upstream = 0
    CG_met_downstream = 0
    CHG_unmet_CDS = 0
    CHG_unmet_UTR = 0
    CHG_unmet_intron = 0
    CHG_unmet_non_genic = 0
    CHG_unmet_upstream = 0
    CHG_unmet_downstream = 0
    CHG_rather_unmet_CDS = 0
    CHG_rather_unmet_UTR = 0
    CHG_rather_unmet_intron = 0
    CHG_rather_unmet_non_genic = 0
    CHG_rather_unmet_upstream = 0
    CHG_rather_unmet_downstream = 0
    CHG_50met_CDS = 0
    CHG_50met_UTR = 0
    CHG_50met_intron = 0
    CHG_50met_non_genic = 0
    CHG_50met_upstream = 0
    CHG_50met_downstream = 0
    CHG_rather_met_CDS = 0
    CHG_rather_met_UTR = 0
    CHG_rather_met_intron = 0
    CHG_rather_met_non_genic = 0
    CHG_rather_met_upstream = 0
    CHG_rather_met_downstream = 0
    CHG_met_CDS = 0
    CHG_met_UTR = 0
    CHG_met_intron = 0
    CHG_met_non_genic = 0
    CHG_met_upstream = 0
    CHG_met_downstream = 0
    CHH_unmet_CDS = 0
    CHH_unmet_UTR = 0
    CHH_unmet_intron = 0
    CHH_unmet_non_genic = 0
    CHH_unmet_upstream = 0
    CHH_unmet_downstream = 0
    CHH_rather_unmet_CDS = 0
    CHH_rather_unmet_UTR = 0
    CHH_rather_unmet_intron = 0
    CHH_rather_unmet_non_genic = 0
    CHH_rather_unmet_upstream = 0
    CHH_rather_unmet_downstream = 0
    CHH_50met_CDS = 0
    CHH_50met_UTR = 0
    CHH_50met_intron = 0
    CHH_50met_non_genic = 0
    CHH_50met_upstream = 0
    CHH_50met_downstream = 0
    CHH_rather_met_CDS = 0
    CHH_rather_met_UTR = 0
    CHH_rather_met_intron = 0
    CHH_rather_met_non_genic = 0
    CHH_rather_met_upstream = 0
    CHH_rather_met_downstream = 0
    CHH_met_CDS = 0
    CHH_met_UTR = 0
    CHH_met_intron = 0
    CHH_met_non_genic = 0
    CHH_met_upstream = 0
    CHH_met_downstream = 0

    # (1) + (2)
    CG_CDS = 0
    CHG_CDS = 0
    CHH_CDS = 0
    CG_UTR = 0
    CHG_UTR = 0
    CHH_UTR = 0
    CG_exon = 0
    CHG_exon = 0
    CHH_exon = 0
    CG_intron = 0
    CHG_intron = 0
    CHH_intron = 0
    CG_genic = 0
    CHG_genic = 0
    CHH_genic = 0
    CG_non_genic = 0
    CHG_non_genic = 0
    CHH_non_genic = 0
    CG_upstream = 0
    CHG_upstream = 0
    CHH_upstream = 0
    CG_downstream = 0
    CHG_downstream = 0
    CHH_downstream = 0


    # (1) + (3)
    CG_unmet = 0
    CHG_unmet = 0
    CHH_unmet = 0
    CG_rather_unmet = 0
    CHG_rather_unmet = 0
    CHH_rather_unmet = 0
    CG_50met = 0
    CHG_50met = 0
    CHH_50met = 0
    CG_rather_met = 0
    CHG_rather_met = 0
    CHH_rather_met = 0
    CG_met = 0
    CHG_met = 0
    CHH_met = 0

    # (2) + (3)
    c_unmet_CDS = 0
    c_unmet_UTR = 0
    c_unmet_exon = 0
    c_unmet_intron = 0
    c_unmet_genic = 0
    c_unmet_non_genic = 0
    c_unmet_upstream = 0
    c_unmet_downstream = 0
    c_rather_unmet_CDS = 0
    c_rather_unmet_UTR = 0
    c_rather_unmet_exon = 0
    c_rather_unmet_intron = 0
    c_rather_unmet_genic = 0
    c_rather_unmet_non_genic = 0
    c_rather_unmet_upstream = 0
    c_rather_unmet_downstream = 0
    c_50met_CDS = 0
    c_50met_UTR = 0
    c_50met_exon = 0
    c_50met_intron = 0
    c_50met_genic = 0
    c_50met_non_genic = 0
    c_50met_upstream = 0
    c_50met_downstream = 0
    c_rather_met_CDS = 0
    c_rather_met_UTR = 0
    c_rather_met_exon = 0
    c_rather_met_intron = 0
    c_rather_met_genic = 0
    c_rather_met_non_genic = 0
    c_rather_met_upstream = 0
    c_rather_met_downstream = 0
    c_met_CDS = 0
    c_met_UTR = 0
    c_met_exon = 0
    c_met_intron = 0
    c_met_genic = 0
    c_met_non_genic = 0
    c_met_upstream = 0
    c_met_downstream = 0

    # (1)
    CG_fraction = 0
    CHG_fraction = 0
    CHH_fraction = 0

    # (2)
    c_CDS = 0
    c_UTR = 0
    c_exon = 0
    c_intron = 0
    c_genic = 0
    c_non_genic = 0
    c_upstream = 0
    c_downstream = 0

    # (3)
    c_unmet = 0
    c_rather_unmet = 0
    c_50met = 0
    c_rather_met = 0
    c_met = 0

    # variables for output of CG values in CG relation -> (1) + (2) + (3)
    CG_100_unmet_CDS = 0
    CG_100_unmet_UTR = 0
    CG_100_unmet_exon = 0
    CG_100_unmet_intron = 0
    CG_100_unmet_genic = 0
    CG_100_unmet_non_genic = 0
    CG_100_unmet_upstream = 0
    CG_100_unmet_downstream = 0
    CG_100_rather_unmet_CDS = 0
    CG_100_rather_unmet_UTR = 0
    CG_100_rather_unmet_exon = 0
    CG_100_rather_unmet_intron = 0
    CG_100_rather_unmet_genic = 0
    CG_100_rather_unmet_non_genic = 0
    CG_100_rather_unmet_upstream = 0
    CG_100_rather_unmet_downstream = 0
    CG_100_50met_CDS = 0
    CG_100_50met_UTR = 0
    CG_100_50met_exon = 0
    CG_100_50met_intron = 0
    CG_100_50met_genic = 0
    CG_100_50met_non_genic = 0
    CG_100_50met_upstream = 0
    CG_100_50met_downstream = 0
    CG_100_rather_met_CDS = 0
    CG_100_rather_met_UTR = 0
    CG_100_rather_met_exon = 0
    CG_100_rather_met_intron = 0
    CG_100_rather_met_genic = 0
    CG_100_rather_met_non_genic = 0
    CG_100_rather_met_upstream = 0
    CG_100_rather_met_downstream = 0
    CG_100_met_CDS = 0
    CG_100_met_UTR = 0
    CG_100_met_exon = 0
    CG_100_met_intron = 0
    CG_100_met_genic = 0
    CG_100_met_non_genic = 0
    CG_100_met_upstream = 0
    CG_100_met_downstream = 0
    CHG_100_unmet_CDS = 0
    CHG_100_unmet_UTR = 0
    CHG_100_unmet_exon = 0
    CHG_100_unmet_intron = 0
    CHG_100_unmet_genic = 0
    CHG_100_unmet_non_genic = 0
    CHG_100_unmet_upstream = 0
    CHG_100_unmet_downstream = 0
    CHG_100_rather_unmet_CDS = 0
    CHG_100_rather_unmet_UTR = 0
    CHG_100_rather_unmet_exon = 0
    CHG_100_rather_unmet_intron = 0
    CHG_100_rather_unmet_genic = 0
    CHG_100_rather_unmet_non_genic = 0
    CHG_100_rather_unmet_upstream = 0
    CHG_100_rather_unmet_downstream = 0
    CHG_100_50met_CDS = 0
    CHG_100_50met_UTR = 0
    CHG_100_50met_exon = 0
    CHG_100_50met_intron = 0
    CHG_100_50met_genic = 0
    CHG_100_50met_non_genic = 0
    CHG_100_50met_upstream = 0
    CHG_100_50met_downstream = 0
    CHG_100_rather_met_CDS = 0
    CHG_100_rather_met_UTR = 0
    CHG_100_rather_met_exon = 0
    CHG_100_rather_met_intron = 0
    CHG_100_rather_met_genic = 0
    CHG_100_rather_met_non_genic = 0
    CHG_100_rather_met_upstream = 0
    CHG_100_rather_met_downstream = 0
    CHG_100_met_CDS = 0
    CHG_100_met_UTR = 0
    CHG_100_met_exon = 0
    CHG_100_met_intron = 0
    CHG_100_met_genic = 0
    CHG_100_met_non_genic = 0
    CHG_100_met_upstream = 0
    CHG_100_met_downstream = 0
    CHH_100_unmet_CDS = 0
    CHH_100_unmet_UTR = 0
    CHH_100_unmet_exon = 0
    CHH_100_unmet_intron = 0
    CHH_100_unmet_genic = 0
    CHH_100_unmet_non_genic = 0
    CHH_100_unmet_upstream = 0
    CHH_100_unmet_downstream = 0
    CHH_100_rather_unmet_CDS = 0
    CHH_100_rather_unmet_UTR = 0
    CHH_100_rather_unmet_exon = 0
    CHH_100_rather_unmet_intron = 0
    CHH_100_rather_unmet_genic = 0
    CHH_100_rather_unmet_non_genic = 0
    CHH_100_rather_unmet_upstream = 0
    CHH_100_rather_unmet_downstream = 0
    CHH_100_50met_CDS = 0
    CHH_100_50met_UTR = 0
    CHH_100_50met_exon = 0
    CHH_100_50met_intron = 0
    CHH_100_50met_genic = 0
    CHH_100_50met_non_genic = 0
    CHH_100_50met_upstream = 0
    CHH_100_50met_downstream = 0
    CHH_100_rather_met_CDS = 0
    CHH_100_rather_met_UTR = 0
    CHH_100_rather_met_exon = 0
    CHH_100_rather_met_intron = 0
    CHH_100_rather_met_genic = 0
    CHH_100_rather_met_non_genic = 0
    CHH_100_rather_met_upstream = 0
    CHH_100_rather_met_downstream = 0
    CHH_100_met_CDS = 0
    CHH_100_met_UTR = 0
    CHH_100_met_exon = 0
    CHH_100_met_intron = 0
    CHH_100_met_genic = 0
    CHH_100_met_non_genic = 0
    CHH_100_met_upstream = 0
    CHH_100_met_downstream = 0

    # variables for output of CG values in CG relation -> (1) + (2)
    CG_100_CDS = 0
    CG_100_UTR = 0
    CG_100_exon = 0
    CG_100_intron = 0
    CG_100_genic = 0
    CG_100_non_genic = 0
    CG_100_upstream = 0
    CG_100_downstream = 0
    CHG_100_CDS = 0
    CHG_100_UTR = 0
    CHG_100_exon = 0
    CHG_100_intron = 0
    CHG_100_genic = 0
    CHG_100_non_genic = 0
    CHG_100_upstream = 0
    CHG_100_downstream = 0
    CHH_100_CDS = 0
    CHH_100_UTR = 0
    CHH_100_exon = 0
    CHH_100_intron = 0
    CHH_100_genic = 0
    CHH_100_non_genic = 0
    CHH_100_upstream = 0
    CHH_100_downstream = 0

    # variables for output of CG values in CG relation -> (1) + (3)
    CG_100_unmet = 0
    CG_100_rather_unmet = 0
    CG_100_50met = 0
    CG_100_rather_met = 0
    CG_100_met = 0
    CHG_100_unmet = 0
    CHG_100_rather_unmet = 0
    CHG_100_50met = 0
    CHG_100_rather_met = 0
    CHG_100_met = 0
    CHH_100_unmet = 0
    CHH_100_rather_unmet = 0
    CHH_100_50met = 0
    CHH_100_rather_met = 0
    CHH_100_met = 0

    print('-> variables were created')

############################################################################################################################
############################################################################################################################

# list of lists of lists for storing all 5mer-counter data
counts5mer_list = []
# read in all 5mer files
for file in fileCounts5mer_list:
    with open(file, 'r') as file5mer:
        # get the whole list of lists
        entry5mer = list(csv.reader(file5mer, delimiter='\t'))
        # convert string to int
        for idx, line in enumerate(entry5mer):
            for udx, count in enumerate(line):
                if udx != 0:
                    entry5mer[idx][udx] = int(count)
        # safe list in a list of lists of lists
        counts5mer_list.append(entry5mer)
# print(counts5mer_list)

# list of lists of lists for storing all chr-counter data
countsChr_list = []
# read in all chr files
for file in fileCountsChr_list:
    with open(file, 'r') as fileChr:
        # get the whole list of lists
        entryChr = list(csv.reader(fileChr, delimiter='\t'))
        # convert string to int
        for idx, line in enumerate(entryChr):
            for udx, count in enumerate(line):
                if udx != 0:
                    entryChr[idx][udx] = int(count)
        # safe list in a list of lists of lists
        countsChr_list.append(entryChr)
# print(countsChr_list)

############################################################################################################################
############################################################################################################################

# create empty list of lists with all methylation contexts
for context in contexts:
    current_5_mer_list = [0] * 61
    current_5_mer_list[0] = context
    all_5_mers.append(current_5_mer_list)    

# calculate sums necessary
for idx, kmer in enumerate(all_5_mers):
    for udx, value in enumerate(kmer):
        if udx != 0:
            for edx, countlist in enumerate(counts5mer_list):
                all_5_mers[idx][udx] += counts5mer_list[edx][idx][udx]
                c_counter += counts5mer_list[edx][idx][udx]

############################################################################################################################

# create list for all chromosomes and tigs
for core in countsChr_list:    
    for countsChr in core:
        chr_tigs_distribution.append(countsChr)
        # print(countsChr)
        # print('#################')

############################################################################################################################
############################################################################################################################

# context (1)       region (2)          methylation status
# |                 |                   |
# CG_fraction...    c_CDS...            c_unmet...

# (1) + (2)             CG_CDS...
# (1) + (3)             CG_unmet...
# (2) + (3)             c_unmet_CDS...
# (1) + (2) + (3)       CG_unmet_CDS...

############################################################################################################################

# calculate realative frequencie & calculate values of interest -> context and exon, intron ... and methylation status
for k_mer in all_5_mers:
    # calculate realtive freqneuncies
    k_mer[2] = k_mer[1]/c_counter*100
    k_mer[4] = k_mer[3]/c_counter*100
    k_mer[6] = k_mer[5]/c_counter*100
    k_mer[8] = k_mer[7]/c_counter*100
    k_mer[10] = k_mer[9]/c_counter*100
    k_mer[12] = k_mer[11]/c_counter*100
    k_mer[14] = k_mer[13]/c_counter*100
    k_mer[16] = k_mer[15]/c_counter*100
    k_mer[18] = k_mer[17]/c_counter*100
    k_mer[20] = k_mer[19]/c_counter*100
    k_mer[22] = k_mer[21]/c_counter*100
    k_mer[24] = k_mer[23]/c_counter*100
    k_mer[26] = k_mer[25]/c_counter*100
    k_mer[28] = k_mer[27]/c_counter*100
    k_mer[30] = k_mer[29]/c_counter*100
    k_mer[32] = k_mer[31]/c_counter*100
    k_mer[34] = k_mer[33]/c_counter*100
    k_mer[36] = k_mer[35]/c_counter*100
    k_mer[38] = k_mer[37]/c_counter*100
    k_mer[40] = k_mer[39]/c_counter*100
    k_mer[42] = k_mer[41]/c_counter*100
    k_mer[44] = k_mer[43]/c_counter*100
    k_mer[46] = k_mer[45]/c_counter*100
    k_mer[48] = k_mer[47]/c_counter*100
    k_mer[50] = k_mer[49]/c_counter*100
    k_mer[52] = k_mer[51]/c_counter*100
    k_mer[54] = k_mer[53]/c_counter*100
    k_mer[56] = k_mer[55]/c_counter*100
    k_mer[58] = k_mer[57]/c_counter*100
    k_mer[60] = k_mer[59]/c_counter*100

    # print(k_mer)
    
    # calculate values of interest -> 5-mers and exons, introns in CG, CHG and CHH context
    if k_mer[0][2:4] in CG_context:
        # exon, intron ... in CG -> (1) + (2) + (3)
        CG_unmet_CDS += k_mer[2]
        CG_unmet_UTR += k_mer[12]
        CG_unmet_intron += k_mer[22]
        CG_unmet_non_genic += k_mer[32]
        CG_unmet_upstream += k_mer[42]
        CG_unmet_downstream += k_mer[52]
        CG_rather_unmet_CDS += k_mer[4]
        CG_rather_unmet_UTR += k_mer[14]
        CG_rather_unmet_intron += k_mer[24]
        CG_rather_unmet_non_genic += k_mer[34]
        CG_rather_unmet_upstream += k_mer[44]
        CG_rather_unmet_downstream += k_mer[54]
        CG_50met_CDS += k_mer[6]
        CG_50met_UTR += k_mer[16]
        CG_50met_intron += k_mer[26]
        CG_50met_non_genic += k_mer[36]
        CG_50met_upstream += k_mer[46]
        CG_50met_downstream += k_mer[56]
        CG_rather_met_CDS += k_mer[8]
        CG_rather_met_UTR += k_mer[18]
        CG_rather_met_intron += k_mer[28]
        CG_rather_met_non_genic += k_mer[38]
        CG_rather_met_upstream += k_mer[48]
        CG_rather_met_downstream += k_mer[58]
        CG_met_CDS += k_mer[10]
        CG_met_UTR += k_mer[20]
        CG_met_intron += k_mer[30]
        CG_met_non_genic += k_mer[40]
        CG_met_upstream += k_mer[50]
        CG_met_downstream += k_mer[60]
    elif k_mer[0][2:] in CHG_context:
        # exon, intron ... in CHG -> (1) + (2) + (3)
        CHG_unmet_CDS += k_mer[2]
        CHG_unmet_UTR += k_mer[12]
        CHG_unmet_intron += k_mer[22]
        CHG_unmet_non_genic += k_mer[32]
        CHG_unmet_upstream += k_mer[42]
        CHG_unmet_downstream += k_mer[52]
        CHG_rather_unmet_CDS += k_mer[4]
        CHG_rather_unmet_UTR += k_mer[14]
        CHG_rather_unmet_intron += k_mer[24]
        CHG_rather_unmet_non_genic += k_mer[34]
        CHG_rather_unmet_upstream += k_mer[44]
        CHG_rather_unmet_downstream += k_mer[54]
        CHG_50met_CDS += k_mer[6]
        CHG_50met_UTR += k_mer[16]
        CHG_50met_intron += k_mer[26]
        CHG_50met_non_genic += k_mer[36]
        CHG_50met_upstream += k_mer[46]
        CHG_50met_downstream += k_mer[56]
        CHG_rather_met_CDS += k_mer[8]
        CHG_rather_met_UTR += k_mer[18]
        CHG_rather_met_intron += k_mer[28]
        CHG_rather_met_non_genic += k_mer[38]
        CHG_rather_met_upstream += k_mer[48]
        CHG_rather_met_downstream += k_mer[58]
        CHG_met_CDS += k_mer[10]
        CHG_met_UTR += k_mer[20]
        CHG_met_intron += k_mer[30]
        CHG_met_non_genic += k_mer[40]
        CHG_met_upstream += k_mer[50]
        CHG_met_downstream += k_mer[60]
    elif k_mer[0][2:] in CHH_context:
        # exon, intron ... in CHH -> (1) + (2) + (3)
        CHH_unmet_CDS += k_mer[2]
        CHH_unmet_UTR += k_mer[12]
        CHH_unmet_intron += k_mer[22]
        CHH_unmet_non_genic += k_mer[32]
        CHH_unmet_upstream += k_mer[42]
        CHH_unmet_downstream += k_mer[52]
        CHH_rather_unmet_CDS += k_mer[4]
        CHH_rather_unmet_UTR += k_mer[14]
        CHH_rather_unmet_intron += k_mer[24]
        CHH_rather_unmet_non_genic += k_mer[34]
        CHH_rather_unmet_upstream += k_mer[44]
        CHH_rather_unmet_downstream += k_mer[54]
        CHH_50met_CDS += k_mer[6]
        CHH_50met_UTR += k_mer[16]
        CHH_50met_intron += k_mer[26]
        CHH_50met_non_genic += k_mer[36]
        CHH_50met_upstream += k_mer[46]
        CHH_50met_downstream += k_mer[56]
        CHH_rather_met_CDS += k_mer[8]
        CHH_rather_met_UTR += k_mer[18]
        CHH_rather_met_intron += k_mer[28]
        CHH_rather_met_non_genic += k_mer[38]
        CHH_rather_met_upstream += k_mer[48]
        CHH_rather_met_downstream += k_mer[58]
        CHH_met_CDS += k_mer[10]
        CHH_met_UTR += k_mer[20]
        CHH_met_intron += k_mer[30]
        CHH_met_non_genic += k_mer[40]
        CHH_met_upstream += k_mer[50]
        CHH_met_downstream += k_mer[60]

print('-> data arrays were filled with relative data')

# calculate proportion of CDS, ... in CG, CHG and CHH context -> (1) + (2)
CG_CDS = CG_unmet_CDS + CG_rather_unmet_CDS + CG_50met_CDS + CG_rather_met_CDS + CG_met_CDS
CHG_CDS = CHG_unmet_CDS + CHG_rather_unmet_CDS + CHG_50met_CDS + CHG_rather_met_CDS + CHG_met_CDS
CHH_CDS = CHH_unmet_CDS + CHH_rather_unmet_CDS + CHH_50met_CDS + CHH_rather_met_CDS + CHH_met_CDS
CG_UTR = CG_unmet_UTR + CG_rather_unmet_UTR + CG_50met_UTR + CG_rather_met_UTR + CG_met_UTR
CHG_UTR = CHG_unmet_UTR + CHG_rather_unmet_UTR + CHG_50met_UTR + CHG_rather_met_UTR + CHG_met_UTR
CHH_UTR = CHH_unmet_UTR + CHH_rather_unmet_UTR + CHH_50met_UTR + CHH_rather_met_UTR + CHH_met_UTR
CG_exon = CG_CDS + CG_UTR
CHG_exon = CHG_CDS + CHG_UTR
CHH_exon = CHH_CDS + CHH_UTR
CG_intron = CG_unmet_intron + CG_rather_unmet_intron + CG_50met_intron + CG_rather_met_intron + CG_met_intron
CHG_intron = CHG_unmet_intron + CHG_rather_unmet_intron + CHG_50met_intron + CHG_rather_met_intron + CHG_met_intron
CHH_intron = CHH_unmet_intron + CHH_rather_unmet_intron + CHH_50met_intron + CHH_rather_met_intron + CHH_met_intron
CG_genic = CG_exon + CG_intron
CHG_genic = CHG_exon + CHG_intron
CHH_genic = CHH_exon + CHH_intron
CG_non_genic = CG_unmet_non_genic + CG_rather_unmet_non_genic + CG_50met_non_genic + CG_rather_met_non_genic + CG_met_non_genic
CHG_non_genic = CHG_unmet_non_genic + CHG_rather_unmet_non_genic + CHG_50met_non_genic + CHG_rather_met_non_genic + CHG_met_non_genic
CHH_non_genic = CHH_unmet_non_genic + CHH_rather_unmet_non_genic + CHH_50met_non_genic + CHH_rather_met_non_genic + CHH_met_non_genic
CG_upstream = CG_unmet_upstream + CG_rather_unmet_upstream + CG_50met_upstream + CG_rather_met_upstream + CG_met_upstream
CHG_upstream = CHG_unmet_upstream + CHG_rather_unmet_upstream + CHG_50met_upstream + CHG_rather_met_upstream + CHG_met_upstream
CHH_upstream = CHH_unmet_upstream + CHH_rather_unmet_upstream + CHH_50met_upstream + CHH_rather_met_upstream + CHH_met_upstream
CG_downstream = CG_unmet_downstream + CG_rather_unmet_downstream + CG_50met_downstream + CG_rather_met_downstream + CG_met_downstream
CHG_downstream = CHG_unmet_downstream + CHG_rather_unmet_downstream + CHG_50met_downstream + CHG_rather_met_downstream + CHG_met_downstream
CHH_downstream = CHH_unmet_downstream + CHH_rather_unmet_downstream + CHH_50met_downstream + CHH_rather_met_downstream + CHH_met_downstream

# calculation of proportions of CG, ... in unmet, rather_unmet, ... -> (1) + (3)
CG_unmet = CG_unmet_CDS + CG_unmet_UTR + CG_unmet_intron + CG_unmet_non_genic + CG_unmet_upstream + CG_unmet_downstream
CHG_unmet = CHG_unmet_CDS + CHG_unmet_UTR + CHG_unmet_intron + CHG_unmet_non_genic + CHG_unmet_upstream + CHG_unmet_downstream
CHH_unmet = CHH_unmet_CDS + CHH_unmet_UTR + CHH_unmet_intron + CHH_unmet_non_genic + CHH_unmet_upstream + CHH_unmet_downstream
CG_rather_unmet = CG_rather_unmet_CDS + CG_rather_unmet_UTR + CG_rather_unmet_intron + CG_rather_unmet_non_genic + CG_rather_unmet_upstream + CG_rather_unmet_downstream
CHG_rather_unmet = CHG_rather_unmet_CDS + CHG_rather_unmet_UTR + CHG_rather_unmet_intron + CHG_rather_unmet_non_genic + CHG_rather_unmet_upstream + CHG_rather_unmet_downstream
CHH_rather_unmet = CHH_rather_unmet_CDS + CHH_rather_unmet_UTR + CHH_rather_unmet_intron + CHH_rather_unmet_non_genic + CHH_rather_unmet_upstream + CHH_rather_unmet_downstream
CG_50met = CG_50met_CDS + CG_50met_UTR + CG_50met_intron + CG_50met_non_genic + CG_50met_upstream + CG_50met_downstream
CHG_50met = CHG_50met_CDS + CHG_50met_UTR + CHG_50met_intron + CHG_50met_non_genic + CHG_50met_upstream + CHG_50met_downstream
CHH_50met = CHH_50met_CDS + CHH_50met_UTR + CHH_50met_intron + CHH_50met_non_genic + CHH_50met_upstream + CHH_50met_downstream
CG_rather_met = CG_rather_met_CDS + CG_rather_met_UTR + CG_rather_met_intron + CG_rather_met_non_genic + CG_rather_met_upstream + CG_rather_met_downstream
CHG_rather_met = CHG_rather_met_CDS + CHG_rather_met_UTR + CHG_rather_met_intron + CHG_rather_met_non_genic + CHG_rather_met_upstream + CHG_rather_met_downstream
CHH_rather_met = CHH_rather_met_CDS + CHH_rather_met_UTR + CHH_rather_met_intron + CHH_rather_met_non_genic + CHH_rather_met_upstream + CHH_rather_met_downstream
CG_met = CG_met_CDS + CG_met_UTR + CG_met_intron + CG_met_non_genic + CG_met_upstream + CG_met_downstream
CHG_met = CHG_met_CDS + CHG_met_UTR + CHG_met_intron + CHG_met_non_genic + CHG_met_upstream + CHG_met_downstream
CHH_met = CHH_met_CDS + CHH_met_UTR + CHH_met_intron + CHH_met_non_genic + CHH_met_upstream + CHH_met_downstream

# calculation of proportions of unmet, ... CDS, ... -> (2) + (3) 
c_unmet_CDS = CG_unmet_CDS + CHG_unmet_CDS + CHH_unmet_CDS
c_unmet_UTR = CG_unmet_UTR + CHG_unmet_UTR + CHH_unmet_UTR
c_unmet_exon = c_unmet_CDS + c_unmet_UTR
c_unmet_intron = CG_unmet_intron + CHG_unmet_intron + CHH_unmet_intron
c_unmet_genic = c_unmet_exon + c_unmet_intron
c_unmet_non_genic = CG_unmet_non_genic + CHG_unmet_non_genic + CHH_unmet_non_genic
c_unmet_upstream = CG_unmet_upstream + CHG_unmet_upstream + CHH_unmet_upstream
c_unmet_downstream = CG_unmet_downstream + CHG_unmet_downstream + CHH_unmet_downstream
c_rather_unmet_CDS = CG_rather_unmet_CDS + CHG_rather_unmet_CDS + CHH_rather_unmet_CDS
c_rather_unmet_UTR = CG_rather_unmet_UTR + CHG_rather_unmet_UTR + CHH_rather_unmet_UTR
c_rather_unmet_exon = c_rather_unmet_CDS + c_rather_unmet_UTR
c_rather_unmet_intron = CG_rather_unmet_intron + CHG_rather_unmet_intron + CHH_rather_unmet_intron
c_rather_unmet_genic = c_rather_unmet_exon + c_rather_unmet_intron
c_rather_unmet_non_genic = CG_rather_unmet_non_genic + CHG_rather_unmet_non_genic + CHH_rather_unmet_non_genic
c_rather_unmet_upstream = CG_rather_unmet_upstream + CHG_rather_unmet_upstream + CHH_rather_unmet_upstream
c_rather_unmet_downstream = CG_rather_unmet_downstream + CHG_rather_unmet_downstream + CHH_rather_unmet_downstream
c_50met_CDS = CG_50met_CDS + CHG_50met_CDS + CHH_50met_CDS
c_50met_UTR = CG_50met_UTR + CHG_50met_UTR + CHH_50met_UTR
c_50met_exon = c_50met_CDS + c_50met_UTR
c_50met_intron = CG_50met_intron + CHG_50met_intron + CHH_50met_intron
c_50met_genic = c_50met_exon + c_50met_intron
c_50met_non_genic = CG_50met_non_genic + CHG_50met_non_genic + CHH_50met_non_genic
c_50met_upstream = CG_50met_upstream + CHG_50met_upstream + CHH_50met_upstream
c_50met_downstream = CG_50met_downstream + CHG_50met_downstream + CHH_50met_downstream
c_rather_met_CDS = CG_rather_met_CDS + CHG_rather_met_CDS + CHH_rather_met_CDS
c_rather_met_UTR = CG_rather_met_UTR + CHG_rather_met_UTR + CHH_rather_met_UTR
c_rather_met_exon = c_rather_met_CDS + c_rather_met_UTR
c_rather_met_intron = CG_rather_met_intron + CHG_rather_met_intron + CHH_rather_met_intron
c_rather_met_genic = c_rather_met_exon + c_rather_met_intron
c_rather_met_non_genic = CG_rather_met_non_genic + CHG_rather_met_non_genic + CHH_rather_met_non_genic
c_rather_met_upstream = CG_rather_met_upstream + CHG_rather_met_upstream + CHH_rather_met_upstream
c_rather_met_downstream = CG_rather_met_downstream + CHG_rather_met_downstream + CHH_rather_met_downstream
c_met_CDS = CG_met_CDS + CHG_met_CDS + CHH_met_CDS
c_met_UTR = CG_met_UTR + CHG_met_UTR + CHH_met_UTR
c_met_exon = c_met_CDS + c_met_UTR
c_met_intron = CG_met_intron + CHG_met_intron + CHH_met_intron
c_met_genic = c_met_exon + c_met_intron
c_met_non_genic = CG_met_non_genic + CHG_met_non_genic + CHH_met_non_genic
c_met_upstream = CG_met_upstream + CHG_met_upstream + CHH_met_upstream
c_met_downstream = CG_met_downstream + CHG_met_downstream + CHH_met_downstream

# calculate proportion of CG, CHG and CHH context genome wide -> (1)
CG_fraction = CG_unmet + CG_rather_unmet + CG_50met + CG_rather_met + CG_met
CHG_fraction = CHG_unmet + CHG_rather_unmet + CHG_50met + CHG_rather_met + CHG_met
CHH_fraction = CHH_unmet + CHH_rather_unmet + CHH_50met + CHH_rather_met + CHH_met

# calculate proportion of CDS, ... genome wide  -> (2)
c_CDS = CG_CDS + CHG_CDS + CHH_CDS
c_UTR =  CG_UTR + CHG_UTR + CHH_UTR
c_exon = CG_exon + CHG_exon + CHH_exon
c_intron = CG_intron + CHG_intron + CHH_intron
c_genic = CG_genic + CHG_genic + CHH_genic
c_non_genic = CG_non_genic + CHG_non_genic + CHH_non_genic
c_upstream = CG_upstream + CHG_upstream + CHH_upstream
c_downstream = CG_downstream + CHG_downstream + CHH_downstream

# calculate proportion of unmet, ... genome wide -> (3)
c_unmet = CG_unmet + CHG_unmet + CHH_unmet
c_rather_unmet = CG_rather_unmet + CHG_rather_unmet + CHH_rather_unmet
c_50met = CG_50met + CHG_50met + CHH_50met
c_rather_met = CG_rather_met + CHG_rather_met + CHH_rather_met
c_met = CG_met + CHG_met + CHH_met

print('-> values of interest were calculated')

############################################################################################################################

# calculation for proportion of meth, ... in each context for output (in relation to context) -> (1) + (2) + (3) / (1) + (2) / (1) + (3)
if CG_fraction != 0:
    # (1) + (2) + (3)
    CG_100_unmet_CDS = CG_unmet_CDS / CG_fraction * 100
    CG_100_unmet_UTR = CG_unmet_UTR / CG_fraction * 100
    CG_100_unmet_exon = CG_100_unmet_CDS + CG_100_unmet_UTR
    CG_100_unmet_intron = CG_unmet_intron  / CG_fraction * 100
    CG_100_unmet_genic = CG_100_unmet_exon + CG_100_unmet_intron
    CG_100_unmet_non_genic = CG_unmet_non_genic / CG_fraction * 100
    CG_100_unmet_upstream = CG_unmet_upstream / CG_fraction * 100
    CG_100_unmet_downstream = CG_unmet_downstream / CG_fraction * 100
    CG_100_rather_unmet_CDS = CG_rather_unmet_CDS / CG_fraction * 100
    CG_100_rather_unmet_UTR = CG_rather_unmet_UTR / CG_fraction * 100
    CG_100_rather_unmet_exon = CG_100_rather_unmet_CDS + CG_100_rather_unmet_UTR
    CG_100_rather_unmet_intron = CG_rather_unmet_intron  / CG_fraction * 100
    CG_100_rather_unmet_genic = CG_100_rather_unmet_exon + CG_100_rather_unmet_intron
    CG_100_rather_unmet_non_genic = CG_rather_unmet_non_genic / CG_fraction * 100
    CG_100_rather_unmet_upstream = CG_rather_unmet_upstream / CG_fraction * 100
    CG_100_rather_unmet_downstream = CG_rather_unmet_downstream / CG_fraction * 100
    CG_100_50met_CDS = CG_50met_CDS / CG_fraction * 100
    CG_100_50met_UTR = CG_50met_UTR / CG_fraction * 100
    CG_100_50met_exon = CG_100_50met_CDS + CG_100_50met_UTR
    CG_100_50met_intron = CG_50met_intron  / CG_fraction * 100
    CG_100_50met_genic = CG_100_50met_exon + CG_100_50met_intron
    CG_100_50met_non_genic = CG_50met_non_genic / CG_fraction * 100
    CG_100_50met_upstream = CG_50met_upstream / CG_fraction * 100
    CG_100_50met_downstream = CG_50met_downstream / CG_fraction * 100
    CG_100_rather_met_CDS = CG_rather_met_CDS / CG_fraction * 100
    CG_100_rather_met_UTR = CG_rather_met_UTR / CG_fraction * 100
    CG_100_rather_met_exon = CG_100_rather_met_CDS + CG_100_rather_met_UTR
    CG_100_rather_met_intron = CG_rather_met_intron  / CG_fraction * 100
    CG_100_rather_met_genic = CG_100_rather_met_exon + CG_100_rather_met_intron
    CG_100_rather_met_non_genic = CG_rather_met_non_genic / CG_fraction * 100
    CG_100_rather_met_upstream = CG_rather_met_upstream / CG_fraction * 100
    CG_100_rather_met_downstream = CG_rather_met_downstream / CG_fraction * 100
    CG_100_met_CDS = CG_met_CDS / CG_fraction * 100
    CG_100_met_UTR = CG_met_UTR / CG_fraction * 100
    CG_100_met_exon = CG_100_met_CDS + CG_100_met_UTR
    CG_100_met_intron = CG_met_intron  / CG_fraction * 100
    CG_100_met_genic = CG_100_met_exon + CG_100_met_intron
    CG_100_met_non_genic = CG_met_non_genic / CG_fraction * 100
    CG_100_met_upstream = CG_met_upstream / CG_fraction * 100
    CG_100_met_downstream = CG_met_downstream / CG_fraction * 100
    # (1) + (2)
    CG_100_CDS = CG_CDS / CG_fraction * 100
    CG_100_UTR = CG_UTR / CG_fraction * 100
    CG_100_exon = CG_exon / CG_fraction * 100
    CG_100_intron = CG_intron / CG_fraction * 100
    CG_100_genic = CG_genic / CG_fraction * 100
    CG_100_non_genic = CG_non_genic / CG_fraction * 100
    CG_100_upstream = CG_upstream / CG_fraction * 100
    CG_100_downstream = CG_downstream / CG_fraction * 100
    # (1) + (3)
    CG_100_unmet = CG_unmet / CG_fraction * 100
    CG_100_rather_unmet = CG_rather_unmet / CG_fraction * 100
    CG_100_50met = CG_50met / CG_fraction * 100
    CG_100_rather_met = CG_rather_met / CG_fraction * 100
    CG_100_met = CG_met / CG_fraction * 100
if CHG_fraction != 0:
    # (1) + (2) + (3)
    CHG_100_unmet_CDS = CHG_unmet_CDS / CHG_fraction * 100
    CHG_100_unmet_UTR = CHG_unmet_UTR / CHG_fraction * 100
    CHG_100_unmet_exon = CHG_100_unmet_CDS + CHG_100_unmet_UTR
    CHG_100_unmet_intron = CHG_unmet_intron  / CHG_fraction * 100
    CHG_100_unmet_genic = CHG_100_unmet_exon + CHG_100_unmet_intron
    CHG_100_unmet_non_genic = CHG_unmet_non_genic / CHG_fraction * 100
    CHG_100_unmet_upstream = CHG_unmet_upstream / CHG_fraction * 100
    CHG_100_unmet_downstream = CHG_unmet_downstream / CHG_fraction * 100
    CHG_100_rather_unmet_CDS = CHG_rather_unmet_CDS / CHG_fraction * 100
    CHG_100_rather_unmet_UTR = CHG_rather_unmet_UTR / CHG_fraction * 100
    CHG_100_rather_unmet_exon = CHG_100_rather_unmet_CDS + CHG_100_rather_unmet_UTR
    CHG_100_rather_unmet_intron = CHG_rather_unmet_intron  / CHG_fraction * 100
    CHG_100_rather_unmet_genic = CHG_100_rather_unmet_exon + CHG_100_rather_unmet_intron
    CHG_100_rather_unmet_non_genic = CHG_rather_unmet_non_genic / CHG_fraction * 100
    CHG_100_rather_unmet_upstream = CHG_rather_unmet_upstream / CHG_fraction * 100
    CHG_100_rather_unmet_downstream = CHG_rather_unmet_downstream / CHG_fraction * 100
    CHG_100_50met_CDS = CHG_50met_CDS / CHG_fraction * 100
    CHG_100_50met_UTR = CHG_50met_UTR / CHG_fraction * 100
    CHG_100_50met_exon = CHG_100_50met_CDS + CHG_100_50met_UTR
    CHG_100_50met_intron = CHG_50met_intron  / CHG_fraction * 100
    CHG_100_50met_genic = CHG_100_50met_exon + CHG_100_50met_intron
    CHG_100_50met_non_genic = CHG_50met_non_genic / CHG_fraction * 100
    CHG_100_50met_upstream = CHG_50met_upstream / CHG_fraction * 100
    CHG_100_50met_downstream = CHG_50met_downstream / CHG_fraction * 100
    CHG_100_rather_met_CDS = CHG_rather_met_CDS / CHG_fraction * 100
    CHG_100_rather_met_UTR = CHG_rather_met_UTR / CHG_fraction * 100
    CHG_100_rather_met_exon = CHG_100_rather_met_CDS + CHG_100_rather_met_UTR
    CHG_100_rather_met_intron = CHG_rather_met_intron  / CHG_fraction * 100
    CHG_100_rather_met_genic = CHG_100_rather_met_exon + CHG_100_rather_met_intron
    CHG_100_rather_met_non_genic = CHG_rather_met_non_genic / CHG_fraction * 100
    CHG_100_rather_met_upstream = CHG_rather_met_upstream / CHG_fraction * 100
    CHG_100_rather_met_downstream = CHG_rather_met_downstream / CHG_fraction * 100
    CHG_100_met_CDS = CHG_met_CDS / CHG_fraction * 100
    CHG_100_met_UTR = CHG_met_UTR / CHG_fraction * 100
    CHG_100_met_exon = CHG_100_met_CDS + CHG_100_met_UTR
    CHG_100_met_intron = CHG_met_intron  / CHG_fraction * 100
    CHG_100_met_genic = CHG_100_met_exon + CHG_100_met_intron
    CHG_100_met_non_genic = CHG_met_non_genic / CHG_fraction * 100
    CHG_100_met_upstream = CHG_met_upstream / CHG_fraction * 100
    CHG_100_met_downstream = CHG_met_downstream / CHG_fraction * 100
    # (1) + (2)
    CHG_100_CDS = CHG_CDS / CHG_fraction * 100
    CHG_100_UTR = CHG_UTR / CHG_fraction * 100
    CHG_100_exon = CHG_exon / CHG_fraction * 100
    CHG_100_intron = CHG_intron / CHG_fraction * 100
    CHG_100_genic = CHG_genic / CHG_fraction * 100
    CHG_100_non_genic = CHG_non_genic / CHG_fraction * 100
    CHG_100_upstream = CHG_upstream / CHG_fraction * 100
    CHG_100_downstream = CHG_downstream / CHG_fraction * 100
    # (1) + (3)
    CHG_100_unmet = CHG_unmet / CHG_fraction * 100
    CHG_100_rather_unmet = CHG_rather_unmet / CHG_fraction * 100
    CHG_100_50met = CHG_50met / CHG_fraction * 100
    CHG_100_rather_met = CHG_rather_met / CHG_fraction * 100
    CHG_100_met = CHG_met / CHG_fraction * 100
if CHH_fraction != 0:
    # (1) + (2) + (3)
    CHH_100_unmet_CDS = CHH_unmet_CDS / CHH_fraction * 100
    CHH_100_unmet_UTR = CHH_unmet_UTR / CHH_fraction * 100
    CHH_100_unmet_exon = CHH_100_unmet_CDS + CHH_100_unmet_UTR
    CHH_100_unmet_intron = CHH_unmet_intron  / CHH_fraction * 100
    CHH_100_unmet_genic = CHH_100_unmet_exon + CHH_100_unmet_intron
    CHH_100_unmet_non_genic = CHH_unmet_non_genic / CHH_fraction * 100
    CHH_100_unmet_upstream = CHH_unmet_upstream / CHH_fraction * 100
    CHH_100_unmet_downstream = CHH_unmet_downstream / CHH_fraction * 100
    CHH_100_rather_unmet_CDS = CHH_rather_unmet_CDS / CHH_fraction * 100
    CHH_100_rather_unmet_UTR = CHH_rather_unmet_UTR / CHH_fraction * 100
    CHH_100_rather_unmet_exon = CHH_100_rather_unmet_CDS + CHH_100_rather_unmet_UTR
    CHH_100_rather_unmet_intron = CHH_rather_unmet_intron  / CHH_fraction * 100
    CHH_100_rather_unmet_genic = CHH_100_rather_unmet_exon + CHH_100_rather_unmet_intron
    CHH_100_rather_unmet_non_genic = CHH_rather_unmet_non_genic / CHH_fraction * 100
    CHH_100_rather_unmet_upstream = CHH_rather_unmet_upstream / CHH_fraction * 100
    CHH_100_rather_unmet_downstream = CHH_rather_unmet_downstream / CHH_fraction * 100
    CHH_100_50met_CDS = CHH_50met_CDS / CHH_fraction * 100
    CHH_100_50met_UTR = CHH_50met_UTR / CHH_fraction * 100
    CHH_100_50met_exon = CHH_100_50met_CDS + CHH_100_50met_UTR
    CHH_100_50met_intron = CHH_50met_intron  / CHH_fraction * 100
    CHH_100_50met_genic = CHH_100_50met_exon + CHH_100_50met_intron
    CHH_100_50met_non_genic = CHH_50met_non_genic / CHH_fraction * 100
    CHH_100_50met_upstream = CHH_50met_upstream / CHH_fraction * 100
    CHH_100_50met_downstream = CHH_50met_downstream / CHH_fraction * 100
    CHH_100_rather_met_CDS = CHH_rather_met_CDS / CHH_fraction * 100
    CHH_100_rather_met_UTR = CHH_rather_met_UTR / CHH_fraction * 100
    CHH_100_rather_met_exon = CHH_100_rather_met_CDS + CHH_100_rather_met_UTR
    CHH_100_rather_met_intron = CHH_rather_met_intron  / CHH_fraction * 100
    CHH_100_rather_met_genic = CHH_100_rather_met_exon + CHH_100_rather_met_intron
    CHH_100_rather_met_non_genic = CHH_rather_met_non_genic / CHH_fraction * 100
    CHH_100_rather_met_upstream = CHH_rather_met_upstream / CHH_fraction * 100
    CHH_100_rather_met_downstream = CHH_rather_met_downstream / CHH_fraction * 100
    CHH_100_met_CDS = CHH_met_CDS / CHH_fraction * 100
    CHH_100_met_UTR = CHH_met_UTR / CHH_fraction * 100
    CHH_100_met_exon = CHH_100_met_CDS + CHH_100_met_UTR
    CHH_100_met_intron = CHH_met_intron  / CHH_fraction * 100
    CHH_100_met_genic = CHH_100_met_exon + CHH_100_met_intron
    CHH_100_met_non_genic = CHH_met_non_genic / CHH_fraction * 100
    CHH_100_met_upstream = CHH_met_upstream / CHH_fraction * 100
    CHH_100_met_downstream = CHH_met_downstream / CHH_fraction * 100
    # (1) + (2)
    CHH_100_CDS = CHH_CDS / CHH_fraction * 100
    CHH_100_UTR = CHH_UTR / CHH_fraction * 100
    CHH_100_exon = CHH_exon / CHH_fraction * 100
    CHH_100_intron = CHH_intron / CHH_fraction * 100
    CHH_100_genic = CHH_genic / CHH_fraction * 100
    CHH_100_non_genic = CHH_non_genic / CHH_fraction * 100
    CHH_100_upstream = CHH_upstream / CHH_fraction * 100
    CHH_100_downstream = CHH_downstream / CHH_fraction * 100
    # (1) + (3)
    CHH_100_unmet = CHH_unmet / CHH_fraction * 100
    CHH_100_rather_unmet = CHH_rather_unmet / CHH_fraction * 100
    CHH_100_50met = CHH_50met / CHH_fraction * 100
    CHH_100_rather_met = CHH_rather_met / CHH_fraction * 100
    CHH_100_met = CHH_met / CHH_fraction * 100

print('-> calculations for context output were generated')

############################################################################################################################

# calculate realative frequencie & calculate values of interest -> chr_tigss
for chr_tigs in chr_tigs_distribution:
    # calculate realtive freqneuncies
    chr_tigs[2] = chr_tigs[1]/c_counter*100
    chr_tigs[4] = chr_tigs[3]/c_counter*100
    chr_tigs[6] = chr_tigs[5]/c_counter*100
    chr_tigs[8] = chr_tigs[7]/c_counter*100
    chr_tigs[10] = chr_tigs[9]/c_counter*100
    chr_tigs[12] = chr_tigs[11]/c_counter*100
    chr_tigs[14] = chr_tigs[13]/c_counter*100
    chr_tigs[16] = chr_tigs[15]/c_counter*100
    chr_tigs[18] = chr_tigs[17]/c_counter*100
    chr_tigs[20] = chr_tigs[19]/c_counter*100
    chr_tigs[22] = chr_tigs[21]/c_counter*100
    chr_tigs[24] = chr_tigs[23]/c_counter*100
    chr_tigs[26] = chr_tigs[25]/c_counter*100
    chr_tigs[28] = chr_tigs[27]/c_counter*100
    chr_tigs[30] = chr_tigs[29]/c_counter*100

print('-> calculations for chromosome/contig output were generated')

############################################################################################################################

# write the results in a file
with open(fileLOG_3_1, "w") as out3_1:
    # return values of interest
    out3_1.write('-> tool\n')
    out3_1.write(__file__ + '\n')
    out3_1.write('\n')

    out3_1.write('-> description\n')
    out3_1.write('in the third of three steps the data counts are merged and result values of interst are calculated\n')
    out3_1.write('\n')

    out3_1.write('-> 5mer specific count inputs\n')
    for counts in fileCounts5mer_list:
        out3_1.write(counts + '\n')
    out3_1.write('\n')

    out3_1.write('-> chromosome/contig specific count inputs\n')
    for counts in fileCountsChr_list:
        out3_1.write(counts + '\n')
    out3_1.write('\n')
    
    out3_1.write('-> log and output\n')
    out3_1.write(fileLOG_3_1 + '\n')
    out3_1.write('\n')

    out3_1.write('------------------------------------------------------------------------------\n')
    out3_1.write('\n')

    out3_1.write('-> note\n')
    out3_1.write('UTR: means all cytosines, that are in a CDS but not in an exon \n')
    out3_1.write('non genic: means all cytosines, that are not related to a gene (not in or within a 1000 bp zone around a gene) \n')
    out3_1.write('to calculate all cytosines that are not in a genes add all upstream and downstream Cs to non genic\n')
    out3_1.write('\n')

    out3_1.write('------------------------------------------------------------------------------\n')
    out3_1.write('\n')

    out3_1.write('-> analysed cytosines\n')
    out3_1.write('\n')
    out3_1.write('total number of cytosines: ' + str(c_counter) + '\n')
    out3_1.write('\n')
    out3_1.write('highly unmethylated cytosines: ' + str(round(c_unmet, 2)) + ' %\n')
    out3_1.write('rather unmethylated cytosines: ' + str(round(c_rather_unmet, 2)) + ' %\n')
    out3_1.write('50/50-methylated cytosines: ' + str(round(c_50met, 2)) + ' %\n')
    out3_1.write('rather methylated cytosines: ' + str(round(c_rather_met, 2)) + ' %\n')
    out3_1.write('highly methylated cytosines: ' + str(round(c_met, 2)) + ' %\n')
    out3_1.write('\n')

    out3_1.write('cytosines in CDSs: ' + str(round(c_CDS, 2)) + ' %\n')
    out3_1.write('highly unmethylated cytosines in CDSs: ' + str(round(c_unmet_CDS, 2)) + ' %\n')
    out3_1.write('rather unmethylated cytosines in CDSs: ' + str(round(c_rather_unmet_CDS, 2)) + ' %\n')
    out3_1.write('50/50-methylated cytosines in CDSs: ' + str(round(c_50met_CDS, 2)) + ' %\n')
    out3_1.write('rather methylated cytosines in CDSs: ' + str(round(c_rather_met_CDS, 2)) + ' %\n')
    out3_1.write('highly methylated cytosines in CDSs: ' + str(round(c_met_CDS, 2)) + ' %\n')
    out3_1.write('\n')
    out3_1.write('cytosines in UTRs: ' + str(round(c_UTR, 2)) + ' %\n')
    out3_1.write('highly unmethylated cytosines in UTRs: ' + str(round(c_unmet_UTR, 2)) + ' %\n')
    out3_1.write('rather unmethylated cytosines in UTRs: ' + str(round(c_rather_unmet_UTR, 2)) + ' %\n')
    out3_1.write('50/50-methylated cytosines in UTRs: ' + str(round(c_50met_UTR, 2)) + ' %\n')
    out3_1.write('rather methylated cytosines in UTRs: ' + str(round(c_rather_met_UTR, 2)) + ' %\n')
    out3_1.write('highly methylated cytosines in UTRs: ' + str(round(c_met_UTR, 2)) + ' %\n')
    out3_1.write('\n')
    out3_1.write('cytosines in exons: ' + str(round(c_exon, 2)) + ' %\n')
    out3_1.write('highly unmethylated cytosines in exons: ' + str(round(c_unmet_exon, 2)) + ' %\n')
    out3_1.write('rather unmethylated cytosines in exons: ' + str(round(c_rather_unmet_exon, 2)) + ' %\n')
    out3_1.write('50/50-methylated cytosines in exons: ' + str(round(c_50met_exon, 2)) + ' %\n')
    out3_1.write('rather methylated cytosines in exons: ' + str(round(c_rather_met_exon, 2)) + ' %\n')
    out3_1.write('highly methylated cytosines in exons: ' + str(round(c_met_exon, 2)) + ' %\n')
    out3_1.write('\n')
    out3_1.write('cytosines in introns: ' + str(round(c_intron, 2)) + ' %\n')
    out3_1.write('highly unmethylated cytosines in introns: ' + str(round(c_unmet_intron, 2)) + ' %\n')
    out3_1.write('rather unmethylated cytosines in introns: ' + str(round(c_rather_unmet_intron, 2)) + ' %\n')
    out3_1.write('50/50-methylated cytosines in introns: ' + str(round(c_50met_intron, 2)) + ' %\n')
    out3_1.write('rather methylated cytosines in introns: ' + str(round(c_rather_met_intron, 2)) + ' %\n')
    out3_1.write('highly methylated cytosines in introns: ' + str(round(c_met_intron, 2)) + ' %\n')
    out3_1.write('\n')
    out3_1.write('cytosines in genic regions: ' + str(round(c_genic, 2)) + ' %\n')
    out3_1.write('highly unmethylated cytosines in genic regions: ' + str(round(c_unmet_genic, 2)) + ' %\n')
    out3_1.write('rather unmethylated cytosines in genic regions: ' + str(round(c_rather_unmet_genic, 2)) + ' %\n')
    out3_1.write('50/50-methylated cytosines in genic regions: ' + str(round(c_50met_genic, 2)) + ' %\n')
    out3_1.write('rather methylated cytosines in genic regions: ' + str(round(c_rather_met_genic, 2)) + ' %\n')
    out3_1.write('highly methylated cytosines in genic regions: ' + str(round(c_met_genic, 2)) + ' %\n')
    out3_1.write('\n')
    out3_1.write('cytosines in 1kb upstream: ' + str(round(c_upstream, 2)) + ' %\n')
    out3_1.write('highly unmethylated cytosines in 1kb upstream: ' + str(round(c_unmet_upstream, 2)) + ' %\n')
    out3_1.write('rather unmethylated cytosines in 1kb upstream: ' + str(round(c_rather_unmet_upstream, 2)) + ' %\n')
    out3_1.write('50/50-methylated cytosines in 1kb upstream: ' + str(round(c_50met_upstream, 2)) + ' %\n')
    out3_1.write('rather methylated cytosines in 1kb upstream: ' + str(round(c_rather_met_upstream, 2)) + ' %\n')
    out3_1.write('highly methylated cytosines in 1kb upstream: ' + str(round(c_met_upstream, 2)) + ' %\n')
    out3_1.write('\n')
    out3_1.write('cytosines in 1kb downstream: ' + str(round(c_downstream, 2)) + ' %\n')
    out3_1.write('highly unmethylated cytosines in 1kb downstream: ' + str(round(c_unmet_downstream, 2)) + ' %\n')
    out3_1.write('rather unmethylated cytosines in 1kb downstream: ' + str(round(c_rather_unmet_downstream, 2)) + ' %\n')
    out3_1.write('50/50-methylated cytosines in 1kb downstream: ' + str(round(c_50met_downstream, 2)) + ' %\n')
    out3_1.write('rather methylated cytosines in 1kb downstream: ' + str(round(c_rather_met_downstream, 2)) + ' %\n')
    out3_1.write('highly methylated cytosines in 1kb downstream: ' + str(round(c_met_downstream, 2)) + ' %\n')
    out3_1.write('\n')
    out3_1.write('cytosines in non genic regions: ' + str(round(c_non_genic, 2)) + ' %\n')
    out3_1.write('highly unmethylated cytosines in non genic regions: ' + str(round(c_unmet_non_genic, 2)) + ' %\n')
    out3_1.write('rather unmethylated cytosines in non genic regions: ' + str(round(c_rather_unmet_non_genic, 2)) + ' %\n')
    out3_1.write('50/50-methylated cytosines in non genic regions: ' + str(round(c_50met_non_genic, 2)) + ' %\n')
    out3_1.write('rather methylated cytosines in non genic regions: ' + str(round(c_rather_met_non_genic, 2)) + ' %\n')
    out3_1.write('highly methylated cytosines in non genic regions: ' + str(round(c_met_non_genic, 2)) + ' %\n')
    out3_1.write('\n')

    out3_1.write('------------------------------------------------------------------------------\n')
    out3_1.write('\n')
    
    out3_1.write('-> CG\n')
    out3_1.write('\n')
    out3_1.write('CGs: ' + str(round(CG_fraction, 2)) + ' %\n')
    out3_1.write('\n')
    out3_1.write('highly unmethylated CGs: ' + str(round(CG_100_unmet, 2)) + ' %\n')
    out3_1.write('rather unmethylated CGs: ' + str(round(CG_100_rather_unmet, 2)) + ' %\n')
    out3_1.write('50/50-methylated CGs: ' + str(round(CG_100_50met, 2)) + ' %\n')
    out3_1.write('rather methylated CGs: ' + str(round(CG_100_rather_met, 2)) + ' %\n')
    out3_1.write('highly methylated CGs: ' + str(round(CG_100_met, 2)) + ' %\n')
    out3_1.write('\n')
    out3_1.write('CGs in CDSs: ' + str(round(CG_100_CDS, 2)) + ' %\n')
    out3_1.write('highly unmethylated CGs in CDSs: ' + str(round(CG_100_unmet_CDS, 2)) + ' %\n')
    out3_1.write('rather unmethylated CGs in CDSs: ' + str(round(CG_100_rather_unmet_CDS, 2)) + ' %\n')
    out3_1.write('50/50-methylated CGs in CDSs: ' + str(round(CG_100_50met_CDS, 2)) + ' %\n')
    out3_1.write('rather methylated CGs in CDSs: ' + str(round(CG_100_rather_met_CDS, 2)) + ' %\n')
    out3_1.write('highly methylated CGs in CDSs: ' + str(round(CG_100_met_CDS, 2)) + ' %\n')
    out3_1.write('\n')
    out3_1.write('CGs in UTRs: ' + str(round(CG_100_UTR, 2)) + ' %\n')
    out3_1.write('highly unmethylated CGs in UTRs: ' + str(round(CG_100_unmet_UTR, 2)) + ' %\n')
    out3_1.write('rather unmethylated CGs in UTRs: ' + str(round(CG_100_rather_unmet_UTR, 2)) + ' %\n')
    out3_1.write('50/50-methylated CGs in UTRs: ' + str(round(CG_100_50met_UTR, 2)) + ' %\n')
    out3_1.write('rather methylated CGs in UTRs: ' + str(round(CG_100_rather_met_UTR, 2)) + ' %\n')
    out3_1.write('highly methylated CGs in UTRs: ' + str(round(CG_100_met_UTR, 2)) + ' %\n')
    out3_1.write('\n')
    out3_1.write('CGs in exons: ' + str(round(CG_100_exon, 2)) + ' %\n')
    out3_1.write('highly unmethylated CGs in exons: ' + str(round(CG_100_unmet_exon, 2)) + ' %\n')
    out3_1.write('rather unmethylated CGs in exons: ' + str(round(CG_100_rather_unmet_exon, 2)) + ' %\n')
    out3_1.write('50/50-methylated CGs in exons: ' + str(round(CG_100_50met_exon, 2)) + ' %\n')
    out3_1.write('rather methylated CGs in exons: ' + str(round(CG_100_rather_met_exon, 2)) + ' %\n')
    out3_1.write('highly methylated CGs in exons: ' + str(round(CG_100_met_exon, 2)) + ' %\n')
    out3_1.write('\n')
    out3_1.write('CGs in introns: ' + str(round(CG_100_intron, 2)) + ' %\n')
    out3_1.write('highly unmethylated CGs in introns: ' + str(round(CG_100_unmet_intron, 2)) + ' %\n')
    out3_1.write('rather unmethylated CGs in introns: ' + str(round(CG_100_rather_unmet_intron, 2)) + ' %\n')
    out3_1.write('50/50-methylated CGs in introns: ' + str(round(CG_100_50met_intron, 2)) + ' %\n')
    out3_1.write('rather methylated CGs in introns: ' + str(round(CG_100_rather_met_intron, 2)) + ' %\n')
    out3_1.write('highly methylated CGs in introns: ' + str(round(CG_100_met_intron, 2)) + ' %\n')
    out3_1.write('\n')
    out3_1.write('CGs in genic regions: ' + str(round(CG_100_genic, 2)) + ' %\n')
    out3_1.write('highly unmethylated CGs in genic regions: ' + str(round(CG_100_unmet_genic, 2)) + ' %\n')
    out3_1.write('rather unmethylated CGs in genic regions: ' + str(round(CG_100_rather_unmet_genic, 2)) + ' %\n')
    out3_1.write('50/50-methylated CGs in genic regions: ' + str(round(CG_100_50met_genic, 2)) + ' %\n')
    out3_1.write('rather methylated CGs in genic regions: ' + str(round(CG_100_rather_met_genic, 2)) + ' %\n')
    out3_1.write('highly methylated CGs in genic regions: ' + str(round(CG_100_met_genic, 2)) + ' %\n')
    out3_1.write('\n')
    out3_1.write('CGs in 1kb upstream: ' + str(round(CG_100_upstream, 2)) + ' %\n')
    out3_1.write('highly unmethylated CGs in 1kb upstream: ' + str(round(CG_100_unmet_upstream, 2)) + ' %\n')
    out3_1.write('rather unmethylated CGs in 1kb upstream: ' + str(round(CG_100_rather_unmet_upstream, 2)) + ' %\n')
    out3_1.write('50/50-methylated CGs in 1kb upstream: ' + str(round(CG_100_50met_upstream, 2)) + ' %\n')
    out3_1.write('rather methylated CGs in 1kb upstream: ' + str(round(CG_100_rather_met_upstream, 2)) + ' %\n')
    out3_1.write('highly methylated CGs in 1kb upstream: ' + str(round(CG_100_met_upstream, 2)) + ' %\n')
    out3_1.write('\n')
    out3_1.write('CGs in 1kb downstream: ' + str(round(CG_100_downstream, 2)) + ' %\n')
    out3_1.write('highly unmethylated CGs in 1kb downstream: ' + str(round(CG_100_unmet_downstream, 2)) + ' %\n')
    out3_1.write('rather unmethylated CGs in 1kb downstream: ' + str(round(CG_100_rather_unmet_downstream, 2)) + ' %\n')
    out3_1.write('50/50-methylated CGs in 1kb downstream: ' + str(round(CG_100_50met_downstream, 2)) + ' %\n')
    out3_1.write('rather methylated CGs in 1kb downstream: ' + str(round(CG_100_rather_met_downstream, 2)) + ' %\n')
    out3_1.write('highly methylated CGs in 1kb downstream: ' + str(round(CG_100_met_downstream, 2)) + ' %\n')
    out3_1.write('\n')
    out3_1.write('CGs in non genic regions: ' + str(round(CG_100_non_genic, 2)) + ' %\n')
    out3_1.write('highly unmethylated CGs in non genic regions: ' + str(round(CG_100_unmet_non_genic, 2)) + ' %\n')
    out3_1.write('rather unmethylated CGs in non genic regions: ' + str(round(CG_100_rather_unmet_non_genic, 2)) + ' %\n')
    out3_1.write('50/50-methylated CGs in non genic regions: ' + str(round(CG_100_50met_non_genic, 2)) + ' %\n')
    out3_1.write('rather methylated CGs in non genic regions: ' + str(round(CG_100_rather_met_non_genic, 2)) + ' %\n')
    out3_1.write('highly methylated CGs in non genic regions: ' + str(round(CG_100_met_non_genic, 2)) + ' %\n')
    out3_1.write('\n')

    out3_1.write('------------------------------------------------------------------------------\n')
    out3_1.write('\n')
 
    out3_1.write('-> CHG\n')
    out3_1.write('\n')
    out3_1.write('CHGs: ' + str(round(CHG_fraction, 2)) + ' %\n')
    out3_1.write('\n')
    out3_1.write('highly unmethylated CHGs: ' + str(round(CHG_100_unmet, 2)) + ' %\n')
    out3_1.write('rather unmethylated CHGs: ' + str(round(CHG_100_rather_unmet, 2)) + ' %\n')
    out3_1.write('50/50-methylated CHGs: ' + str(round(CHG_100_50met, 2)) + ' %\n')
    out3_1.write('rather methylated CHGs: ' + str(round(CHG_100_rather_met, 2)) + ' %\n')
    out3_1.write('highly methylated CHGs: ' + str(round(CHG_100_met, 2)) + ' %\n')
    out3_1.write('\n')
    out3_1.write('CHGs in CDSs: ' + str(round(CHG_100_CDS, 2)) + ' %\n')
    out3_1.write('highly unmethylated CHGs in CDSs: ' + str(round(CHG_100_unmet_CDS, 2)) + ' %\n')
    out3_1.write('rather unmethylated CHGs in CDSs: ' + str(round(CHG_100_rather_unmet_CDS, 2)) + ' %\n')
    out3_1.write('50/50-methylated CHGs in CDSs: ' + str(round(CHG_100_50met_CDS, 2)) + ' %\n')
    out3_1.write('rather methylated CHGs in CDSs: ' + str(round(CHG_100_rather_met_CDS, 2)) + ' %\n')
    out3_1.write('highly methylated CHGs in CDSs: ' + str(round(CHG_100_met_CDS, 2)) + ' %\n')
    out3_1.write('\n')
    out3_1.write('CHGs in UTRs: ' + str(round(CHG_100_UTR, 2)) + ' %\n')
    out3_1.write('highly unmethylated CHGs in UTRs: ' + str(round(CHG_100_unmet_UTR, 2)) + ' %\n')
    out3_1.write('rather unmethylated CHGs in UTRs: ' + str(round(CHG_100_rather_unmet_UTR, 2)) + ' %\n')
    out3_1.write('50/50-methylated CHGs in UTRs: ' + str(round(CHG_100_50met_UTR, 2)) + ' %\n')
    out3_1.write('rather methylated CHGs in UTRs: ' + str(round(CHG_100_rather_met_UTR, 2)) + ' %\n')
    out3_1.write('highly methylated CHGs in UTRs: ' + str(round(CHG_100_met_UTR, 2)) + ' %\n')
    out3_1.write('\n')
    out3_1.write('CHGs in exons: ' + str(round(CHG_100_exon, 2)) + ' %\n')
    out3_1.write('highly unmethylated CHGs in exons: ' + str(round(CHG_100_unmet_exon, 2)) + ' %\n')
    out3_1.write('rather unmethylated CHGs in exons: ' + str(round(CHG_100_rather_unmet_exon, 2)) + ' %\n')
    out3_1.write('50/50-methylated CHGs in exons: ' + str(round(CHG_100_50met_exon, 2)) + ' %\n')
    out3_1.write('rather methylated CHGs in exons: ' + str(round(CHG_100_rather_met_exon, 2)) + ' %\n')
    out3_1.write('highly methylated CHGs in exons: ' + str(round(CHG_100_met_exon, 2)) + ' %\n')
    out3_1.write('\n')
    out3_1.write('CHGs in introns: ' + str(round(CHG_100_intron, 2)) + ' %\n')
    out3_1.write('highly unmethylated CHGs in introns: ' + str(round(CHG_100_unmet_intron, 2)) + ' %\n')
    out3_1.write('rather unmethylated CHGs in introns: ' + str(round(CHG_100_rather_unmet_intron, 2)) + ' %\n')
    out3_1.write('50/50-methylated CHGs in introns: ' + str(round(CHG_100_50met_intron, 2)) + ' %\n')
    out3_1.write('rather methylated CHGs in introns: ' + str(round(CHG_100_rather_met_intron, 2)) + ' %\n')
    out3_1.write('highly methylated CHGs in introns: ' + str(round(CHG_100_met_intron, 2)) + ' %\n')
    out3_1.write('\n')
    out3_1.write('CHGs in genic regions: ' + str(round(CHG_100_genic, 2)) + ' %\n')
    out3_1.write('highly unmethylated CHGs in genic regions: ' + str(round(CHG_100_unmet_genic, 2)) + ' %\n')
    out3_1.write('rather unmethylated CHGs in genic regions: ' + str(round(CHG_100_rather_unmet_genic, 2)) + ' %\n')
    out3_1.write('50/50-methylated CHGs in genic regions: ' + str(round(CHG_100_50met_genic, 2)) + ' %\n')
    out3_1.write('rather methylated CHGs in genic regions: ' + str(round(CHG_100_rather_met_genic, 2)) + ' %\n')
    out3_1.write('highly methylated CHGs in genic regions: ' + str(round(CHG_100_met_genic, 2)) + ' %\n')
    out3_1.write('\n')
    out3_1.write('CHGs in 1kb upstream: ' + str(round(CHG_100_upstream, 2)) + ' %\n')
    out3_1.write('highly unmethylated CHGs in 1kb upstream: ' + str(round(CHG_100_unmet_upstream, 2)) + ' %\n')
    out3_1.write('rather unmethylated CHGs in 1kb upstream: ' + str(round(CHG_100_rather_unmet_upstream, 2)) + ' %\n')
    out3_1.write('50/50-methylated CHGs in 1kb upstream: ' + str(round(CHG_100_50met_upstream, 2)) + ' %\n')
    out3_1.write('rather methylated CHGs in 1kb upstream: ' + str(round(CHG_100_rather_met_upstream, 2)) + ' %\n')
    out3_1.write('highly methylated CHGs in 1kb upstream: ' + str(round(CHG_100_met_upstream, 2)) + ' %\n')
    out3_1.write('\n')
    out3_1.write('CHGs in 1kb downstream: ' + str(round(CHG_100_downstream, 2)) + ' %\n')
    out3_1.write('highly unmethylated CHGs in 1kb downstream: ' + str(round(CHG_100_unmet_downstream, 2)) + ' %\n')
    out3_1.write('rather unmethylated CHGs in 1kb downstream: ' + str(round(CHG_100_rather_unmet_downstream, 2)) + ' %\n')
    out3_1.write('50/50-methylated CHGs in 1kb downstream: ' + str(round(CHG_100_50met_downstream, 2)) + ' %\n')
    out3_1.write('rather methylated CHGs in 1kb downstream: ' + str(round(CHG_100_rather_met_downstream, 2)) + ' %\n')
    out3_1.write('highly methylated CHGs in 1kb downstream: ' + str(round(CHG_100_met_downstream, 2)) + ' %\n')
    out3_1.write('\n')
    out3_1.write('CHGs in non genic regions: ' + str(round(CHG_100_non_genic, 2)) + ' %\n')
    out3_1.write('highly unmethylated CHGs in non genic regions: ' + str(round(CHG_100_unmet_non_genic, 2)) + ' %\n')
    out3_1.write('rather unmethylated CHGs in non genic regions: ' + str(round(CHG_100_rather_unmet_non_genic, 2)) + ' %\n')
    out3_1.write('50/50-methylated CHGs in non genic regions: ' + str(round(CHG_100_50met_non_genic, 2)) + ' %\n')
    out3_1.write('rather methylated CHGs in non genic regions: ' + str(round(CHG_100_rather_met_non_genic, 2)) + ' %\n')
    out3_1.write('highly methylated CHGs in non genic regions: ' + str(round(CHG_100_met_non_genic, 2)) + ' %\n')
    out3_1.write('\n')

    out3_1.write('------------------------------------------------------------------------------\n')
    out3_1.write('\n')

    out3_1.write('-> CHH\n')
    out3_1.write('\n')
    out3_1.write('CHHs: ' + str(round(CHH_fraction, 2)) + ' %\n')
    out3_1.write('\n')
    out3_1.write('highly unmethylated CHHs: ' + str(round(CHH_100_unmet, 2)) + ' %\n')
    out3_1.write('rather unmethylated CHHs: ' + str(round(CHH_100_rather_unmet, 2)) + ' %\n')
    out3_1.write('50/50-methylated CHHs: ' + str(round(CHH_100_50met, 2)) + ' %\n')
    out3_1.write('rather methylated CHHs: ' + str(round(CHH_100_rather_met, 2)) + ' %\n')
    out3_1.write('highly methylated CHHs: ' + str(round(CHH_100_met, 2)) + ' %\n')
    out3_1.write('\n')
    out3_1.write('CHHs in CDSs: ' + str(round(CHH_100_CDS, 2)) + ' %\n')
    out3_1.write('highly unmethylated CHHs in CDSs: ' + str(round(CHH_100_unmet_CDS, 2)) + ' %\n')
    out3_1.write('rather unmethylated CHHs in CDSs: ' + str(round(CHH_100_rather_unmet_CDS, 2)) + ' %\n')
    out3_1.write('50/50-methylated CHHs in CDSs: ' + str(round(CHH_100_50met_CDS, 2)) + ' %\n')
    out3_1.write('rather methylated CHHs in CDSs: ' + str(round(CHH_100_rather_met_CDS, 2)) + ' %\n')
    out3_1.write('highly methylated CHHs in CDSs: ' + str(round(CHH_100_met_CDS, 2)) + ' %\n')
    out3_1.write('\n')
    out3_1.write('CHHs in UTRs: ' + str(round(CHH_100_UTR, 2)) + ' %\n')
    out3_1.write('highly unmethylated CHHs in UTRs: ' + str(round(CHH_100_unmet_UTR, 2)) + ' %\n')
    out3_1.write('rather unmethylated CHHs in UTRs: ' + str(round(CHH_100_rather_unmet_UTR, 2)) + ' %\n')
    out3_1.write('50/50-methylated CHHs in UTRs: ' + str(round(CHH_100_50met_UTR, 2)) + ' %\n')
    out3_1.write('rather methylated CHHs in UTRs: ' + str(round(CHH_100_rather_met_UTR, 2)) + ' %\n')
    out3_1.write('highly methylated CHHs in UTRs: ' + str(round(CHH_100_met_UTR, 2)) + ' %\n')
    out3_1.write('\n')
    out3_1.write('CHHs in exons: ' + str(round(CHH_100_exon, 2)) + ' %\n')
    out3_1.write('highly unmethylated CHHs in exons: ' + str(round(CHH_100_unmet_exon, 2)) + ' %\n')
    out3_1.write('rather unmethylated CHHs in exons: ' + str(round(CHH_100_rather_unmet_exon, 2)) + ' %\n')
    out3_1.write('50/50-methylated CHHs in exons: ' + str(round(CHH_100_50met_exon, 2)) + ' %\n')
    out3_1.write('rather methylated CHHs in exons: ' + str(round(CHH_100_rather_met_exon, 2)) + ' %\n')
    out3_1.write('highly methylated CHHs in exons: ' + str(round(CHH_100_met_exon, 2)) + ' %\n')
    out3_1.write('\n')
    out3_1.write('CHHs in introns: ' + str(round(CHH_100_intron, 2)) + ' %\n')
    out3_1.write('highly unmethylated CHHs in introns: ' + str(round(CHH_100_unmet_intron, 2)) + ' %\n')
    out3_1.write('rather unmethylated CHHs in introns: ' + str(round(CHH_100_rather_unmet_intron, 2)) + ' %\n')
    out3_1.write('50/50-methylated CHHs in introns: ' + str(round(CHH_100_50met_intron, 2)) + ' %\n')
    out3_1.write('rather methylated CHHs in introns: ' + str(round(CHH_100_rather_met_intron, 2)) + ' %\n')
    out3_1.write('highly methylated CHHs in introns: ' + str(round(CHH_100_met_intron, 2)) + ' %\n')
    out3_1.write('\n')
    out3_1.write('CHHs in genic regions: ' + str(round(CHH_100_genic, 2)) + ' %\n')
    out3_1.write('highly unmethylated CHHs in genic regions: ' + str(round(CHH_100_unmet_genic, 2)) + ' %\n')
    out3_1.write('rather unmethylated CHHs in genic regions: ' + str(round(CHH_100_rather_unmet_genic, 2)) + ' %\n')
    out3_1.write('50/50-methylated CHHs in genic regions: ' + str(round(CHH_100_50met_genic, 2)) + ' %\n')
    out3_1.write('rather methylated CHHs in genic regions: ' + str(round(CHH_100_rather_met_genic, 2)) + ' %\n')
    out3_1.write('highly methylated CHHs in genic regions: ' + str(round(CHH_100_met_genic, 2)) + ' %\n')
    out3_1.write('\n')
    out3_1.write('CHHs in 1kb upstream: ' + str(round(CHH_100_upstream, 2)) + ' %\n')
    out3_1.write('highly unmethylated CHHs in 1kb upstream: ' + str(round(CHH_100_unmet_upstream, 2)) + ' %\n')
    out3_1.write('rather unmethylated CHHs in 1kb upstream: ' + str(round(CHH_100_rather_unmet_upstream, 2)) + ' %\n')
    out3_1.write('50/50-methylated CHHs in 1kb upstream: ' + str(round(CHH_100_50met_upstream, 2)) + ' %\n')
    out3_1.write('rather methylated CHHs in 1kb upstream: ' + str(round(CHH_100_rather_met_upstream, 2)) + ' %\n')
    out3_1.write('highly methylated CHHs in 1kb upstream: ' + str(round(CHH_100_met_upstream, 2)) + ' %\n')
    out3_1.write('\n')
    out3_1.write('CHHs in 1kb downstream: ' + str(round(CHH_100_downstream, 2)) + ' %\n')
    out3_1.write('highly unmethylated CHHs in 1kb downstream: ' + str(round(CHH_100_unmet_downstream, 2)) + ' %\n')
    out3_1.write('rather unmethylated CHHs in 1kb downstream: ' + str(round(CHH_100_rather_unmet_downstream, 2)) + ' %\n')
    out3_1.write('50/50-methylated CHHs in 1kb downstream: ' + str(round(CHH_100_50met_downstream, 2)) + ' %\n')
    out3_1.write('rather methylated CHHs in 1kb downstream: ' + str(round(CHH_100_rather_met_downstream, 2)) + ' %\n')
    out3_1.write('highly methylated CHHs in 1kb downstream: ' + str(round(CHH_100_met_downstream, 2)) + ' %\n')
    out3_1.write('\n')
    out3_1.write('CHHs in non genic regions: ' + str(round(CHH_100_non_genic, 2)) + ' %\n')
    out3_1.write('highly unmethylated CHHs in non genic regions: ' + str(round(CHH_100_unmet_non_genic, 2)) + ' %\n')
    out3_1.write('rather unmethylated CHHs in non genic regions: ' + str(round(CHH_100_rather_unmet_non_genic, 2)) + ' %\n')
    out3_1.write('50/50-methylated CHHs in non genic regions: ' + str(round(CHH_100_50met_non_genic, 2)) + ' %\n')
    out3_1.write('rather methylated CHHs in non genic regions: ' + str(round(CHH_100_rather_met_non_genic, 2)) + ' %\n')
    out3_1.write('highly methylated CHHs in non genic regions: ' + str(round(CHH_100_met_non_genic, 2)) + ' %\n')
    out3_1.write('\n')
    
with open(fileLOG_3_2, "w") as out3_2:

    out3_2.write('[number_highly_unmethylated, %_highly_unmethylated, ..., number_highly_methylated, %_highly_methylated] in CDS, UTRs, introns, 1kb upstream, 1kb downstream and non_genic regions\n')
    for entry in all_5_mers:
        # out3_2.write(str(entry) + '\n')

        # Output-Alternative -> besser lesbar
        out3_2.write(str(entry[0]) + '\n')
        out3_2.write(str(entry[1:11]) + '\n')
        out3_2.write(str(entry[11:21]) + '\n')
        out3_2.write(str(entry[21:31]) + '\n')
        out3_2.write(str(entry[41:51]) + '\n')
        out3_2.write(str(entry[51:]) + '\n')
        out3_2.write(str(entry[31:41]) + '\n')
        out3_2.write('\n')

with open(fileLOG_3_3, "w") as out3_3:

    out3_3.write('[number_highly_unmethylated, %_highly_unmethylated, ..., number_highly_methylated, %_highly_methylated] in CG, CHG and CHH context\n')
    out3_3.write('\n')
    for entry in chr_tigs_distribution:
        # out3_3.write(str(entry) + '\n')
        
        # Output-Alternative -> besser lesbar
        out3_3.write(str(entry[0]) + '\n')
        out3_3.write(str(entry[1:11]) + '\n')
        out3_3.write(str(entry[11:21]) + '\n')
        out3_3.write(str(entry[21:]) + '\n')
        out3_3.write('\n')

    # print('all values were printed to result file sucessfully')

print('-> outputs were generated')
print('-> script 3 STARTS')
print('Time needed: ' + str(datetime.now() - startTime))

############################################################################################################################



#############################################################################################################################