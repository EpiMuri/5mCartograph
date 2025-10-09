# 5mCartograph_REPEATS_smk_3_summarize_calculate_and_output.py -> summarizes data of subsets and calculates data of interest

import re
from datetime import datetime
import numpy as np
import pandas as pd
import statistics
import csv
import sys

# start time
startTime = datetime.now()
print('-> script 3 STARTS')

############################################################################################################################
############################################################################################################################

# inputs
fileCounts5mer_chr1 = sys.argv[1]
fileCounts5mer_chr2 = sys.argv[2]
fileCounts5mer_chr3 = sys.argv[3]
fileCounts5mer_chr4 = sys.argv[4]
fileCounts5mer_chr5 = sys.argv[5]
fileCounts5mer_chr6 = sys.argv[6]
fileCounts5mer_chr7 = sys.argv[7]
fileCounts5mer_chr8 = sys.argv[8]
fileCounts5mer_chr9 = sys.argv[9]
fileCounts5mer_tigs = sys.argv[10]
fileCounts5mer_list = [fileCounts5mer_chr1, fileCounts5mer_chr2, fileCounts5mer_chr3, fileCounts5mer_chr4, fileCounts5mer_chr5, fileCounts5mer_chr6, fileCounts5mer_chr7, fileCounts5mer_chr8, fileCounts5mer_chr9, fileCounts5mer_tigs]
fileCountsChr_chr1 = sys.argv[11]
fileCountsChr_chr2 = sys.argv[12]
fileCountsChr_chr3 = sys.argv[13]
fileCountsChr_chr4 = sys.argv[14]
fileCountsChr_chr5 = sys.argv[15]
fileCountsChr_chr6 = sys.argv[16]
fileCountsChr_chr7 = sys.argv[17]
fileCountsChr_chr8 = sys.argv[18]
fileCountsChr_chr9 = sys.argv[19]
fileCountsChr_tigs = sys.argv[20]
fileCountsChr_list = [fileCountsChr_chr1, fileCountsChr_chr2, fileCountsChr_chr3, fileCountsChr_chr4, fileCountsChr_chr5, fileCountsChr_chr6, fileCountsChr_chr7, fileCountsChr_chr8, fileCountsChr_chr9, fileCountsChr_tigs]

# create location string for output 
file_out = sys.argv[1].rsplit('/', 2)[0]

# file path log output (also includes results)
fileLOG_3 = file_out + '/5mCartograph_REPEATS_out3/5mCartograph_REPEATS_smk_out3_log_and_results.txt'
print(fileLOG_3)

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
CHH_context = ['CAA', 'CAC', 'CAT', 'CCA', 'CCC', 'CCT', 'CTA', 'CTC', 'CTT']
CHG_context =  ['CAG', 'CCG', 'CTG']

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
    CG_unmet_repeat = 0
    CG_unmet_upstream = 0
    CG_unmet_downstream = 0
    CG_unmet_non_repeat = 0
    CG_rather_unmet_repeat = 0
    CG_rather_unmet_upstream = 0
    CG_rather_unmet_downstream = 0
    CG_rather_unmet_non_repeat = 0
    CG_50met_repeat = 0
    CG_50met_upstream = 0
    CG_50met_downstream = 0
    CG_50met_non_repeat = 0
    CG_rather_met_repeat = 0
    CG_rather_met_upstream = 0
    CG_rather_met_downstream = 0
    CG_rather_met_non_repeat = 0
    CG_met_repeat = 0
    CG_met_upstream = 0
    CG_met_downstream = 0
    CG_met_non_repeat = 0
    CHH_unmet_repeat = 0
    CHH_unmet_upstream = 0
    CHH_unmet_downstream = 0
    CHH_unmet_non_repeat = 0
    CHH_rather_unmet_repeat = 0
    CHH_rather_unmet_upstream = 0
    CHH_rather_unmet_downstream = 0
    CHH_rather_unmet_non_repeat = 0
    CHH_50met_repeat = 0
    CHH_50met_upstream = 0
    CHH_50met_downstream = 0
    CHH_50met_non_repeat = 0
    CHH_rather_met_repeat = 0
    CHH_rather_met_upstream = 0
    CHH_rather_met_downstream = 0
    CHH_rather_met_non_repeat = 0
    CHH_met_repeat = 0
    CHH_met_upstream = 0
    CHH_met_downstream = 0
    CHH_met_non_repeat = 0
    CHG_unmet_repeat = 0
    CHG_unmet_upstream = 0
    CHG_unmet_downstream = 0
    CHG_unmet_non_repeat = 0
    CHG_rather_unmet_repeat = 0
    CHG_rather_unmet_upstream = 0
    CHG_rather_unmet_downstream = 0
    CHG_rather_unmet_non_repeat = 0
    CHG_50met_repeat = 0
    CHG_50met_upstream = 0
    CHG_50met_downstream = 0
    CHG_50met_non_repeat = 0
    CHG_rather_met_repeat = 0
    CHG_rather_met_upstream = 0
    CHG_rather_met_downstream = 0
    CHG_rather_met_non_repeat = 0
    CHG_met_repeat = 0
    CHG_met_upstream = 0
    CHG_met_downstream = 0
    CHG_met_non_repeat = 0

    # (1) + (2)
    CG_repeat = 0
    CHH_repeat = 0
    CHG_repeat = 0
    CG_upstream = 0
    CHH_upstream = 0
    CHG_upstream = 0
    CG_downstream = 0
    CHH_downstream = 0
    CHG_downstream = 0
    CG_non_repeat = 0
    CHH_non_repeat = 0
    CHG_non_repeat = 0

    # (1) + (3)
    CG_unmet = 0
    CHH_unmet = 0
    CHG_unmet = 0
    CG_rather_unmet = 0
    CHH_rather_unmet = 0
    CHG_rather_unmet = 0
    CG_50met = 0
    CHH_50met = 0
    CHG_50met = 0
    CG_rather_met = 0
    CHH_rather_met = 0
    CHG_rather_met = 0
    CG_met = 0
    CHH_met = 0
    CHG_met = 0

    # (2) + (3)
    c_unmet_repeat = 0
    c_unmet_upstream = 0
    c_unmet_downstream = 0
    c_unmet_non_repeat = 0
    c_rather_unmet_repeat = 0
    c_rather_unmet_upstream = 0
    c_rather_unmet_downstream = 0
    c_rather_unmet_non_repeat = 0
    c_50met_repeat = 0
    c_50met_upstream = 0
    c_50met_downstream = 0
    c_50met_non_repeat = 0
    c_rather_met_repeat = 0
    c_rather_met_upstream = 0
    c_rather_met_downstream = 0
    c_rather_met_non_repeat = 0
    c_met_repeat = 0
    c_met_upstream = 0
    c_met_downstream = 0
    c_met_non_repeat = 0

    # (1)
    CG_fraction = 0
    CHH_fraction = 0
    CHG_fraction = 0

    # (2)
    c_repeat = 0
    c_upstream = 0
    c_downstream = 0
    c_non_repeat = 0

    # (3)
    c_unmet = 0
    c_rather_unmet = 0
    c_50met = 0
    c_rather_met = 0
    c_met = 0

    # variables for output of CG values in CG relation -> (1) + (2) + (3)
    CG_100_unmet_repeat = 0
    CG_100_unmet_upstream = 0
    CG_100_unmet_downstream = 0
    CG_100_unmet_non_repeat = 0
    CG_100_rather_unmet_repeat = 0
    CG_100_rather_unmet_upstream = 0
    CG_100_rather_unmet_downstream = 0
    CG_100_rather_unmet_non_repeat = 0
    CG_100_50met_repeat = 0
    CG_100_50met_upstream = 0
    CG_100_50met_downstream = 0
    CG_100_50met_non_repeat = 0
    CG_100_rather_met_repeat = 0
    CG_100_rather_met_upstream = 0
    CG_100_rather_met_downstream = 0
    CG_100_rather_met_non_repeat = 0
    CG_100_met_repeat = 0
    CG_100_met_upstream = 0
    CG_100_met_downstream = 0
    CG_100_met_non_repeat = 0
    CHH_100_unmet_repeat = 0
    CHH_100_unmet_upstream = 0
    CHH_100_unmet_downstream = 0
    CHH_100_unmet_non_repeat = 0
    CHH_100_rather_unmet_repeat = 0
    CHH_100_rather_unmet_upstream = 0
    CHH_100_rather_unmet_downstream = 0
    CHH_100_rather_unmet_non_repeat = 0
    CHH_100_50met_repeat = 0
    CHH_100_50met_upstream = 0
    CHH_100_50met_downstream = 0
    CHH_100_50met_non_repeat = 0
    CHH_100_rather_met_repeat = 0
    CHH_100_rather_met_upstream = 0
    CHH_100_rather_met_downstream = 0
    CHH_100_rather_met_non_repeat = 0
    CHH_100_met_repeat = 0
    CHH_100_met_upstream = 0
    CHH_100_met_downstream = 0
    CHH_100_met_non_repeat = 0
    CHG_100_unmet_repeat = 0
    CHG_100_unmet_upstream = 0
    CHG_100_unmet_downstream = 0
    CHG_100_unmet_non_repeat = 0
    CHG_100_rather_unmet_repeat = 0
    CHG_100_rather_unmet_upstream = 0
    CHG_100_rather_unmet_downstream = 0
    CHG_100_rather_unmet_non_repeat = 0
    CHG_100_50met_repeat = 0
    CHG_100_50met_upstream = 0
    CHG_100_50met_downstream = 0
    CHG_100_50met_non_repeat = 0
    CHG_100_rather_met_repeat = 0
    CHG_100_rather_met_upstream = 0
    CHG_100_rather_met_downstream = 0
    CHG_100_rather_met_non_repeat = 0
    CHG_100_met_repeat = 0
    CHG_100_met_upstream = 0
    CHG_100_met_downstream = 0
    CHG_100_met_non_repeat = 0

    # variables for output of CG values in CG relation -> (1) + (2)
    CG_100_repeat = 0
    CG_100_upstream = 0
    CG_100_downstream = 0
    CG_100_non_repeat = 0
    CHH_100_repeat = 0
    CHH_100_upstream = 0
    CHH_100_downstream = 0
    CHH_100_non_repeat = 0
    CHG_100_repeat = 0
    CHG_100_upstream = 0
    CHG_100_downstream = 0
    CHG_100_non_repeat = 0

    # variables for output of CG values in CG relation -> (1) + (3)
    CG_100_unmet = 0
    CG_100_rather_unmet = 0
    CG_100_50met = 0
    CG_100_rather_met = 0
    CG_100_met = 0
    CHH_100_unmet = 0
    CHH_100_rather_unmet = 0
    CHH_100_50met = 0
    CHH_100_rather_met = 0
    CHH_100_met = 0
    CHG_100_unmet = 0
    CHG_100_rather_unmet = 0
    CHG_100_50met = 0
    CHG_100_rather_met = 0
    CHG_100_met = 0

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
# CG_fraction...    c_repeat...            c_unmet...

# (1) + (2)             CG_repeat...
# (1) + (3)             CG_unmet...
# (2) + (3)             c_unmet_repeat...
# (1) + (2) + (3)       CG_unmet_repeat...

############################################################################################################################

# calculate realative frequencie & calculate values of interest -> context and exon, downstream ... and methylation status
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

    # print(k_mer)
    
    # calculate values of interest -> 5-mers and exons, downstreams in CG, CHH and CHG context
    if k_mer[0][2:4] in CG_context:
        # exon, downstream ... in CG -> (1) + (2) + (3)
        CG_unmet_repeat += k_mer[2]
        CG_unmet_upstream += k_mer[12]
        CG_unmet_downstream += k_mer[22]
        CG_unmet_non_repeat += k_mer[32]
        CG_rather_unmet_repeat += k_mer[4]
        CG_rather_unmet_upstream += k_mer[14]
        CG_rather_unmet_downstream += k_mer[24]
        CG_rather_unmet_non_repeat += k_mer[34]
        CG_50met_repeat += k_mer[6]
        CG_50met_upstream += k_mer[16]
        CG_50met_downstream += k_mer[26]
        CG_50met_non_repeat += k_mer[36]
        CG_rather_met_repeat += k_mer[8]
        CG_rather_met_upstream += k_mer[18]
        CG_rather_met_downstream += k_mer[28]
        CG_rather_met_non_repeat += k_mer[38]
        CG_met_repeat += k_mer[10]
        CG_met_upstream += k_mer[20]
        CG_met_downstream += k_mer[30]
        CG_met_non_repeat += k_mer[40]
    elif k_mer[0][2:] in CHH_context:
        # exon, downstream ... in CHH -> (1) + (2) + (3)
        CHH_unmet_repeat += k_mer[2]
        CHH_unmet_upstream += k_mer[12]
        CHH_unmet_downstream += k_mer[22]
        CHH_unmet_non_repeat += k_mer[32]
        CHH_rather_unmet_repeat += k_mer[4]
        CHH_rather_unmet_upstream += k_mer[14]
        CHH_rather_unmet_downstream += k_mer[24]
        CHH_rather_unmet_non_repeat += k_mer[34]
        CHH_50met_repeat += k_mer[6]
        CHH_50met_upstream += k_mer[16]
        CHH_50met_downstream += k_mer[26]
        CHH_50met_non_repeat += k_mer[36]
        CHH_rather_met_repeat += k_mer[8]
        CHH_rather_met_upstream += k_mer[18]
        CHH_rather_met_downstream += k_mer[28]
        CHH_rather_met_non_repeat += k_mer[38]
        CHH_met_repeat += k_mer[10]
        CHH_met_upstream += k_mer[20]
        CHH_met_downstream += k_mer[30]
        CHH_met_non_repeat += k_mer[40]
    elif k_mer[0][2:] in CHG_context:
        # exon, downstream ... in CHG -> (1) + (2) + (3)
        CHG_unmet_repeat += k_mer[2]
        CHG_unmet_upstream += k_mer[12]
        CHG_unmet_downstream += k_mer[22]
        CHG_unmet_non_repeat += k_mer[32]
        CHG_rather_unmet_repeat += k_mer[4]
        CHG_rather_unmet_upstream += k_mer[14]
        CHG_rather_unmet_downstream += k_mer[24]
        CHG_rather_unmet_non_repeat += k_mer[34]
        CHG_50met_repeat += k_mer[6]
        CHG_50met_upstream += k_mer[16]
        CHG_50met_downstream += k_mer[26]
        CHG_50met_non_repeat += k_mer[36]
        CHG_rather_met_repeat += k_mer[8]
        CHG_rather_met_upstream += k_mer[18]
        CHG_rather_met_downstream += k_mer[28]
        CHG_rather_met_non_repeat += k_mer[38]
        CHG_met_repeat += k_mer[10]
        CHG_met_upstream += k_mer[20]
        CHG_met_downstream += k_mer[30]
        CHG_met_non_repeat += k_mer[40]

print('-> data arrays were filled with relative data')

# calculate proportion of ... in CG, CHH and CHG context -> (1) + (2)
CG_repeat = CG_unmet_repeat + CG_rather_unmet_repeat + CG_50met_repeat + CG_rather_met_repeat + CG_met_repeat
CHH_repeat = CHH_unmet_repeat + CHH_rather_unmet_repeat + CHH_50met_repeat + CHH_rather_met_repeat + CHH_met_repeat
CHG_repeat = CHG_unmet_repeat + CHG_rather_unmet_repeat + CHG_50met_repeat + CHG_rather_met_repeat + CHG_met_repeat
CG_upstream = CG_unmet_upstream + CG_rather_unmet_upstream + CG_50met_upstream + CG_rather_met_upstream + CG_met_upstream
CHH_upstream = CHH_unmet_upstream + CHH_rather_unmet_upstream + CHH_50met_upstream + CHH_rather_met_upstream + CHH_met_upstream
CHG_upstream = CHG_unmet_upstream + CHG_rather_unmet_upstream + CHG_50met_upstream + CHG_rather_met_upstream + CHG_met_upstream
CG_downstream = CG_unmet_downstream + CG_rather_unmet_downstream + CG_50met_downstream + CG_rather_met_downstream + CG_met_downstream
CHH_downstream = CHH_unmet_downstream + CHH_rather_unmet_downstream + CHH_50met_downstream + CHH_rather_met_downstream + CHH_met_downstream
CHG_downstream = CHG_unmet_downstream + CHG_rather_unmet_downstream + CHG_50met_downstream + CHG_rather_met_downstream + CHG_met_downstream
CG_non_repeat = CG_unmet_non_repeat + CG_rather_unmet_non_repeat + CG_50met_non_repeat + CG_rather_met_non_repeat + CG_met_non_repeat
CHH_non_repeat = CHH_unmet_non_repeat + CHH_rather_unmet_non_repeat + CHH_50met_non_repeat + CHH_rather_met_non_repeat + CHH_met_non_repeat
CHG_non_repeat = CHG_unmet_non_repeat + CHG_rather_unmet_non_repeat + CHG_50met_non_repeat + CHG_rather_met_non_repeat + CHG_met_non_repeat

# calculation of proportions of CG, ... in unmet, rather_unmet, ... -> (1) + (3)
CG_unmet = CG_unmet_repeat + CG_unmet_upstream + CG_unmet_downstream + CG_unmet_non_repeat
CHH_unmet = CHH_unmet_repeat + CHH_unmet_upstream + CHH_unmet_downstream + CHH_unmet_non_repeat
CHG_unmet = CHG_unmet_repeat + CHG_unmet_upstream + CHG_unmet_downstream + CHG_unmet_non_repeat
CG_rather_unmet = CG_rather_unmet_repeat + CG_rather_unmet_upstream + CG_rather_unmet_downstream + CG_rather_unmet_non_repeat
CHH_rather_unmet = CHH_rather_unmet_repeat + CHH_rather_unmet_upstream + CHH_rather_unmet_downstream + CHH_rather_unmet_non_repeat
CHG_rather_unmet = CHG_rather_unmet_repeat + CHG_rather_unmet_upstream + CHG_rather_unmet_downstream + CHG_rather_unmet_non_repeat
CG_50met = CG_50met_repeat + CG_50met_upstream + CG_50met_downstream + CG_50met_non_repeat
CHH_50met = CHH_50met_repeat + CHH_50met_upstream + CHH_50met_downstream + CHH_50met_non_repeat
CHG_50met = CHG_50met_repeat + CHG_50met_upstream + CHG_50met_downstream + CHG_50met_non_repeat
CG_rather_met = CG_rather_met_repeat + CG_rather_met_upstream + CG_rather_met_downstream + CG_rather_met_non_repeat
CHH_rather_met = CHH_rather_met_repeat + CHH_rather_met_upstream + CHH_rather_met_downstream + CHH_rather_met_non_repeat
CHG_rather_met = CHG_rather_met_repeat + CHG_rather_met_upstream + CHG_rather_met_downstream + CHG_rather_met_non_repeat
CG_met = CG_met_repeat + CG_met_upstream + CG_met_downstream + CG_met_non_repeat
CHH_met = CHH_met_repeat + CHH_met_upstream + CHH_met_downstream + CHH_met_non_repeat
CHG_met = CHG_met_repeat + CHG_met_upstream + CHG_met_downstream + CHG_met_non_repeat

# calculation of proportions of unmet, ... repeat, ... -> (2) + (3) 
c_unmet_repeat = CG_unmet_repeat + CHH_unmet_repeat + CHG_unmet_repeat
c_unmet_upstream = CG_unmet_upstream + CHH_unmet_upstream + CHG_unmet_upstream
c_unmet_downstream = CG_unmet_downstream + CHH_unmet_downstream + CHG_unmet_downstream
c_unmet_non_repeat = CG_unmet_non_repeat + CHH_unmet_non_repeat + CHG_unmet_non_repeat
c_rather_unmet_repeat = CG_rather_unmet_repeat + CHH_rather_unmet_repeat + CHG_rather_unmet_repeat
c_rather_unmet_upstream = CG_rather_unmet_upstream + CHH_rather_unmet_upstream + CHG_rather_unmet_upstream
c_rather_unmet_downstream = CG_rather_unmet_downstream + CHH_rather_unmet_downstream + CHG_rather_unmet_downstream
c_rather_unmet_non_repeat = CG_rather_unmet_non_repeat + CHH_rather_unmet_non_repeat + CHG_rather_unmet_non_repeat
c_50met_repeat = CG_50met_repeat + CHH_50met_repeat + CHG_50met_repeat
c_50met_upstream = CG_50met_upstream + CHH_50met_upstream + CHG_50met_upstream
c_50met_downstream = CG_50met_downstream + CHH_50met_downstream + CHG_50met_downstream
c_50met_non_repeat = CG_50met_non_repeat + CHH_50met_non_repeat + CHG_50met_non_repeat
c_rather_met_repeat = CG_rather_met_repeat + CHH_rather_met_repeat + CHG_rather_met_repeat
c_rather_met_upstream = CG_rather_met_upstream + CHH_rather_met_upstream + CHG_rather_met_upstream
c_rather_met_downstream = CG_rather_met_downstream + CHH_rather_met_downstream + CHG_rather_met_downstream
c_rather_met_non_repeat = CG_rather_met_non_repeat + CHH_rather_met_non_repeat + CHG_rather_met_non_repeat
c_met_repeat = CG_met_repeat + CHH_met_repeat + CHG_met_repeat
c_met_upstream = CG_met_upstream + CHH_met_upstream + CHG_met_upstream
c_met_downstream = CG_met_downstream + CHH_met_downstream + CHG_met_downstream
c_met_non_repeat = CG_met_non_repeat + CHH_met_non_repeat + CHG_met_non_repeat

# calculate proportion of CG, CHH and CHG context genome wide -> (1)
CG_fraction = CG_unmet + CG_rather_unmet + CG_50met + CG_rather_met + CG_met
CHH_fraction = CHH_unmet + CHH_rather_unmet + CHH_50met + CHH_rather_met + CHH_met
CHG_fraction = CHG_unmet + CHG_rather_unmet + CHG_50met + CHG_rather_met + CHG_met

# calculate proportion of repeat, ... genome wide  -> (2)
c_repeat = CG_repeat + CHH_repeat + CHG_repeat
c_upstream =  CG_upstream + CHH_upstream + CHG_upstream
c_downstream = CG_downstream + CHH_downstream + CHG_downstream
c_non_repeat = CG_non_repeat + CHH_non_repeat + CHG_non_repeat

# calculate proportion of unmet, ... genome wide -> (3)
c_unmet = CG_unmet + CHH_unmet + CHG_unmet
c_rather_unmet = CG_rather_unmet + CHH_rather_unmet + CHG_rather_unmet
c_50met = CG_50met + CHH_50met + CHG_50met
c_rather_met = CG_rather_met + CHH_rather_met + CHG_rather_met
c_met = CG_met + CHH_met + CHG_met

print('-> values of interest were calculated')

############################################################################################################################

# calculation for proportion of meth, ... in each context for output (in relation to context) -> (1) + (2) + (3) / (1) + (2) / (1) + (3)
if CG_fraction != 0:
    # (1) + (2) + (3)
    CG_100_unmet_repeat = CG_unmet_repeat / CG_fraction * 100
    CG_100_unmet_upstream = CG_unmet_upstream / CG_fraction * 100
    CG_100_unmet_downstream = CG_unmet_downstream  / CG_fraction * 100
    CG_100_unmet_non_repeat = CG_unmet_non_repeat / CG_fraction * 100
    CG_100_rather_unmet_repeat = CG_rather_unmet_repeat / CG_fraction * 100
    CG_100_rather_unmet_upstream = CG_rather_unmet_upstream / CG_fraction * 100
    CG_100_rather_unmet_downstream = CG_rather_unmet_downstream  / CG_fraction * 100
    CG_100_rather_unmet_non_repeat = CG_rather_unmet_non_repeat / CG_fraction * 100
    CG_100_50met_repeat = CG_50met_repeat / CG_fraction * 100
    CG_100_50met_upstream = CG_50met_upstream / CG_fraction * 100
    CG_100_50met_downstream = CG_50met_downstream  / CG_fraction * 100
    CG_100_50met_non_repeat = CG_50met_non_repeat / CG_fraction * 100
    CG_100_rather_met_repeat = CG_rather_met_repeat / CG_fraction * 100
    CG_100_rather_met_upstream = CG_rather_met_upstream / CG_fraction * 100
    CG_100_rather_met_downstream = CG_rather_met_downstream  / CG_fraction * 100
    CG_100_rather_met_non_repeat = CG_rather_met_non_repeat / CG_fraction * 100
    CG_100_met_repeat = CG_met_repeat / CG_fraction * 100
    CG_100_met_upstream = CG_met_upstream / CG_fraction * 100
    CG_100_met_downstream = CG_met_downstream  / CG_fraction * 100
    CG_100_met_non_repeat = CG_met_non_repeat / CG_fraction * 100
    # (1) + (2)
    CG_100_repeat = CG_repeat / CG_fraction * 100
    CG_100_upstream = CG_upstream / CG_fraction * 100
    CG_100_downstream = CG_downstream / CG_fraction * 100
    CG_100_non_repeat = CG_non_repeat / CG_fraction * 100
    # (1) + (3)
    CG_100_unmet = CG_unmet / CG_fraction * 100
    CG_100_rather_unmet = CG_rather_unmet / CG_fraction * 100
    CG_100_50met = CG_50met / CG_fraction * 100
    CG_100_rather_met = CG_rather_met / CG_fraction * 100
    CG_100_met = CG_met / CG_fraction * 100
if CHH_fraction != 0:
    # (1) + (2) + (3)
    CHH_100_unmet_repeat = CHH_unmet_repeat / CHH_fraction * 100
    CHH_100_unmet_upstream = CHH_unmet_upstream / CHH_fraction * 100
    CHH_100_unmet_downstream = CHH_unmet_downstream  / CHH_fraction * 100
    CHH_100_unmet_non_repeat = CHH_unmet_non_repeat / CHH_fraction * 100
    CHH_100_rather_unmet_repeat = CHH_rather_unmet_repeat / CHH_fraction * 100
    CHH_100_rather_unmet_upstream = CHH_rather_unmet_upstream / CHH_fraction * 100
    CHH_100_rather_unmet_downstream = CHH_rather_unmet_downstream  / CHH_fraction * 100
    CHH_100_rather_unmet_non_repeat = CHH_rather_unmet_non_repeat / CHH_fraction * 100
    CHH_100_50met_repeat = CHH_50met_repeat / CHH_fraction * 100
    CHH_100_50met_upstream = CHH_50met_upstream / CHH_fraction * 100
    CHH_100_50met_downstream = CHH_50met_downstream  / CHH_fraction * 100
    CHH_100_50met_non_repeat = CHH_50met_non_repeat / CHH_fraction * 100
    CHH_100_rather_met_repeat = CHH_rather_met_repeat / CHH_fraction * 100
    CHH_100_rather_met_upstream = CHH_rather_met_upstream / CHH_fraction * 100
    CHH_100_rather_met_downstream = CHH_rather_met_downstream  / CHH_fraction * 100
    CHH_100_rather_met_non_repeat = CHH_rather_met_non_repeat / CHH_fraction * 100
    CHH_100_met_repeat = CHH_met_repeat / CHH_fraction * 100
    CHH_100_met_upstream = CHH_met_upstream / CHH_fraction * 100
    CHH_100_met_downstream = CHH_met_downstream  / CHH_fraction * 100
    CHH_100_met_non_repeat = CHH_met_non_repeat / CHH_fraction * 100
    # (1) + (2)
    CHH_100_repeat = CHH_repeat / CHH_fraction * 100
    CHH_100_upstream = CHH_upstream / CHH_fraction * 100
    CHH_100_downstream = CHH_downstream / CHH_fraction * 100
    CHH_100_non_repeat = CHH_non_repeat / CHH_fraction * 100
    # (1) + (3)
    CHH_100_unmet = CHH_unmet / CHH_fraction * 100
    CHH_100_rather_unmet = CHH_rather_unmet / CHH_fraction * 100
    CHH_100_50met = CHH_50met / CHH_fraction * 100
    CHH_100_rather_met = CHH_rather_met / CHH_fraction * 100
    CHH_100_met = CHH_met / CHH_fraction * 100
if CHG_fraction != 0:
    # (1) + (2) + (3)
    CHG_100_unmet_repeat = CHG_unmet_repeat / CHG_fraction * 100
    CHG_100_unmet_upstream = CHG_unmet_upstream / CHG_fraction * 100
    CHG_100_unmet_downstream = CHG_unmet_downstream  / CHG_fraction * 100
    CHG_100_unmet_non_repeat = CHG_unmet_non_repeat / CHG_fraction * 100
    CHG_100_rather_unmet_repeat = CHG_rather_unmet_repeat / CHG_fraction * 100
    CHG_100_rather_unmet_upstream = CHG_rather_unmet_upstream / CHG_fraction * 100
    CHG_100_rather_unmet_downstream = CHG_rather_unmet_downstream  / CHG_fraction * 100
    CHG_100_rather_unmet_non_repeat = CHG_rather_unmet_non_repeat / CHG_fraction * 100
    CHG_100_50met_repeat = CHG_50met_repeat / CHG_fraction * 100
    CHG_100_50met_upstream = CHG_50met_upstream / CHG_fraction * 100
    CHG_100_50met_downstream = CHG_50met_downstream  / CHG_fraction * 100
    CHG_100_50met_non_repeat = CHG_50met_non_repeat / CHG_fraction * 100
    CHG_100_rather_met_repeat = CHG_rather_met_repeat / CHG_fraction * 100
    CHG_100_rather_met_upstream = CHG_rather_met_upstream / CHG_fraction * 100
    CHG_100_rather_met_downstream = CHG_rather_met_downstream  / CHG_fraction * 100
    CHG_100_rather_met_non_repeat = CHG_rather_met_non_repeat / CHG_fraction * 100
    CHG_100_met_repeat = CHG_met_repeat / CHG_fraction * 100
    CHG_100_met_upstream = CHG_met_upstream / CHG_fraction * 100
    CHG_100_met_downstream = CHG_met_downstream  / CHG_fraction * 100
    CHG_100_met_non_repeat = CHG_met_non_repeat / CHG_fraction * 100
    # (1) + (2)
    CHG_100_repeat = CHG_repeat / CHG_fraction * 100
    CHG_100_upstream = CHG_upstream / CHG_fraction * 100
    CHG_100_downstream = CHG_downstream / CHG_fraction * 100
    CHG_100_non_repeat = CHG_non_repeat / CHG_fraction * 100
    # (1) + (3)
    CHG_100_unmet = CHG_unmet / CHG_fraction * 100
    CHG_100_rather_unmet = CHG_rather_unmet / CHG_fraction * 100
    CHG_100_50met = CHG_50met / CHG_fraction * 100
    CHG_100_rather_met = CHG_rather_met / CHG_fraction * 100
    CHG_100_met = CHG_met / CHG_fraction * 100

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
with open(fileLOG_3, "w") as out3:
    # return values of interest
    out3.write('-> tool\n')
    out3.write(__file__ + '\n')
    out3.write('\n')

    out3.write('-> description\n')
    out3.write('in the third of three steps the data counts are merged and result values of interst are calculated\n')
    out3.write('\n')

    out3.write('-> 5mer specific count inputs\n')
    for counts in fileCounts5mer_list:
        out3.write(counts + '\n')
    out3.write('\n')

    out3.write('-> chromosome/contig specific count inputs\n')
    for counts in fileCountsChr_list:
        out3.write(counts + '\n')
    out3.write('\n')
    
    out3.write('-> log and output\n')
    out3.write(fileLOG_3 + '\n')
    out3.write('\n')

    out3.write('------------------------------------------------------------------------------\n')
    out3.write('\n')

    out3.write('-> note 15.12.2024\n')
    out3.write('non repeat: means all cytosines, that are not related to a repeat (not in or within a 1000 bp zone around a repeat) \n')
    out3.write('to calculate all cytosines that are not in a repeat add all upstream and downstream Cs to non repeat\n')
    out3.write('\n')

    out3.write('------------------------------------------------------------------------------\n')
    out3.write('\n')

    out3.write('-> analysed cytosines\n')
    out3.write('\n')
    out3.write('total number of cytosines: ' + str(c_counter) + '\n')
    out3.write('\n')
    out3.write('unmethylated cytosines: ' + str(round(c_unmet, 2)) + ' %\n')
    out3.write('rather unmethylated cytosines: ' + str(round(c_rather_unmet, 2)) + ' %\n')
    out3.write('50/50-methylated cytosines: ' + str(round(c_50met, 2)) + ' %\n')
    out3.write('rather methylated cytosines: ' + str(round(c_rather_met, 2)) + ' %\n')
    out3.write('methylated cytosines: ' + str(round(c_met, 2)) + ' %\n')
    out3.write('\n')

    out3.write('cytosines in repeats: ' + str(round(c_repeat, 2)) + ' %\n')
    out3.write('unmethylated cytosines in repeats: ' + str(round(c_unmet_repeat, 2)) + ' %\n')
    out3.write('rather unmethylated cytosines in repeats: ' + str(round(c_rather_unmet_repeat, 2)) + ' %\n')
    out3.write('50/50-methylated cytosines in repeats: ' + str(round(c_50met_repeat, 2)) + ' %\n')
    out3.write('rather methylated cytosines in repeats: ' + str(round(c_rather_met_repeat, 2)) + ' %\n')
    out3.write('methylated cytosines in repeats: ' + str(round(c_met_repeat, 2)) + ' %\n')
    out3.write('\n')
    out3.write('cytosines in upstream: ' + str(round(c_upstream, 2)) + ' %\n')
    out3.write('unmethylated cytosines in upstream: ' + str(round(c_unmet_upstream, 2)) + ' %\n')
    out3.write('rather unmethylated cytosines in upstreams: ' + str(round(c_rather_unmet_upstream, 2)) + ' %\n')
    out3.write('50/50-methylated cytosines in upstream: ' + str(round(c_50met_upstream, 2)) + ' %\n')
    out3.write('rather methylated cytosines in upstream: ' + str(round(c_rather_met_upstream, 2)) + ' %\n')
    out3.write('methylated cytosines in upstream: ' + str(round(c_met_upstream, 2)) + ' %\n')
    out3.write('\n')
    out3.write('cytosines in downstream: ' + str(round(c_downstream, 2)) + ' %\n')
    out3.write('unmethylated cytosines in downstream: ' + str(round(c_unmet_downstream, 2)) + ' %\n')
    out3.write('rather unmethylated cytosines in downstream: ' + str(round(c_rather_unmet_downstream, 2)) + ' %\n')
    out3.write('50/50-methylated cytosines in downstream: ' + str(round(c_50met_downstream, 2)) + ' %\n')
    out3.write('rather methylated cytosines in downstream: ' + str(round(c_rather_met_downstream, 2)) + ' %\n')
    out3.write('methylated cytosines in downstream: ' + str(round(c_met_downstream, 2)) + ' %\n')
    out3.write('\n')
    out3.write('cytosines in non repeat regions: ' + str(round(c_non_repeat, 2)) + ' %\n')
    out3.write('unmethylated cytosines in non repeat regions: ' + str(round(c_unmet_non_repeat, 2)) + ' %\n')
    out3.write('rather unmethylated cytosines in non repeat regions: ' + str(round(c_rather_unmet_non_repeat, 2)) + ' %\n')
    out3.write('50/50-methylated cytosines in non repeat regions: ' + str(round(c_50met_non_repeat, 2)) + ' %\n')
    out3.write('rather methylated cytosines in non repeat regions: ' + str(round(c_rather_met_non_repeat, 2)) + ' %\n')
    out3.write('methylated cytosines in non repeat regions: ' + str(round(c_met_non_repeat, 2)) + ' %\n')
    out3.write('\n')

    out3.write('------------------------------------------------------------------------------\n')
    out3.write('\n')
    
    out3.write('-> CG\n')
    out3.write('\n')
    out3.write('CGs: ' + str(round(CG_fraction, 2)) + ' %\n')
    out3.write('\n')
    out3.write('unmethylated CGs: ' + str(round(CG_100_unmet, 2)) + ' %\n')
    out3.write('rather unmethylated CGs: ' + str(round(CG_100_rather_unmet, 2)) + ' %\n')
    out3.write('50/50-methylated CGs: ' + str(round(CG_100_50met, 2)) + ' %\n')
    out3.write('rather methylated CGs: ' + str(round(CG_100_rather_met, 2)) + ' %\n')
    out3.write('methylated CGs: ' + str(round(CG_100_met, 2)) + ' %\n')
    out3.write('\n')
    out3.write('CGs in repeats: ' + str(round(CG_100_repeat, 2)) + ' %\n')
    out3.write('unmethylated CGs in repeats: ' + str(round(CG_100_unmet_repeat, 2)) + ' %\n')
    out3.write('rather unmethylated CGs in repeats: ' + str(round(CG_100_rather_unmet_repeat, 2)) + ' %\n')
    out3.write('50/50-methylated CGs in repeats: ' + str(round(CG_100_50met_repeat, 2)) + ' %\n')
    out3.write('rather methylated CGs in repeats: ' + str(round(CG_100_rather_met_repeat, 2)) + ' %\n')
    out3.write('methylated CGs in repeats: ' + str(round(CG_100_met_repeat, 2)) + ' %\n')
    out3.write('\n')
    out3.write('CGs in upstream: ' + str(round(CG_100_upstream, 2)) + ' %\n')
    out3.write('unmethylated CGs in upstream: ' + str(round(CG_100_unmet_upstream, 2)) + ' %\n')
    out3.write('rather unmethylated CGs in upstream: ' + str(round(CG_100_rather_unmet_upstream, 2)) + ' %\n')
    out3.write('50/50-methylated CGs in upstream: ' + str(round(CG_100_50met_upstream, 2)) + ' %\n')
    out3.write('rather methylated CGs in upstream: ' + str(round(CG_100_rather_met_upstream, 2)) + ' %\n')
    out3.write('methylated CGs in upstream: ' + str(round(CG_100_met_upstream, 2)) + ' %\n')
    out3.write('\n')
    out3.write('CGs in downstream: ' + str(round(CG_100_downstream, 2)) + ' %\n')
    out3.write('unmethylated CGs in downstream: ' + str(round(CG_100_unmet_downstream, 2)) + ' %\n')
    out3.write('rather unmethylated CGs in downstream: ' + str(round(CG_100_rather_unmet_downstream, 2)) + ' %\n')
    out3.write('50/50-methylated CGs in downstream: ' + str(round(CG_100_50met_downstream, 2)) + ' %\n')
    out3.write('rather methylated CGs in downstream: ' + str(round(CG_100_rather_met_downstream, 2)) + ' %\n')
    out3.write('methylated CGs in downstream: ' + str(round(CG_100_met_downstream, 2)) + ' %\n')
    out3.write('\n')
    out3.write('CGs in non repeat regions: ' + str(round(CG_100_non_repeat, 2)) + ' %\n')
    out3.write('unmethylated CGs in non repeat regions: ' + str(round(CG_100_unmet_non_repeat, 2)) + ' %\n')
    out3.write('rather unmethylated CGs in non repeat regions: ' + str(round(CG_100_rather_unmet_non_repeat, 2)) + ' %\n')
    out3.write('50/50-methylated CGs in non repeat regions: ' + str(round(CG_100_50met_non_repeat, 2)) + ' %\n')
    out3.write('rather methylated CGs in non repeat regions: ' + str(round(CG_100_rather_met_non_repeat, 2)) + ' %\n')
    out3.write('methylated CGs in non repeat regions: ' + str(round(CG_100_met_non_repeat, 2)) + ' %\n')
    out3.write('\n')

    out3.write('------------------------------------------------------------------------------\n')
    out3.write('\n')
 
    out3.write('-> CHH\n')
    out3.write('\n')
    out3.write('CHHs: ' + str(round(CHH_fraction, 2)) + ' %\n')
    out3.write('\n')
    out3.write('unmethylated CHHs: ' + str(round(CHH_100_unmet, 2)) + ' %\n')
    out3.write('rather unmethylated CHHs: ' + str(round(CHH_100_rather_unmet, 2)) + ' %\n')
    out3.write('50/50-methylated CHHs: ' + str(round(CHH_100_50met, 2)) + ' %\n')
    out3.write('rather methylated CHHs: ' + str(round(CHH_100_rather_met, 2)) + ' %\n')
    out3.write('methylated CHHs: ' + str(round(CHH_100_met, 2)) + ' %\n')
    out3.write('\n')
    out3.write('CHHs in repeats: ' + str(round(CHH_100_repeat, 2)) + ' %\n')
    out3.write('unmethylated CHHs in repeats: ' + str(round(CHH_100_unmet_repeat, 2)) + ' %\n')
    out3.write('rather unmethylated CHHs in repeats: ' + str(round(CHH_100_rather_unmet_repeat, 2)) + ' %\n')
    out3.write('50/50-methylated CHHs in repeats: ' + str(round(CHH_100_50met_repeat, 2)) + ' %\n')
    out3.write('rather methylated CHHs in repeats: ' + str(round(CHH_100_rather_met_repeat, 2)) + ' %\n')
    out3.write('methylated CHHs in repeats: ' + str(round(CHH_100_met_repeat, 2)) + ' %\n')
    out3.write('\n')
    out3.write('CHHs in upstream: ' + str(round(CHH_100_upstream, 2)) + ' %\n')
    out3.write('unmethylated CHHs in upstream: ' + str(round(CHH_100_unmet_upstream, 2)) + ' %\n')
    out3.write('rather unmethylated CHHs in upstream: ' + str(round(CHH_100_rather_unmet_upstream, 2)) + ' %\n')
    out3.write('50/50-methylated CHHs in upstream: ' + str(round(CHH_100_50met_upstream, 2)) + ' %\n')
    out3.write('rather methylated CHHs in upstream: ' + str(round(CHH_100_rather_met_upstream, 2)) + ' %\n')
    out3.write('methylated CHHs in upstream: ' + str(round(CHH_100_met_upstream, 2)) + ' %\n')
    out3.write('\n')
    out3.write('CHHs in downstream: ' + str(round(CHH_100_downstream, 2)) + ' %\n')
    out3.write('unmethylated CHHs in downstream: ' + str(round(CHH_100_unmet_downstream, 2)) + ' %\n')
    out3.write('rather unmethylated CHHs in downstream: ' + str(round(CHH_100_rather_unmet_downstream, 2)) + ' %\n')
    out3.write('50/50-methylated CHHs in downstream: ' + str(round(CHH_100_50met_downstream, 2)) + ' %\n')
    out3.write('rather methylated CHHs in downstream: ' + str(round(CHH_100_rather_met_downstream, 2)) + ' %\n')
    out3.write('methylated CHHs in downstream: ' + str(round(CHH_100_met_downstream, 2)) + ' %\n')
    out3.write('\n')
    out3.write('CHHs in non repeat regions: ' + str(round(CHH_100_non_repeat, 2)) + ' %\n')
    out3.write('unmethylated CHHs in non repeat regions: ' + str(round(CHH_100_unmet_non_repeat, 2)) + ' %\n')
    out3.write('rather unmethylated CHHs in non repeat regions: ' + str(round(CHH_100_rather_unmet_non_repeat, 2)) + ' %\n')
    out3.write('50/50-methylated CHHs in non repeat regions: ' + str(round(CHH_100_50met_non_repeat, 2)) + ' %\n')
    out3.write('rather methylated CHHs in non repeat regions: ' + str(round(CHH_100_rather_met_non_repeat, 2)) + ' %\n')
    out3.write('methylated CHHs in non repeat regions: ' + str(round(CHH_100_met_non_repeat, 2)) + ' %\n')
    out3.write('\n')

    out3.write('------------------------------------------------------------------------------\n')
    out3.write('\n')

    out3.write('-> CHG\n')
    out3.write('\n')
    out3.write('CHGs: ' + str(round(CHG_fraction, 2)) + ' %\n')
    out3.write('\n')
    out3.write('unmethylated CHGs: ' + str(round(CHG_100_unmet, 2)) + ' %\n')
    out3.write('rather unmethylated CHGs: ' + str(round(CHG_100_rather_unmet, 2)) + ' %\n')
    out3.write('50/50-methylated CHGs: ' + str(round(CHG_100_50met, 2)) + ' %\n')
    out3.write('rather methylated CHGs: ' + str(round(CHG_100_rather_met, 2)) + ' %\n')
    out3.write('methylated CHGs: ' + str(round(CHG_100_met, 2)) + ' %\n')
    out3.write('\n')
    out3.write('CHGs in repeats: ' + str(round(CHG_100_repeat, 2)) + ' %\n')
    out3.write('unmethylated CHGs in repeats: ' + str(round(CHG_100_unmet_repeat, 2)) + ' %\n')
    out3.write('rather unmethylated CHGs in repeats: ' + str(round(CHG_100_rather_unmet_repeat, 2)) + ' %\n')
    out3.write('50/50-methylated CHGs in repeats: ' + str(round(CHG_100_50met_repeat, 2)) + ' %\n')
    out3.write('rather methylated CHGs in repeats: ' + str(round(CHG_100_rather_met_repeat, 2)) + ' %\n')
    out3.write('methylated CHGs in repeats: ' + str(round(CHG_100_met_repeat, 2)) + ' %\n')
    out3.write('\n')
    out3.write('CHGs in upstream: ' + str(round(CHG_100_upstream, 2)) + ' %\n')
    out3.write('unmethylated CHGs in upstream: ' + str(round(CHG_100_unmet_upstream, 2)) + ' %\n')
    out3.write('rather unmethylated CHGs in upstream: ' + str(round(CHG_100_rather_unmet_upstream, 2)) + ' %\n')
    out3.write('50/50-methylated CHGs in upstream: ' + str(round(CHG_100_50met_upstream, 2)) + ' %\n')
    out3.write('rather methylated CHGs in upstream: ' + str(round(CHG_100_rather_met_upstream, 2)) + ' %\n')
    out3.write('methylated CHGs in upstream: ' + str(round(CHG_100_met_upstream, 2)) + ' %\n')
    out3.write('\n')
    out3.write('CHGs in downstream: ' + str(round(CHG_100_downstream, 2)) + ' %\n')
    out3.write('unmethylated CHGs in downstream: ' + str(round(CHG_100_unmet_downstream, 2)) + ' %\n')
    out3.write('rather unmethylated CHGs in downstream: ' + str(round(CHG_100_rather_unmet_downstream, 2)) + ' %\n')
    out3.write('50/50-methylated CHGs in downstream: ' + str(round(CHG_100_50met_downstream, 2)) + ' %\n')
    out3.write('rather methylated CHGs in downstream: ' + str(round(CHG_100_rather_met_downstream, 2)) + ' %\n')
    out3.write('methylated CHGs in downstream: ' + str(round(CHG_100_met_downstream, 2)) + ' %\n')
    out3.write('\n')
    out3.write('CHGs in non repeat regions: ' + str(round(CHG_100_non_repeat, 2)) + ' %\n')
    out3.write('unmethylated CHGs in non repeat regions: ' + str(round(CHG_100_unmet_non_repeat, 2)) + ' %\n')
    out3.write('rather unmethylated CHGs in non repeat regions: ' + str(round(CHG_100_rather_unmet_non_repeat, 2)) + ' %\n')
    out3.write('50/50-methylated CHGs in non repeat regions: ' + str(round(CHG_100_50met_non_repeat, 2)) + ' %\n')
    out3.write('rather methylated CHGs in non repeat regions: ' + str(round(CHG_100_rather_met_non_repeat, 2)) + ' %\n')
    out3.write('methylated CHGs in non repeat regions: ' + str(round(CHG_100_met_non_repeat, 2)) + ' %\n')
    out3.write('\n')
    
    out3.write('------------------------------------------------------------------------------\n')
    out3.write('\n')

    out3.write('[number_unmethylated, %_unmethylated, ..., number_methylated, %_methylated] in repeat, 1 kb upstream, 1 kb downstream and non_repeat regions\n')
    for entry in all_5_mers:
        # out3.write(str(entry) + '\n')

        # Output-Alternative -> besser lesbar
        out3.write(str(entry[0]) + '\n')
        out3.write(str(entry[1:11]) + '\n')
        out3.write(str(entry[11:21]) + '\n')
        out3.write(str(entry[21:31]) + '\n')
        out3.write(str(entry[31:41]) + '\n')
        out3.write('\n')

    out3.write('------------------------------------------------------------------------------\n')
    out3.write('\n')

    for chr_tigs in chr_tigs_distribution:
        out3.write(str(chr_tigs[0]) + '\n')
        out3.write('unmethylated Cs on ' + str(chr_tigs[0]) + ' (of all Cs analysed): ' + str(round((chr_tigs[2] + chr_tigs[12] + chr_tigs[22]), 2)) + ' %\n')
        out3.write('| ' + str(round(chr_tigs[2], 2)) + ' % CG | ' + str(round(chr_tigs[12], 2)) + ' % CHH | ' + str(round(chr_tigs[22], 2)) + ' % CHG |\n')

        out3.write('rather unmethylated Cs on ' + str(chr_tigs[0]) + ' (of all Cs analysed): ' + str(round((chr_tigs[4] + chr_tigs[14] + chr_tigs[24]), 2)) + ' %\n')
        out3.write('| ' + str(round(chr_tigs[4], 2)) + ' % CG | ' + str(round(chr_tigs[14], 2)) + ' % CHH | ' + str(round(chr_tigs[24], 2)) + ' % CHG |\n')

        out3.write('50/50-methylated Cs on ' + str(chr_tigs[0]) + ' (of all Cs analysed): ' + str(round((chr_tigs[6] + chr_tigs[16] + chr_tigs[26]), 2)) + ' %\n')
        out3.write('| ' + str(round(chr_tigs[6], 2)) + ' % CG | ' + str(round(chr_tigs[16], 2)) + ' % CHH | ' + str(round(chr_tigs[26], 2)) + ' % CHG |\n')

        out3.write('rather methylated Cs on ' + str(chr_tigs[0]) + ' (of all Cs analysed): ' + str(round((chr_tigs[8] + chr_tigs[18] + chr_tigs[28]), 2)) + ' %\n')
        out3.write('| ' + str(round(chr_tigs[8], 2)) + ' % CG | ' + str(round(chr_tigs[18], 2)) + ' % CHH | ' + str(round(chr_tigs[28], 2)) + ' % CHG |\n')

        out3.write('methylated Cs on ' + str(chr_tigs[0]) + ' (of all Cs analysed): ' + str(round((chr_tigs[10] + chr_tigs[20] + chr_tigs[30]), 2)) + ' %\n')
        out3.write('| ' + str(round(chr_tigs[10], 2)) + ' % CG | ' + str(round(chr_tigs[20], 2)) + ' % CHH | ' + str(round(chr_tigs[30], 2)) + ' % CHG |\n')
        out3.write('\n')
   
    out3.write('------------------------------------------------------------------------------\n')
    out3.write('\n')

    out3.write('[number_unmethylated, %_unmethylated, ..., number_methylated, %_methylated] in CG, CHH and CHG context\n')
    out3.write('\n')
    for entry in chr_tigs_distribution:
        # out3.write(str(entry) + '\n')
        
        # Output-Alternative -> besser lesbar
        out3.write(str(entry[0]) + '\n')
        out3.write(str(entry[1:11]) + '\n')
        out3.write(str(entry[11:21]) + '\n')
        out3.write(str(entry[21:]) + '\n')
        out3.write('\n')

    out3.write('\n')

    # print('all values were printed to result file sucessfully')

print('-> outputs were generated')
print('-> script 3 STARTS')
print('Time needed: ' + str(datetime.now() - startTime))

############################################################################################################################



#############################################################################################################################