# 5mCartograph_REPEATS_smk_2_chr_specific_counts.py -> calculate 5mer and chromosome specific data (the chromosomes are processed individually, while the contigs run in the 10th slot)


import re
from datetime import datetime
import pandas as pd
import statistics
import sys
import os

# start time
startTime = datetime.now()
print('-> script 2 STARTS')

############################################################################################################################
############################################################################################################################

# python3 [script used] [partial DSP result file path] [partial repeat_annotation file path] [folder_name] [subfolder_name]

# get chromosome for file names
chrom = re.split('_', sys.argv[1])[-1].strip()
chrom = re.split('\.', chrom)[0]

# DSP result file path for one chromosome / all tigs
fileDSPout = sys.argv[1]

# repeat_annotation file path for one chromosome / all tigs
fileAnno = sys.argv[2]

# create location string for output 
file_out = sys.argv[3] + sys.argv[4]

# file path output
fileCounts5mer = file_out + '5mCartograph_REPEATS_out2_' + chrom + '/5mCartograph_REPEATS_out2_counts_per_5mer_' + chrom + '.txt'
fileCountsChr = file_out + '5mCartograph_REPEATS_out2_' + chrom + '/5mCartograph_REPEATS_out2_counts_per_chromosome_' + chrom + '.txt'

# file path log output
fileLOG_2 = file_out + '5mCartograph_REPEATS_out2_' + chrom + '/5mCartograph_REPEATS_out2_log_' + chrom + '.txt'

############################################################################################################################
############################################################################################################################

# create result folder
os.makedirs(file_out + '5mCartograph_REPEATS_out2_' + chrom + '/', exist_ok=True)

############################################################################################################################
############################################################################################################################

# variables

# collect all 5-mers and how often they are unmethylated and methylated
# [[5-mer, #_unmethylated, %_unmethylated, #_methylated, %_methylated], [...], ...]
all_5_mers = []
chr_tigs_distribution = []
unique_chr_tigss = []
current_5_mer_index = None
current_chr_tigs_index = None

# number of columns in DeepSignal Plant output 
columns = 11
# index of column with chr_tigs
chr_tigs_column = 0
# index of column with k-mer
k_mer_column = 10
# index of column with C position
C_position_column = 1
# index of called_label -> methylation probality
met_prob_column = 5
# index of number of reads
read_count_column = 8
# variable for average probability
current_average_met_probability = 0
# for storage of current k_mer
k_mer_sequence = ''
# for storage and creation of new k-mer entries
current_5_mer_list = None
# for storage and creation of new chr_tigs entries
current_chr_tigs_list = None

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

print('-> variables were created')

############################################################################################################################
############################################################################################################################

# write log head
with open(fileLOG_2, 'w') as outLog:
    outLog.write('-> tool\n')
    outLog.write(__file__ + '\n')
    outLog.write('\n')

    outLog.write('-> description\n')
    outLog.write('in the second of three steps DeepSignal-plant output is analysed per chromosome/contigs\n')
    outLog.write('the cytosines are counted in different categories including their methylation probability and location\n')
    outLog.write('\n')

    outLog.write('-> DeepSignal-plant output as input\n')
    outLog.write(fileDSPout + '\n')
    outLog.write('\n')

    outLog.write('-> repeat annotation input\n')
    outLog.write(fileAnno + '\n')
    outLog.write('\n')

    outLog.write('-> output\n')
    outLog.write(fileCounts5mer + '\n')
    outLog.write(fileCountsChr + '\n')
    outLog.write('\n')

    outLog.write('->  log output\n')
    outLog.write(fileLOG_2 + '\n')
    outLog.write('\n')

    outLog.write('->  progress: idx is printed each 1000000 lines \n')

############################################################################################################################
############################################################################################################################

# read in file 
with open(fileDSPout, "r") as inDSP:
    # string with elements that contain one row of text file each
    linesC = inDSP.readlines()
print('-> DeepSignal-plant output was read in -> ' + str(len(linesC)) + ' lines')

############################################################################################################################

# read in repeat_annotation file
with open(fileAnno, 'r') as inAnno:
    linesAnno = inAnno.readlines()
print('-> repeat annotation was read in')

############################################################################################################################
############################################################################################################################

# create giant list of lists of lists to manage the anno file

# giant list contains sublists, one sublist for each repeat
# structure giantAnnoList
# [ [...], [...], [...], ... ] 
giantAnnoList = []
# each sublist contains: chromosmome, start, end

# cycles through all lines of fileAnno to find repeats, then all attributes of the repeat and its subunits are collected and stored in giantAnnoList
# cycles through linesAnno
for idx,lineAnno in enumerate(linesAnno):

    # makes sure that data lines are analysed only
    if (lineAnno[0] != '#') and (len(re.split('\t', lineAnno)) >= 9):

        # finds all lines with repeat data
        if re.split('\t', lineAnno)[2].strip() == 'dispersed_repeat' or re.split('\t', lineAnno)[2].strip() == 'repeat':
            
            giantAnnoList.append([re.split('\t', lineAnno)[0].strip(), int(re.split('\t', lineAnno)[3].strip()), int(re.split('\t', lineAnno)[4].strip())])

print('-> giant list of repeats was created')

############################################################################################################################
############################################################################################################################

# create empty list of lists with all methylation contexts
for context in contexts:
    current_5_mer_list = [0] * 61
    current_5_mer_list[0] = context
    all_5_mers.append(current_5_mer_list)            

############################################################################################################################
############################################################################################################################

# calculation of methylation distribution on chromosomes/contigs
# calculation of methylation distribution of 5-mers
# the 5-mers are categoriesed in one ways:
# highly unmethylated, rather unmethylated, 50/50-methylated, rather methylated, highly methylated

# find all C-positions and corresponding k-mers -> go through all lines
for idx, lineC in enumerate(linesC):
    
    # as a progress controll in every 1000000th cycle the cycle number is printed
    if idx % 1000000 == 0:
        # print(idx)
        with open(fileLOG_2, "a") as out2:
            # final return
            out2.write(str(idx) + '\n')
                
    # if line is not empty, line is analysed
    if len(re.split('\t', lineC)) >= columns:

        # calcuate average methylation probability using re.split()
        current_average_met_probability = float(re.split('\t', lineC)[met_prob_column])/float(re.split('\t', lineC)[read_count_column])

        # chromosome/contig (chr_tigs) of the current lineC (C)
        current_chr_tigs = re.split('\t', lineC)[chr_tigs_column].strip()

        # position of the current lineC (C)
        current_C_position = int(re.split('\t', lineC)[C_position_column].strip())

        # k_mer_sequence of the current lineC (C)
        k_mer_sequence = re.split('\t', lineC)[k_mer_column].strip()

        ####################################################################################################################

        # analysis of methylation distribution of 5-mers
        # in different repeattic, non_genic regions

        # count data lines
        c_counter += 1
        
        # get index to of the corresponding context
        current_5_mer_index = contexts.index(k_mer_sequence)
        
        # analyse whether the C is in a repeat (+ upstream + downstream)
        # repeat -> repeat
        # upstream -> upstream
        # downstream -> downstream
        
        # loop to cycle through all repeats of repeat_annotation
        # if C in repeat -> repeat [1-10]
        # if C 1 kb upstream -> upstream [11-20]
        # if C 1 kb downstream -> downstream [21-30]
        # if C is not close to repeat -> non repeat [31-40]

        # to check whether C has already been sorted to a group
        is_in_repeat = False
        is_in_1kb_upstream = False
        is_in_1kb_downstream = False

        # cycles through all repeats in giantAnnoList
        for repeat in giantAnnoList:

            # check whether current C is in current repeat (chromosomes and borders)
            if (repeat[0] == current_chr_tigs) and (current_C_position >= repeat[1]) and (current_C_position <= repeat[2]):

                # is in a repeat!
                # repeat values are later on calculated by addition of CDS, UTR and intron values
                is_in_repeat = True
            
            # check whether current C is in 1 kb before current repeat (upstream)
            elif (repeat[0] == current_chr_tigs) and (current_C_position >= (repeat[1] - 1000)) and (current_C_position < repeat[1]):

                # is in 1 kb upstream!
                is_in_1kb_upstream = True

            # check whether current C is in 1 kb after current repeat (downstream)
            elif (repeat[0] == current_chr_tigs) and (current_C_position > repeat[2]) and (current_C_position <= (repeat[2] + 1000)):
                # is in 1 kb downstream!
                is_in_1kb_downstream = True

        # check whether C is in repeat
        if is_in_repeat == True:

            # discriminate between unmethylated and methylated Cs -> 5-mers
            if current_average_met_probability >= 0 and current_average_met_probability <= 0.2:
                # counts all unmethalted Cs
                all_5_mers[current_5_mer_index][1] += 1
            elif current_average_met_probability > 0.2 and current_average_met_probability <= 0.4:
                # counts all rather unmethylated Cs
                all_5_mers[current_5_mer_index][3] += 1
            elif current_average_met_probability > 0.4 and current_average_met_probability <= 0.6:
                # counts all 50/50-methylated Cs
                all_5_mers[current_5_mer_index][5] += 1
            elif current_average_met_probability > 0.6 and current_average_met_probability <= 0.8:
                # counts all rather methylated Cs
                all_5_mers[current_5_mer_index][7] += 1
            elif current_average_met_probability > 0.8 and current_average_met_probability <= 1:
                # counts all methylated Cs
                all_5_mers[current_5_mer_index][9] += 1

        # check whether C is in 1 kb upstream of repeat
        elif is_in_1kb_upstream == True:

            # discriminate between unmethylated and methylated Cs -> 5-mers
            if current_average_met_probability >= 0 and current_average_met_probability <= 0.2:
                # counts all unmethalted Cs
                all_5_mers[current_5_mer_index][11] += 1
            elif current_average_met_probability > 0.2 and current_average_met_probability <= 0.4:
                # counts all rather unmethylated Cs
                all_5_mers[current_5_mer_index][13] += 1
            elif current_average_met_probability > 0.4 and current_average_met_probability <= 0.6:
                # counts all 50/50-methylated Cs
                all_5_mers[current_5_mer_index][15] += 1
            elif current_average_met_probability > 0.6 and current_average_met_probability <= 0.8:
                # counts all rather methylated Cs
                all_5_mers[current_5_mer_index][17] += 1
            elif current_average_met_probability > 0.8 and current_average_met_probability <= 1:
                # counts all methylated Cs
                all_5_mers[current_5_mer_index][19] += 1 

        # check whether C is in 1 kb downstream of repeat
        elif is_in_1kb_downstream == True:

            # discriminate between unmethylated and methylated Cs -> 5-mers
            if current_average_met_probability >= 0 and current_average_met_probability <= 0.2:
                # counts all unmethalted Cs
                all_5_mers[current_5_mer_index][21] += 1
            elif current_average_met_probability > 0.2 and current_average_met_probability <= 0.4:
                # counts all rather unmethylated Cs
                all_5_mers[current_5_mer_index][23] += 1
            elif current_average_met_probability > 0.4 and current_average_met_probability <= 0.6:
                # counts all 50/50-methylated Cs
                all_5_mers[current_5_mer_index][25] += 1
            elif current_average_met_probability > 0.6 and current_average_met_probability <= 0.8:
                # counts all rather methylated Cs
                all_5_mers[current_5_mer_index][27] += 1
            elif current_average_met_probability > 0.8 and current_average_met_probability <= 1:
                # counts all methylated Cs
                all_5_mers[current_5_mer_index][29] += 1 

        # C is not in repeat or close to it
        else:
            
            # C is not repeat related

            # discriminate between unmethylated and methylated Cs -> 5-mers
            if current_average_met_probability >= 0 and current_average_met_probability <= 0.2:
                # counts all unmethalted Cs
                all_5_mers[current_5_mer_index][31] += 1
            elif current_average_met_probability > 0.2 and current_average_met_probability <= 0.4:
                # counts all rather unmethylated Cs
                all_5_mers[current_5_mer_index][33] += 1
            elif current_average_met_probability > 0.4 and current_average_met_probability <= 0.6:
                # counts all 50/50-methylated Cs
                all_5_mers[current_5_mer_index][35] += 1
            elif current_average_met_probability > 0.6 and current_average_met_probability <= 0.8:
                # counts all rather methylated Cs
                all_5_mers[current_5_mer_index][37] += 1
            elif current_average_met_probability > 0.8 and current_average_met_probability <= 1:
                # counts all methylated Cs
                all_5_mers[current_5_mer_index][39] += 1

        ####################################################################################################################

        # analysis of methylation distribution on chromosomes/contigs

        # add unknown chr_tigs to unique_chr_tigss and chr_tigs_distribution
        # sublist structure: 
        # [0] : chromosome/contig
        # [1-10] : methylation distribution in CG context
        # [11-20] : methylation distribution in CHG context
        # [21-30] : methylation distribution in CHH context
        if current_chr_tigs not in unique_chr_tigss:
            unique_chr_tigss.append(current_chr_tigs)
            current_chr_tigs_list = [0] * 31
            current_chr_tigs_list[0] = current_chr_tigs
            chr_tigs_distribution.append(current_chr_tigs_list)
        current_chr_tigs_index = unique_chr_tigss.index(current_chr_tigs)

        if k_mer_sequence[2:4] in CG_context:
            # discriminate between unmethylated and methylated Cs -> chr_tigs
            if current_average_met_probability >= 0 and current_average_met_probability <= 0.2:
                # counts all unmethalted Cs
                chr_tigs_distribution[current_chr_tigs_index][1] += 1
            elif current_average_met_probability > 0.2 and current_average_met_probability <= 0.4:
                # counts all rather unmethylated Cs
                chr_tigs_distribution[current_chr_tigs_index][3] += 1
            elif current_average_met_probability > 0.4 and current_average_met_probability <= 0.6:
                # counts all 50/50-methylated Cs
                chr_tigs_distribution[current_chr_tigs_index][5] += 1
            elif current_average_met_probability > 0.6 and current_average_met_probability <= 0.8:
                # counts all rather methylated Cs
                chr_tigs_distribution[current_chr_tigs_index][7] += 1
            elif current_average_met_probability > 0.8 and current_average_met_probability <= 1:
                # counts all methylated Cs
                chr_tigs_distribution[current_chr_tigs_index][9] += 1
        elif k_mer_sequence[2:] in CHG_context:
            # discriminate between unmethylated and methylated Cs -> chr_tigs
            if current_average_met_probability >= 0 and current_average_met_probability <= 0.2:
                # counts all unmethalted Cs
                chr_tigs_distribution[current_chr_tigs_index][11] += 1
            elif current_average_met_probability > 0.2 and current_average_met_probability <= 0.4:
                # counts all rather unmethylated Cs
                chr_tigs_distribution[current_chr_tigs_index][13] += 1
            elif current_average_met_probability > 0.4 and current_average_met_probability <= 0.6:
                # counts all 50/50-methylated Cs
                chr_tigs_distribution[current_chr_tigs_index][15] += 1
            elif current_average_met_probability > 0.6 and current_average_met_probability <= 0.8:
                # counts all rather methylated Cs
                chr_tigs_distribution[current_chr_tigs_index][17] += 1
            elif current_average_met_probability > 0.8 and current_average_met_probability <= 1:
                # counts all methylated Cs
                chr_tigs_distribution[current_chr_tigs_index][19] += 1
        elif k_mer_sequence[2:] in CHH_context:
            # discriminate between unmethylated and methylated Cs -> chr_tigs
            if current_average_met_probability >= 0 and current_average_met_probability <= 0.2:
                # counts all unmethalted Cs
                chr_tigs_distribution[current_chr_tigs_index][21] += 1
            elif current_average_met_probability > 0.2 and current_average_met_probability <= 0.4:
                # counts all rather unmethylated Cs
                chr_tigs_distribution[current_chr_tigs_index][23] += 1
            elif current_average_met_probability > 0.4 and current_average_met_probability <= 0.6:
                # counts all 50/50-methylated Cs
                chr_tigs_distribution[current_chr_tigs_index][25] += 1
            elif current_average_met_probability > 0.6 and current_average_met_probability <= 0.8:
                # counts all rather methylated Cs
                chr_tigs_distribution[current_chr_tigs_index][27] += 1
            elif current_average_met_probability > 0.8 and current_average_met_probability <= 1:
                # counts all methylated Cs
                chr_tigs_distribution[current_chr_tigs_index][29] += 1

############################################################################################################################

# print('raw data is stored in all_5_mers')
# sort all_5_mers alphabetically
all_5_mers.sort()
# print('data is sorted')
# sort chr_tigs_distribution alphabetically
chr_tigs_distribution.sort()

print('-> data array was filled with absolute counts')

############################################################################################################################
############################################################################################################################

# data output
with open(fileCounts5mer, 'w') as outCounts5mer:
    for line in all_5_mers:
        outCounts5mer.write('\t'.join(str(x) for x in line) + '\n')
with open(fileCountsChr, 'w') as outCountsChr:
    for idx, line in enumerate(chr_tigs_distribution):
        # if idx != len(chr_tigs_distribution)-1:
            # because of multiple contigs the output must be structured differently for contig data
        outCountsChr.write('\t'.join(str(x) for x in line) +  '\n')
        # else:
            # outCountsChr.write('\t'.join(str(x) for x in line))

############################################################################################################################
############################################################################################################################

print('-> outputs were repeatrated')
print('-> script 2 ENDS')
print('Time needed: ' + str(datetime.now() - startTime))

############################################################################################################################
############################################################################################################################
