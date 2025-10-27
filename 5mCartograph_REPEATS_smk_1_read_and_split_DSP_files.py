# 5mCartograph_REPEATS_smk_1_read_and_split_DSP_files.py -> split and save the DSP output and repeat annotation file
# usable for all applications including TE specific data

import re
from datetime import datetime
import sys
import os

# start time
startTime = datetime.now()

# print('Hello, my name is ...')
print('-> hello, I am 5mCartograph_REPEATS REPEAT - an analysis tool for DeepSignal-plant data sets')
print('-> script 1 STARTS')

############################################################################################################################
############################################################################################################################

# use of the script:
# python3 [script used] [DSP result file path] [repeat annotation file path] [folder_name] [subfolder_name]

# DSP result file path
# format: chr, pos, strand, pos_in_strand, sum_prob_0 (unmet), sum_prob_1 (met), #_meth, #_unmet, coverage, met_frequency, k_mer
fileDSPout_all = sys.argv[1]

# repeat annotation file path
fileRepeatAnno_all = sys.argv[2]

# create location string for output 
file_out = sys.argv[3].strip()

# subfolder 
subfolder = sys.argv[4]

# file path log output
fileLOG_1 = file_out + subfolder + '5mCartograph_REPEATS_out1/5mCartograph_REPEATS_out1_log.txt'

############################################################################################################################
############################################################################################################################

# create result folder
os.makedirs(file_out + subfolder + '5mCartograph_REPEATS_out1/', exist_ok=True)

############################################################################################################################
############################################################################################################################

# write log head
with open(fileLOG_1, "w") as outLog:
    # return values of interest
    outLog.write('-> tool\n')
    outLog.write(__file__ + '\n')
    outLog.write('\n')

    outLog.write('-> description\n')
    outLog.write('in the first of three steps DeepSignal-plant output und repeat annotation are split\n')
    outLog.write('for parallelsized data processing\n')
    outLog.write('\n')

    outLog.write('-> DeepSignal-plant output as input\n')
    outLog.write(fileDSPout_all + '\n')
    outLog.write('\n')

    outLog.write('-> repeat annotation input\n')
    outLog.write(fileRepeatAnno_all + '\n')
    outLog.write('\n')

    outLog.write('->  log output\n')
    outLog.write(fileLOG_1 + '\n')
    outLog.write('\n')

############################################################################################################################
############################################################################################################################

# read in DSP/DS3 file 
with open(fileDSPout_all, "r") as inDSP:
    # string with elements that contain one row of text file each
    linesC = inDSP.readlines()
print('-> DeepSignal-plant output was read in -> ' + str(len(linesC)) + ' lines')

############################################################################################################################

# read in annotation file
with open(fileRepeatAnno_all, 'r') as inAnno:
    linesRepeatAnno = inAnno.readlines()
print('-> repeat annotation was read in')

############################################################################################################################
############################################################################################################################

chromosomes = []
# get list of chromosome names + tig      
for lineAnno in linesRepeatAnno:
    # test whether line holds chr/tig info
    if not lineAnno.startswith('#'):
        # get chrosomome from current line
        lineAnno_chr = re.split('\t', lineAnno)[0].strip()
        # check whether chromsome is already in list and if chromosome is chr
        if 'chr' in lineAnno_chr and lineAnno_chr not in chromosomes:
            chromosomes.append(lineAnno_chr)
# append listname for contigs
chromosomes.append('tig')
# sort list of chromosomes
chromosomes.sort()

# create list for checks
C_in_chr = [False for _ in range(len(chromosomes))]
chr_in_repeat_anno = [False for _ in range(len(chromosomes))]

############################################################################################################################
############################################################################################################################

# empty lists for file locations
fileDSPout_list = []
fileRepeatAnno_list = []
# create file paths for output
for chr in chromosomes:
    fileDSPout_list.append(file_out + subfolder + '5mCartograph_REPEATS_out1/5mCartograph_REPEATS_out1_frequency_' + chr + '.tsv')
    fileRepeatAnno_list.append(file_out + subfolder + '5mCartograph_REPEATS_out1/5mCartograph_REPEATS_out1_annotation_' + chr + '.gff')

############################################################################################################################
############################################################################################################################

# sort the lines into corresponding output files
for lineC in linesC:
    # get chromosome from current C
    lineC_chr = re.split('\t', lineC)[0].strip()
    # find correct list / file
    for idx, chromosome in enumerate(chromosomes):
        # store line in correct file
        if chromosome in lineC_chr:
            with open(fileDSPout_list[idx], 'a') as DSP_out:
                # set controll mechanism True -> there is data from this chromosome
                C_in_chr[idx] = True
                DSP_out.write(lineC)
            break

print('-> DeepSignal-plant output was split')

############################################################################################################################
############################################################################################################################

# sort the annotation lines into corresponding output files
for idx, lineRepeatAnno in enumerate(linesRepeatAnno):
    # test wether line is a real data line
    if not lineRepeatAnno.startswith('#') and lineRepeatAnno != '\n' and lineRepeatAnno != '':
        # get chrosomome from first entry line
        lineRepeatAnno_chr = re.split('\t', lineRepeatAnno)[0].strip()
        # find correct list / file
        for udx, chromosome in enumerate(chromosomes):
            # store line in corret file
            if chromosome in lineRepeatAnno_chr:
                with open(fileRepeatAnno_list[udx], 'a') as anno_out:
                    # set controll mechanism True -> there is data from this chromosome
                    chr_in_repeat_anno[udx] = True
                    anno_out.write(lineRepeatAnno)

print('-> repeat annotation was split')
                        
############################################################################################################################
############################################################################################################################

# create empty files for all chromosomes/tigs that are not in the dataset but exist in sugar beet
for idx, boolean in enumerate(C_in_chr):
    if boolean == False:
        with open(fileDSPout_list[idx], 'w') as DSP_out:
            DSP_out.write('')

for idx, boolean in enumerate(chr_in_repeat_anno):
    if boolean == False:
        with open(fileRepeatAnno_list[idx], 'w') as anno_out:
            anno_out.write('')

############################################################################################################################
############################################################################################################################

# write output file names to log file
with open(fileLOG_1, "a") as outLog:
    # print only the file names for split linesC and repeat annotation files that were used
    outLog.write('-> output\n')
    outLog.write('-> split DeepSignal-plant output\n')
    for idx, boolean in enumerate(C_in_chr):
        if boolean == True:
            outLog.write(fileDSPout_list[idx] + '\n')
        else:
            outLog.write(fileDSPout_list[idx] + ' -> no data\n')
    outLog.write('-> split repeat annotation\n')
    for idx, boolean in enumerate(chr_in_repeat_anno):
        if boolean == True:
            outLog.write(fileRepeatAnno_list[idx] + '\n')
        else:
            outLog.write(fileRepeatAnno_list[idx] + ' -> no data\n')
    outLog.write('\n')

    outLog.write('Number of lines: ' + str(len(linesC)) + '\n')
    outLog.write('\n')

print('-> script 1 ENDS')
print('Time needed: ' + str(datetime.now() - startTime))

############################################################################################################################
#############################################################################################################################