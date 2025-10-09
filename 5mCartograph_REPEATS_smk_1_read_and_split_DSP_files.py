# 5mCartograph_REPEATS_smk_1_read_and_split_DSP_files.py -> split and save the DSP output and repeat annotation file
# usable for all applications including TE specific data

import re
from datetime import datetime
import sys

# start time
startTime = datetime.now()

# print('Hello, my name is ...')
print('-> hello, I am 5mCartograph_REPEATS REPEAT - an analysis tool for DeepSignal-plant data sets')
print('-> script 1 STARTS')

############################################################################################################################
############################################################################################################################

# use of the script:
# python3 [script used] [DSP result file path] [repeat annotation file path] [subfolder]

# DSP result file path
# format: chr, pos, strand, pos_in_strand, sum_prob_0 (unmet), sum_prob_1 (met), #_meth, #_unmet, coverage, met_frequency, k_mer
fileDSPout_all = sys.argv[1]

# repeat annotation file path
fileRepeatAnno_all = sys.argv[2]

# create location string for output 
file_out = sys.argv[1].rsplit('/', 1)[0]

# subfolder 
subfolder = sys.argv[3]

# file path output
fileDSPout_chr1 = file_out + '/' + subfolder + '5mCartograph_REPEATS_out1/5mCartograph_REPEATS_out1_frequency_chr1.tsv'
fileDSPout_chr2 = file_out + '/' + subfolder + '5mCartograph_REPEATS_out1/5mCartograph_REPEATS_out1_frequency_chr2.tsv'
fileDSPout_chr3 = file_out + '/' + subfolder + '5mCartograph_REPEATS_out1/5mCartograph_REPEATS_out1_frequency_chr3.tsv'
fileDSPout_chr4 = file_out + '/' + subfolder + '5mCartograph_REPEATS_out1/5mCartograph_REPEATS_out1_frequency_chr4.tsv'
fileDSPout_chr5 = file_out + '/' + subfolder + '5mCartograph_REPEATS_out1/5mCartograph_REPEATS_out1_frequency_chr5.tsv'
fileDSPout_chr6 = file_out + '/' + subfolder + '5mCartograph_REPEATS_out1/5mCartograph_REPEATS_out1_frequency_chr6.tsv'
fileDSPout_chr7 = file_out + '/' + subfolder + '5mCartograph_REPEATS_out1/5mCartograph_REPEATS_out1_frequency_chr7.tsv'
fileDSPout_chr8 = file_out + '/' + subfolder + '5mCartograph_REPEATS_out1/5mCartograph_REPEATS_out1_frequency_chr8.tsv'
fileDSPout_chr9 = file_out + '/' + subfolder + '5mCartograph_REPEATS_out1/5mCartograph_REPEATS_out1_frequency_chr9.tsv'
fileDSPout_tigs = file_out + '/' + subfolder + '5mCartograph_REPEATS_out1/5mCartograph_REPEATS_out1_frequency_tigs.tsv'
fileDSPout_list = [fileDSPout_chr1, fileDSPout_chr2, fileDSPout_chr3, fileDSPout_chr4, fileDSPout_chr5, fileDSPout_chr6, fileDSPout_chr7, fileDSPout_chr8, fileDSPout_chr9, fileDSPout_tigs]
fileRepeatAnno_chr1 = file_out + '/' + subfolder + '5mCartograph_REPEATS_out1/5mCartograph_REPEATS_out1_repeat_annotation_chr1.gff'
fileRepeatAnno_chr2 = file_out + '/' + subfolder + '5mCartograph_REPEATS_out1/5mCartograph_REPEATS_out1_repeat_annotation_chr2.gff'
fileRepeatAnno_chr3 = file_out + '/' + subfolder + '5mCartograph_REPEATS_out1/5mCartograph_REPEATS_out1_repeat_annotation_chr3.gff'
fileRepeatAnno_chr4 = file_out + '/' + subfolder + '5mCartograph_REPEATS_out1/5mCartograph_REPEATS_out1_repeat_annotation_chr4.gff'
fileRepeatAnno_chr5 = file_out + '/' + subfolder + '5mCartograph_REPEATS_out1/5mCartograph_REPEATS_out1_repeat_annotation_chr5.gff'
fileRepeatAnno_chr6 = file_out + '/' + subfolder + '5mCartograph_REPEATS_out1/5mCartograph_REPEATS_out1_repeat_annotation_chr6.gff'
fileRepeatAnno_chr7 = file_out + '/' + subfolder + '5mCartograph_REPEATS_out1/5mCartograph_REPEATS_out1_repeat_annotation_chr7.gff'
fileRepeatAnno_chr8 = file_out + '/' + subfolder + '5mCartograph_REPEATS_out1/5mCartograph_REPEATS_out1_repeat_annotation_chr8.gff'
fileRepeatAnno_chr9 = file_out + '/' + subfolder + '5mCartograph_REPEATS_out1/5mCartograph_REPEATS_out1_repeat_annotation_chr9.gff'
fileRepeatAnno_tigs = file_out + '/' + subfolder + '5mCartograph_REPEATS_out1/5mCartograph_REPEATS_out1_repeat_annotation_tigs.gff'
fileRepeatAnno_list = [fileRepeatAnno_chr1, fileRepeatAnno_chr2, fileRepeatAnno_chr3, fileRepeatAnno_chr4, fileRepeatAnno_chr5, fileRepeatAnno_chr6, fileRepeatAnno_chr7, fileRepeatAnno_chr8, fileRepeatAnno_chr9, fileRepeatAnno_tigs]

# file path log output
fileLOG_1 = file_out + '/' + subfolder + '5mCartograph_REPEATS_out1/5mCartograph_REPEATS_out1_log.txt'

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

# chromosomes/tigs
chromosomes = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'tig']
C_in_chr = [False, False, False, False, False, False, False, False, False, False]
chr_in_repeat_anno = [False, False, False, False, False, False, False, False, False, False]

############################################################################################################################
############################################################################################################################

# read in file 
with open(fileDSPout_all, "r") as inDSP:
    # string with elements that contain one row of text file each
    linesC = inDSP.readlines()
print('-> DeepSignal-plant output was read in -> ' + str(len(linesC)) + ' lines')

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

# read in annotation file
with open(fileRepeatAnno_all, 'r') as inAnno:
    linesRepeatAnno = inAnno.readlines()
print('-> repeat annotation was read in')

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