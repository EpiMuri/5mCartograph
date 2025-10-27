# 5mCartograph_smk_2_chr_specific_counts.py -> calculate 5mer and chromosome specific data (the chromosomes are processed individually, while the contigs run in the 10th slot)


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

# use of the script:
# python3 [script used] [partial DSP result file path] [partial annotation file path] [output_folder]

# get chromosome for file names
chrom = re.split('_', sys.argv[1])[-1].strip()
chrom = re.split('\.', chrom)[0]

# DSP result file path for one chromosome / all tigs
fileDSPout = sys.argv[1]

# annotation file path for one chromosome / all tigs
fileAnno = sys.argv[2]

# create location string for output 
file_out = sys.argv[3].strip()

# file path output
fileCounts5mer = file_out + '/5mCartograph_out2_' + chrom + '/5mCartograph_out2_counts_per_5mer_' + chrom + '.txt'
fileCountsChr = file_out + '/5mCartograph_out2_' + chrom + '/5mCartograph_out2_counts_per_chromosome_' + chrom + '.txt'

# file path log output
fileLOG_2 = file_out + '/5mCartograph_out2_' + chrom + '/5mCartograph_out2_log_' + chrom + '.txt'

############################################################################################################################
############################################################################################################################

# create result folder
os.makedirs(file_out + '/5mCartograph_out2_' + chrom + '/', exist_ok=True)

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

    outLog.write('-> annotation input\n')
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

# read in annotation file
with open(fileAnno, 'r') as inAnno:
    linesAnno = inAnno.readlines()
print('-> annotation was read in')

############################################################################################################################
############################################################################################################################

# create giant list of lists of lists to manage the anno file

# giant list contains sublists, one sublist for each gene
# structure giantAnnoList
# [ [[], [], ...], [[], [], ...], [[], [], ...], ... ] 
giantAnnoList = []

# each sublist contains lists of the line of the gene, mRNA, exons and CDSs
# structure geneSublist
# [geneline, mRNAline, exonline, CDSline, exonline, CDSline, ... ]
geneSublist = []

# each subsublist contains attributes of the gene/mRNA/exon/CDS
# structure attributeSublist
# [chr, classification (mRNA, exon, CDS), lower border, upper border]
attributeSublist = []
# attributeSublist for gene (needs to be appended last)
attributeSublistGene = []

# cycles through all lines of fileAnno to find genes, then all attributes of the gene and its subunits are collected and stored in giantAnnoList
# cycles through linesAnno
for idx,lineAnno in enumerate(linesAnno):

    # makes sure that data lines are analysed only
    if (lineAnno[0] != '#') and (len(re.split('\t', lineAnno)) == 9):

        # finds all lines with gene data
        if re.split('\t', lineAnno)[2].strip() == 'gene' or re.split('\t', lineAnno)[2].strip() == 'region':

            # cycles through all lines from the current line till the end of the file ...
            for i in range(idx, len(linesAnno)):
                
                # ... until the end of the gene block
                if '##' in linesAnno[i]:
                    
                    # if block ends the current geneSublist is sorted and subjected to giantAnnoList
                    geneSublist.sort()
                    geneSublist.append(attributeSublistGene)
                    giantAnnoList.append(geneSublist)
                    geneSublist = []
                    break
                # if the iteration is not at the end of a block, the attributes of each gene/exon/CDS are stored in attribute Sublist and subjected to geneSublist
                else:
                    # all attributes of the current lineAnno are stored in attributeSublist
                    attributeSublist = [re.split('\t', linesAnno[i])[0].strip(), re.split('\t', linesAnno[i])[2].strip(), int(re.split('\t', linesAnno[i])[3].strip()), int(re.split('\t', linesAnno[i])[4].strip())]
                    
                    # if attributeSublist is for gene it needs to be appended as last element and is therefore stored individually
                    if attributeSublist[1] == 'gene' or  attributeSublist[1] == 'region':
                        # current attributeSublist is subjected to attributeSublistGene
                        attributeSublistGene = attributeSublist
                    # if attribute Sublist is not for gene it can be appended regularly
                    else:
                        # attributeSublist is appended to geneSublist
                        geneSublist.append(attributeSublist)
                    # attributeSublist is cleared so that it can be filled again with the next attributes
                    attributeSublist = []

print('-> giant list of genes was created')

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
# the 5-mers are categoriesed in in two ways:
# unmethylated, rather unmethylated, 50/50-methylated, rather methylated, methylated
# non_genic, intron, CDS and exon - but not CDS

# find all C-positions and corresponding k-mers -> go through all lines
for idx, lineC in enumerate(linesC):
    
    # as a progress controll in every 1000000th cycle the cycle number is printed
    if idx % 1000000 == 0:
        # print(idx)
        with open(fileLOG_2, "a") as out2:
            # final return
            out2.write(str(idx) + '\n')
                
    # if line is not empty, line is analysed
    if len(re.split('\t', lineC)) == columns:

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
        # in different genetic, non_genic regions

        # count data lines
        c_counter += 1
        
        # get index to of the corresponding context
        current_5_mer_index = contexts.index(k_mer_sequence)
        
        # analyse whether the C is in a non_genic region, in an intron, in CDS or in UTR (+ upstream + downstream)
        # non_genic -> non_genic
        # genic -> intron + UTR + CDS
        # exon -> UTR + CDS
        # CDS -> CDS
        # intron -> intron
        # upstream -> upstream
        # downstream -> downstream
        
        # loop to cycle through all genes of annotation
        # if C not in genes -> non_genic [31-40]
        # if C not exons -> intron [21-30]
        # if C not in CDS -> UTR [11-20]
        # if C in CDS -> CDS [1-10]
        # if C 1 kb upstream -> upstream [41-50]
        # if C 1 kb downstream -> downstream [51-60]

        # to check whether C has already been sorted to a group (CDS, UTR, intron or non-genic)
        is_in_gene = False
        is_in_CDS = False
        is_in_exon = False
        is_in_intron = False
        is_in_1kb_upstream = False
        is_in_1kb_downstream = False

        # cycles through all genes in giantAnnoList
        for gene in giantAnnoList:

            # check whether current C is in current gene (chromosomes and borders)
            if (gene[-2][0] == current_chr_tigs) and (current_C_position >= gene[-2][2]) and (current_C_position <= gene[-2][3]):

                # is in a gene!
                # gene values are later on calculated by addition of CDS, UTR and intron values
                is_in_gene = True
                
                # cycles through entries of the current gene to check, whether current C is in exon and CDS
                for entry in gene:
                    
                    # dicrimination between CDS, UTR and intron
                    # exon values are later on calculated by addition of CDS and UTR values

                    # checks whether current C is in current genes CDS
                    if (entry[1] == 'CDS') and (current_C_position >= entry[2]) and (current_C_position <= entry[3]):
                        
                        # is in this cycles CDS!
                        is_in_CDS = True

                    # checks whether current C is in current genes exons if it is not in CDS -> UTR
                    if (entry[1] == 'exon') and (current_C_position >= entry[2]) and (current_C_position <= entry[3]):

                        # is in this cylces Exon!
                        is_in_exon = True
    
                # if current C is in current gene, but not in CDS and exon -> intron
                if (is_in_CDS == False) and (is_in_exon == False):

                    # is in intron!
                    is_in_intron = True
            
            # check whether current C is in 1 kb before current gene (upstream)
            elif (gene[-2][0] == current_chr_tigs) and (current_C_position >= (gene[-2][2] - 1000)) and (current_C_position < gene[-2][2]):

                # is in 1 kb upstream!
                is_in_1kb_upstream = True

            # check whether current C is in 1 kb after current gene (downstream)
            elif (gene[-2][0] == current_chr_tigs) and (current_C_position > gene[-2][3]) and (current_C_position <= (gene[-2][3] + 1000)):
                # is in 1 kb downstream!
                is_in_1kb_downstream = True

        # check whether C is in a CDS
        if is_in_CDS == True:

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

        # check whether C is in a UTR (in exon, but not in CDS)
        elif is_in_exon == True:

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

        # check whether C is in a gene
        elif is_in_intron == True:

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

        # check whether C is in 1 kb upstream of genic region
        elif is_in_1kb_upstream == True:

            # discriminate between unmethylated and methylated Cs -> 5-mers
            if current_average_met_probability >= 0 and current_average_met_probability <= 0.2:
                # counts all unmethalted Cs
                all_5_mers[current_5_mer_index][41] += 1
            elif current_average_met_probability > 0.2 and current_average_met_probability <= 0.4:
                # counts all rather unmethylated Cs
                all_5_mers[current_5_mer_index][43] += 1
            elif current_average_met_probability > 0.4 and current_average_met_probability <= 0.6:
                # counts all 50/50-methylated Cs
                all_5_mers[current_5_mer_index][45] += 1
            elif current_average_met_probability > 0.6 and current_average_met_probability <= 0.8:
                # counts all rather methylated Cs
                all_5_mers[current_5_mer_index][47] += 1
            elif current_average_met_probability > 0.8 and current_average_met_probability <= 1:
                # counts all methylated Cs
                all_5_mers[current_5_mer_index][49] += 1 

        # check whether C is in 1 kb downstream of genic region
        elif is_in_1kb_downstream == True:

            # discriminate between unmethylated and methylated Cs -> 5-mers
            if current_average_met_probability >= 0 and current_average_met_probability <= 0.2:
                # counts all unmethalted Cs
                all_5_mers[current_5_mer_index][51] += 1
            elif current_average_met_probability > 0.2 and current_average_met_probability <= 0.4:
                # counts all rather unmethylated Cs
                all_5_mers[current_5_mer_index][53] += 1
            elif current_average_met_probability > 0.4 and current_average_met_probability <= 0.6:
                # counts all 50/50-methylated Cs
                all_5_mers[current_5_mer_index][55] += 1
            elif current_average_met_probability > 0.6 and current_average_met_probability <= 0.8:
                # counts all rather methylated Cs
                all_5_mers[current_5_mer_index][57] += 1
            elif current_average_met_probability > 0.8 and current_average_met_probability <= 1:
                # counts all methylated Cs
                all_5_mers[current_5_mer_index][59] += 1 

        # check whether C is in a genic region
        elif is_in_gene == False:
            
            # C is in a non genic region

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

# sort all_5_mers alphabetically
all_5_mers.sort()
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
        outCountsChr.write('\t'.join(str(x) for x in line) +  '\n')

############################################################################################################################
############################################################################################################################

print('-> outputs were generated')
print('-> script 2 ENDS')
print('Time needed: ' + str(datetime.now() - startTime))

############################################################################################################################
############################################################################################################################
