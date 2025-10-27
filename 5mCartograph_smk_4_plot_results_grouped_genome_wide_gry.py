# 5mCartograph_plot_results.py -> plots results with matplotlib
# this script is designed for the plotting of genome-wide methylation rates
# gry = group related y axis
# for each group (5mC, CG, CHG, CHH) 100% means all cytosines of this group

# NEEDED: standard 5mCartograph
# /path/to/5mCartograph/5mCartograph_smk_3_summarize_calculate_and_output.py

import re
from datetime import datetime
import statistics
import sys
import numpy as np
import matplotlib.pyplot as plt

############################################################################################################################

# use
# python3 [script] [5mCartograph_result] [title]

############################################################################################################################

# data
result_file = sys.argv[1]

# create location string for output 
folder = sys.argv[1].rsplit('/', 1)[0]

# log file
log_file = folder + '/5mCartograph_plot_results_grouped_genome_wide_gry_log.txt'

# title
title = sys.argv[2]
title = title.replace('_', ' ')

# gets starttime
startTime = datetime.now()

############################################################################################################################

# read 5mCartograph result file (based on DeepSignal-plant output)
with open(result_file, 'r') as rf:
    result_lines = rf.readlines()

############################################################################################################################

# get exon count and gene length for each set
def get_data_of_interest(lines):

    # lists for C, CG, CHH and CHG for 5 values each 
    values = [None]*20
    C_p = None
    CG_p = None
    CHG_p = None
    CHH_p = None

    # cycle through list and get % values that correspond to context
    for idx, line in enumerate(lines):
        # all Cs
        C_p = 100
        if line.startswith('highly unmethylated cytosines: ', 0):
            values[0] = float(line.split(' ')[-2])
        elif line.startswith('rather unmethylated cytosines: ', 0):
            values[1] = float(line.split(' ')[-2])
        elif line.startswith('50/50-methylated cytosines: ', 0):
            values[2] = float(line.split(' ')[-2])
        elif line.startswith('rather methylated cytosines: ', 0):
            values[3] = float(line.split(' ')[-2])
        elif line.startswith('highly methylated cytosines: ', 0):
            values[4] = float(line.split(' ')[-2])
        # CG
        elif line.startswith('CGs: ', 0):
            CG_p = float(line.split(' ')[-2])
        elif line.startswith('highly unmethylated CGs: ', 0):
            values[5] = float(line.split(' ')[-2])
        elif line.startswith('rather unmethylated CGs: ', 0):
            values[6] = float(line.split(' ')[-2])
        elif line.startswith('50/50-methylated CGs: ', 0):
            values[7] = float(line.split(' ')[-2])
        elif line.startswith('rather methylated CGs: ', 0):
            values[8] = float(line.split(' ')[-2])
        elif line.startswith('highly methylated CGs: ', 0):
            values[9] = float(line.split(' ')[-2])
        # CHG
        elif line.startswith('CHGs: ', 0):
            CHG_p = float(line.split(' ')[-2])
        elif line.startswith('highly unmethylated CHGs: ', 0):
            values[10] = float(line.split(' ')[-2])
        elif line.startswith('rather unmethylated CHGs: ', 0):
            values[11] = float(line.split(' ')[-2])
        elif line.startswith('50/50-methylated CHGs: ', 0):
            values[12] = float(line.split(' ')[-2])
        elif line.startswith('rather methylated CHGs: ', 0):
            values[13] = float(line.split(' ')[-2])
        elif line.startswith('highly methylated CHGs: ', 0):
            values[14] = float(line.split(' ')[-2])
        # CHH
        elif line.startswith('CHHs: ', 0):
            CHH_p = float(line.split(' ')[-2])
        elif line.startswith('highly unmethylated CHHs: ', 0):
            values[15] = float(line.split(' ')[-2])
        elif line.startswith('rather unmethylated CHHs: ', 0):
            values[16] = float(line.split(' ')[-2])
        elif line.startswith('50/50-methylated CHHs: ', 0):
            values[17] = float(line.split(' ')[-2])
        elif line.startswith('rather methylated CHHs: ', 0):
            values[18] = float(line.split(' ')[-2])
        elif line.startswith('highly methylated CHHs: ', 0):
            values[19] = float(line.split(' ')[-2])
        

    return values

# get values def
values = get_data_of_interest(result_lines)

# create list for every context/group
C_values = values[0:5]
CG_values = values[5:10]
CHG_values = values[10:15]
CHH_values = values[15:20]
print(C_values)
print(CG_values)
print(CHG_values)
print(CHH_values)

############################################################################################################################

# print maxima to file
with open(log_file, 'w') as fileOut:
    fileOut.write('-> tool\n')
    fileOut.write(__file__ + '\n')
    fileOut.write('\n')
    fileOut.write('-> input\n')
    fileOut.write(result_file + '\n')
    fileOut.write('\n')
    fileOut.write('-> output\n')
    fileOut.write('5mCartograph - genome-wide 5mC, CG, CHG and CHH methylation probability distribution figure (grouped)\n')
    fileOut.write('probabilities are in relation to all cytosines in one context\n')
    fileOut.write(log_file + '\n')
    fileOut.write('\n')
    fileOut.write('-> results\n')
    fileOut.write('-> methylation probability distribution [%]\n')
    fileOut.write('5mC: \t' + str(C_values) + '\n')
    fileOut.write('CG: \t' + str(CG_values) + '\n')
    fileOut.write('CHG: \t' + str(CHG_values) + '\n')
    fileOut.write('CHH: \t' + str(CHH_values) + '\n')
    fileOut.write('\n')

############################################################################################################################

# colums and rows for the plot
meth_prob_label = ['highly unmethylated', 'rather unmethylated', '50/50 methylated', 'rather methylated', 'highly methylated']
colums = ['all C', 'CG', 'CHG', 'CHH']
x_pos = [1, 2, 3, 4]
colors = ['forestgreen', 'limegreen', 'paleturquoise', 'deepskyblue', 'dodgerblue']
width = 0.16

# data sets for the plot
probabilities = np.array([C_values, CG_values, CHG_values, CHH_values])
print(probabilities)

############################################################################################################################

fig, ax = plt.subplots(figsize = (8, 5))
ax.set_xticks(x_pos, colums)

for idx, context in enumerate(probabilities):
    for udx, meth_prob_value in enumerate(context):
        if idx == 0:
            p = ax.bar((x_pos[idx]-0.32+(udx*0.16)), meth_prob_value, width, label=meth_prob_label[udx], color=colors[udx])
        else:
            p = ax.bar((x_pos[idx]-0.32+(udx*0.16)), meth_prob_value, width, color=colors[udx])
        # plot bar heights
        ax.bar_label(p, padding=3, fmt='%.1f', fontsize=8)

ax.set_ylim(top=120)
ax.set_yticks([0, 20, 40, 60, 80, 100])
ax.set_ylabel('proportion of cytosines\nper context [%]', fontsize='x-large')
ax.set_xticklabels(colums, fontsize='x-large')
ax.set_title(title, fontsize='xx-large')
ax.legend(loc='upper right', fontsize='large')
ax.spines[['right', 'top']].set_visible(False)
ax.spines['left'].set_bounds(0, 100)
ax.vlines(1.5, 0, 100, linestyles='dashed', color='black')

# plt.savefig(folder+'/5mCartograph_plot_results_grouped_genome_wide_gry.png', bbox_inches='tight')
plt.savefig(folder+'/5mCartograph_plot_results_grouped_genome_wide_gry_bar_hight.png', bbox_inches='tight')
plt.show()

############################################################################################################################

# returns time needed for comparisons -> 0:00:00.001678 -> large set should be done in about an hour
print(datetime.now() - startTime)
