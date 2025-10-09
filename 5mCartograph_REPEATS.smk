from datetime import datetime

# snakefile for the analysis of 5mC in the whole genome of sugar beet for REPEATS
# the original basic_5mC_analysis python script was split into three subscripts
# the script orginally run 21 days -> now it is about one days due to parallelization
# the tool uses splitting of lines instead of the dataframes, as it runs faster
# it is not able to calculate a score for the quality of the methylation calling yet

# FILES NEEDED:
# 5mCartograph_REPEATS.smk -> parallelization of subscript 2 and coordination of subscripts
# Config_5mCargotraph_....yaml -> configuration of sample specific files (individual for each sample)
# 5mCartograph_REPEATS_smk_1_read_and_split_DSP_files.py -> split and save the DSP output and repeat annotation file
# 5mCartograph_REPEATS_smk_2_chr_specific_counts.py -> calculate 5mer and chromosome specific data (the chromosomes are processed individually, while the contigs run in the 10th slot)
# 5mCartograph_REPEATS_smk_3_summarize_calculate_and_output.py -> summarizes data of subsets and calculates data of interest
# 5mCartograph_REPEATS_smk_4_plot_results.py -> plots results with matplotlib

# tsv needs to be placed in destination directory
# snakefile needs to be used in 'snakemake' environment 
# Command: snakemake all --cores 10 -s 5mCartograph_REPEATS.smk

############################################################################################################################
############################################################################################################################

# harbours specific information about the data sets used
configfile: '/path/to/scripts/config_file'

# variables
folder_name = config['folder_name']
subfolder_name = config['subfolder_name']
tsv_name = config['tsv_name']
anno_name = config['anno_name']
plot_title = config['plot_title']

# for efficiency check
startTime = datetime.now()

############################################################################################################################
############################################################################################################################

rule split:
    input:
        DSP = expand("{folder}{tsv}", folder=folder_name, tsv=tsv_name),
        REPEATS_annotation = expand("{anno}", anno=anno_name)
    output:
        dir_1 = directory("{folder}{subfolder_name}5mCartograph_REPEATS_out1/"),
        DSP_chr1 = "{folder}{subfolder_name}5mCartograph_REPEATS_out1/5mCartograph_REPEATS_out1_frequency_chr1.tsv",
        DSP_chr2 = "{folder}{subfolder_name}5mCartograph_REPEATS_out1/5mCartograph_REPEATS_out1_frequency_chr2.tsv",
        DSP_chr3 = "{folder}{subfolder_name}5mCartograph_REPEATS_out1/5mCartograph_REPEATS_out1_frequency_chr3.tsv",
        DSP_chr4 = "{folder}{subfolder_name}5mCartograph_REPEATS_out1/5mCartograph_REPEATS_out1_frequency_chr4.tsv",
        DSP_chr5 = "{folder}{subfolder_name}5mCartograph_REPEATS_out1/5mCartograph_REPEATS_out1_frequency_chr5.tsv",
        DSP_chr6 = "{folder}{subfolder_name}5mCartograph_REPEATS_out1/5mCartograph_REPEATS_out1_frequency_chr6.tsv",
        DSP_chr7 = "{folder}{subfolder_name}5mCartograph_REPEATS_out1/5mCartograph_REPEATS_out1_frequency_chr7.tsv",
        DSP_chr8 = "{folder}{subfolder_name}5mCartograph_REPEATS_out1/5mCartograph_REPEATS_out1_frequency_chr8.tsv",
        DSP_chr9 = "{folder}{subfolder_name}5mCartograph_REPEATS_out1/5mCartograph_REPEATS_out1_frequency_chr9.tsv",
        DSP_tigs = "{folder}{subfolder_name}5mCartograph_REPEATS_out1/5mCartograph_REPEATS_out1_frequency_tigs.tsv",
        REPEATS_annotation_chr1 = "{folder}{subfolder_name}5mCartograph_REPEATS_out1/5mCartograph_REPEATS_out1_repeat_annotation_chr1.gff",
        REPEATS_annotation_chr2 = "{folder}{subfolder_name}5mCartograph_REPEATS_out1/5mCartograph_REPEATS_out1_repeat_annotation_chr2.gff",
        REPEATS_annotation_chr3 = "{folder}{subfolder_name}5mCartograph_REPEATS_out1/5mCartograph_REPEATS_out1_repeat_annotation_chr3.gff",
        REPEATS_annotation_chr4 = "{folder}{subfolder_name}5mCartograph_REPEATS_out1/5mCartograph_REPEATS_out1_repeat_annotation_chr4.gff",
        REPEATS_annotation_chr5 = "{folder}{subfolder_name}5mCartograph_REPEATS_out1/5mCartograph_REPEATS_out1_repeat_annotation_chr5.gff",
        REPEATS_annotation_chr6 = "{folder}{subfolder_name}5mCartograph_REPEATS_out1/5mCartograph_REPEATS_out1_repeat_annotation_chr6.gff",
        REPEATS_annotation_chr7 = "{folder}{subfolder_name}5mCartograph_REPEATS_out1/5mCartograph_REPEATS_out1_repeat_annotation_chr7.gff",
        REPEATS_annotation_chr8 = "{folder}{subfolder_name}5mCartograph_REPEATS_out1/5mCartograph_REPEATS_out1_repeat_annotation_chr8.gff",
        REPEATS_annotation_chr9 = "{folder}{subfolder_name}5mCartograph_REPEATS_out1/5mCartograph_REPEATS_out1_repeat_annotation_chr9.gff",
        REPEATS_annotation_tigs = "{folder}{subfolder_name}5mCartograph_REPEATS_out1/5mCartograph_REPEATS_out1_repeat_annotation_tigs.gff"
    params:
        subfolder_name = config['subfolder_name']
    shell:
        """
        mkdir -p {output.dir_1}
        python3 /path/to/scripts/5mCartograph_REPEATS/5mCartograph_REPEATS_smk_1_read_and_split_DSP_files.py {input.DSP} {input.REPEATS_annotation} {params.subfolder_name}
        # python3 /path/to/5mCartograph_repeats/scripts/5mCartograph_REPEATS_smk_1_read_and_split_DSP_files.py {input.DSP} {input.REPEATS_annotation} {params.subfolder_name}
        """

############################################################################################################################
############################################################################################################################

rule calculate_chr1:
    input:
        DSP_chr1 = "{folder}{subfolder_name}5mCartograph_REPEATS_out1/5mCartograph_REPEATS_out1_frequency_chr1.tsv",
        REPEATS_annotation_chr1 = "{folder}{subfolder_name}5mCartograph_REPEATS_out1/5mCartograph_REPEATS_out1_repeat_annotation_chr1.gff"
    output:
        dir_2 = directory("{folder}{subfolder_name}5mCartograph_REPEATS_out2_chr1/"),
        result_5mer_chr1 = "{folder}{subfolder_name}5mCartograph_REPEATS_out2_chr1/5mCartograph_REPEATS_out2_counts_per_5mer_chr1.txt",
        result_chr_chr1 = "{folder}{subfolder_name}5mCartograph_REPEATS_out2_chr1/5mCartograph_REPEATS_out2_counts_per_chromosome_chr1.txt"
    shell:
        """
        mkdir -p {output.dir_2}
        python3 /path/to/scripts/5mCartograph_REPEATS/5mCartograph_REPEATS_smk_2_chr_specific_counts.py {input.DSP_chr1} {input.REPEATS_annotation_chr1}
        # python3 /path/to/5mCartograph_repeats/scripts/5mCartograph_REPEATS_smk_2_chr_specific_counts.py {input.DSP_chr1} {input.REPEATS_annotation_chr1}
        """

############################################################################################################################

rule calculate_chr2:
    input:
        DSP_chr2 = "{folder}{subfolder_name}5mCartograph_REPEATS_out1/5mCartograph_REPEATS_out1_frequency_chr2.tsv",
        REPEATS_annotation_chr2 = "{folder}{subfolder_name}5mCartograph_REPEATS_out1/5mCartograph_REPEATS_out1_repeat_annotation_chr2.gff"
    output:
        dir_2 = directory("{folder}{subfolder_name}5mCartograph_REPEATS_out2_chr2/"),
        result_5mer_chr2 = "{folder}{subfolder_name}5mCartograph_REPEATS_out2_chr2/5mCartograph_REPEATS_out2_counts_per_5mer_chr2.txt",
        result_chr_chr2 = "{folder}{subfolder_name}5mCartograph_REPEATS_out2_chr2/5mCartograph_REPEATS_out2_counts_per_chromosome_chr2.txt"
    shell:
        """
        mkdir -p {output.dir_2}
        python3 /path/to/scripts/5mCartograph_REPEATS/5mCartograph_REPEATS_smk_2_chr_specific_counts.py {input.DSP_chr2} {input.REPEATS_annotation_chr2}
        # python3 /path/to/5mCartograph_repeats/scripts/5mCartograph_REPEATS_smk_2_chr_specific_counts.py {input.DSP_chr2} {input.REPEATS_annotation_chr2}
        """

############################################################################################################################

rule calculate_chr3:
    input:
        DSP_chr3 = "{folder}{subfolder_name}5mCartograph_REPEATS_out1/5mCartograph_REPEATS_out1_frequency_chr3.tsv",
        REPEATS_annotation_chr3 = "{folder}{subfolder_name}5mCartograph_REPEATS_out1/5mCartograph_REPEATS_out1_repeat_annotation_chr3.gff"
    output:
        dir_2 = directory("{folder}{subfolder_name}5mCartograph_REPEATS_out2_chr3/"),
        result_5mer_chr3 = "{folder}{subfolder_name}5mCartograph_REPEATS_out2_chr3/5mCartograph_REPEATS_out2_counts_per_5mer_chr3.txt",
        result_chr_chr3 = "{folder}{subfolder_name}5mCartograph_REPEATS_out2_chr3/5mCartograph_REPEATS_out2_counts_per_chromosome_chr3.txt"
    shell:
        """
        mkdir -p {output.dir_2}
        python3 /path/to/scripts/5mCartograph_REPEATS/5mCartograph_REPEATS_smk_2_chr_specific_counts.py {input.DSP_chr3} {input.REPEATS_annotation_chr3}
        # python3 /path/to/5mCartograph_repeats/scripts/5mCartograph_REPEATS_smk_2_chr_specific_counts.py {input.DSP_chr3} {input.REPEATS_annotation_chr3}
        """

############################################################################################################################

rule calculate_chr4:
    input:
        DSP_chr4 = "{folder}{subfolder_name}5mCartograph_REPEATS_out1/5mCartograph_REPEATS_out1_frequency_chr4.tsv",
        REPEATS_annotation_chr4 = "{folder}{subfolder_name}5mCartograph_REPEATS_out1/5mCartograph_REPEATS_out1_repeat_annotation_chr4.gff"
    output:
        dir_2 = directory("{folder}{subfolder_name}5mCartograph_REPEATS_out2_chr4/"),
        result_5mer_chr4 = "{folder}{subfolder_name}5mCartograph_REPEATS_out2_chr4/5mCartograph_REPEATS_out2_counts_per_5mer_chr4.txt",
        result_chr_chr4 = "{folder}{subfolder_name}5mCartograph_REPEATS_out2_chr4/5mCartograph_REPEATS_out2_counts_per_chromosome_chr4.txt"
    shell:
        """
        mkdir -p {output.dir_2}
        python3 /path/to/scripts/5mCartograph_REPEATS/5mCartograph_REPEATS_smk_2_chr_specific_counts.py {input.DSP_chr4} {input.REPEATS_annotation_chr4}
        # python3 /path/to/5mCartograph_repeats/scripts/5mCartograph_REPEATS_smk_2_chr_specific_counts.py {input.DSP_chr4} {input.REPEATS_annotation_chr4}
        """

############################################################################################################################

rule calculate_chr5:
    input:
        DSP_chr5 = "{folder}{subfolder_name}5mCartograph_REPEATS_out1/5mCartograph_REPEATS_out1_frequency_chr5.tsv",
        REPEATS_annotation_chr5 = "{folder}{subfolder_name}5mCartograph_REPEATS_out1/5mCartograph_REPEATS_out1_repeat_annotation_chr5.gff"
    output:
        dir_2 = directory("{folder}{subfolder_name}5mCartograph_REPEATS_out2_chr5/"),
        result_5mer_chr5 = "{folder}{subfolder_name}5mCartograph_REPEATS_out2_chr5/5mCartograph_REPEATS_out2_counts_per_5mer_chr5.txt",
        result_chr_chr5 = "{folder}{subfolder_name}5mCartograph_REPEATS_out2_chr5/5mCartograph_REPEATS_out2_counts_per_chromosome_chr5.txt"
    shell:
        """
        mkdir -p {output.dir_2}
        python3 /path/to/scripts/5mCartograph_REPEATS/5mCartograph_REPEATS_smk_2_chr_specific_counts.py {input.DSP_chr5} {input.REPEATS_annotation_chr5}
        # python3 /path/to/5mCartograph_repeats/scripts/5mCartograph_REPEATS_smk_2_chr_specific_counts.py {input.DSP_chr5} {input.REPEATS_annotation_chr5}
        """

############################################################################################################################

rule calculate_chr6:
    input:
        DSP_chr6 = "{folder}{subfolder_name}5mCartograph_REPEATS_out1/5mCartograph_REPEATS_out1_frequency_chr6.tsv",
        REPEATS_annotation_chr6 = "{folder}{subfolder_name}5mCartograph_REPEATS_out1/5mCartograph_REPEATS_out1_repeat_annotation_chr6.gff"
    output:
        dir_2 = directory("{folder}{subfolder_name}5mCartograph_REPEATS_out2_chr6/"),
        result_5mer_chr6 = "{folder}{subfolder_name}5mCartograph_REPEATS_out2_chr6/5mCartograph_REPEATS_out2_counts_per_5mer_chr6.txt",
        result_chr_chr6 = "{folder}{subfolder_name}5mCartograph_REPEATS_out2_chr6/5mCartograph_REPEATS_out2_counts_per_chromosome_chr6.txt"
    shell:
        """
        mkdir -p {output.dir_2}
        python3 /path/to/scripts/5mCartograph_REPEATS/5mCartograph_REPEATS_smk_2_chr_specific_counts.py {input.DSP_chr6} {input.REPEATS_annotation_chr6}
        # python3 /path/to/5mCartograph_repeats/scripts/5mCartograph_REPEATS_smk_2_chr_specific_counts.py {input.DSP_chr6} {input.REPEATS_annotation_chr6}
        """

############################################################################################################################

rule calculate_chr7:
    input:
        DSP_chr7 = "{folder}{subfolder_name}5mCartograph_REPEATS_out1/5mCartograph_REPEATS_out1_frequency_chr7.tsv",
        REPEATS_annotation_chr7 = "{folder}{subfolder_name}5mCartograph_REPEATS_out1/5mCartograph_REPEATS_out1_repeat_annotation_chr7.gff"
    output:
        dir_2 = directory("{folder}{subfolder_name}5mCartograph_REPEATS_out2_chr7/"),
        result_5mer_chr7 = "{folder}{subfolder_name}5mCartograph_REPEATS_out2_chr7/5mCartograph_REPEATS_out2_counts_per_5mer_chr7.txt",
        result_chr_chr7 = "{folder}{subfolder_name}5mCartograph_REPEATS_out2_chr7/5mCartograph_REPEATS_out2_counts_per_chromosome_chr7.txt"
    shell:
        """
        mkdir -p {output.dir_2}
        python3 /path/to/scripts/5mCartograph_REPEATS/5mCartograph_REPEATS_smk_2_chr_specific_counts.py {input.DSP_chr7} {input.REPEATS_annotation_chr7}
        # python3 /path/to/5mCartograph_repeats/scripts/5mCartograph_REPEATS_smk_2_chr_specific_counts.py {input.DSP_chr7} {input.REPEATS_annotation_chr7}
        """

############################################################################################################################

rule calculate_chr8:
    input:
        DSP_chr8 = "{folder}{subfolder_name}5mCartograph_REPEATS_out1/5mCartograph_REPEATS_out1_frequency_chr8.tsv",
        REPEATS_annotation_chr8 = "{folder}{subfolder_name}5mCartograph_REPEATS_out1/5mCartograph_REPEATS_out1_repeat_annotation_chr8.gff"
    output:
        dir_2 = directory("{folder}{subfolder_name}5mCartograph_REPEATS_out2_chr8/"),
        result_5mer_chr8 = "{folder}{subfolder_name}5mCartograph_REPEATS_out2_chr8/5mCartograph_REPEATS_out2_counts_per_5mer_chr8.txt",
        result_chr_chr8 = "{folder}{subfolder_name}5mCartograph_REPEATS_out2_chr8/5mCartograph_REPEATS_out2_counts_per_chromosome_chr8.txt"
    shell:
        """
        mkdir -p {output.dir_2}
        python3 /path/to/scripts/5mCartograph_REPEATS/5mCartograph_REPEATS_smk_2_chr_specific_counts.py {input.DSP_chr8} {input.REPEATS_annotation_chr8}
        # python3 /path/to/5mCartograph_repeats/scripts/5mCartograph_REPEATS_smk_2_chr_specific_counts.py {input.DSP_chr8} {input.REPEATS_annotation_chr8}
        """

############################################################################################################################

rule calculate_chr9:
    input:
        DSP_chr9 = "{folder}{subfolder_name}5mCartograph_REPEATS_out1/5mCartograph_REPEATS_out1_frequency_chr9.tsv",
        REPEATS_annotation_chr9 = "{folder}{subfolder_name}5mCartograph_REPEATS_out1/5mCartograph_REPEATS_out1_repeat_annotation_chr9.gff"
    output:
        dir_2 = directory("{folder}{subfolder_name}5mCartograph_REPEATS_out2_chr9/"),
        result_5mer_chr9 = "{folder}{subfolder_name}5mCartograph_REPEATS_out2_chr9/5mCartograph_REPEATS_out2_counts_per_5mer_chr9.txt",
        result_chr_chr9 = "{folder}{subfolder_name}5mCartograph_REPEATS_out2_chr9/5mCartograph_REPEATS_out2_counts_per_chromosome_chr9.txt"
    shell:
        """
        mkdir -p {output.dir_2}
        python3 /path/to/scripts/5mCartograph_REPEATS/5mCartograph_REPEATS_smk_2_chr_specific_counts.py {input.DSP_chr9} {input.REPEATS_annotation_chr9}
        # python3 /path/to/5mCartograph_repeats/scripts/5mCartograph_REPEATS_smk_2_chr_specific_counts.py {input.DSP_chr9} {input.REPEATS_annotation_chr9}
        """

############################################################################################################################

rule calculate_tigs:
    input:
        DSP_tigs = "{folder}{subfolder_name}5mCartograph_REPEATS_out1/5mCartograph_REPEATS_out1_frequency_tigs.tsv",
        REPEATS_annotation_tigs = "{folder}{subfolder_name}5mCartograph_REPEATS_out1/5mCartograph_REPEATS_out1_repeat_annotation_tigs.gff"
    output:
        dir_2 = directory("{folder}{subfolder_name}5mCartograph_REPEATS_out2_tigs/"),
        result_5mer_tigs = "{folder}{subfolder_name}5mCartograph_REPEATS_out2_tigs/5mCartograph_REPEATS_out2_counts_per_5mer_tigs.txt",
        result_chr_tigs = "{folder}{subfolder_name}5mCartograph_REPEATS_out2_tigs/5mCartograph_REPEATS_out2_counts_per_chromosome_tigs.txt"
    shell:
        """
        mkdir -p {output.dir_2}
        python3 /path/to/scripts/5mCartograph_REPEATS/5mCartograph_REPEATS_smk_2_chr_specific_counts.py {input.DSP_tigs} {input.REPEATS_annotation_tigs}
        # python3 /path/to/5mCartograph_repeats/scripts/5mCartograph_REPEATS_smk_2_chr_specific_counts.py {input.DSP_tigs} {input.REPEATS_annotation_tigs}
        """

############################################################################################################################
############################################################################################################################

rule result_summary:
    input:
        result_5mer_chr1 = "{folder}{subfolder_name}5mCartograph_REPEATS_out2_chr1/5mCartograph_REPEATS_out2_counts_per_5mer_chr1.txt",
        result_5mer_chr2 = "{folder}{subfolder_name}5mCartograph_REPEATS_out2_chr2/5mCartograph_REPEATS_out2_counts_per_5mer_chr2.txt",
        result_5mer_chr3 = "{folder}{subfolder_name}5mCartograph_REPEATS_out2_chr3/5mCartograph_REPEATS_out2_counts_per_5mer_chr3.txt",
        result_5mer_chr4 = "{folder}{subfolder_name}5mCartograph_REPEATS_out2_chr4/5mCartograph_REPEATS_out2_counts_per_5mer_chr4.txt",
        result_5mer_chr5 = "{folder}{subfolder_name}5mCartograph_REPEATS_out2_chr5/5mCartograph_REPEATS_out2_counts_per_5mer_chr5.txt",
        result_5mer_chr6 = "{folder}{subfolder_name}5mCartograph_REPEATS_out2_chr6/5mCartograph_REPEATS_out2_counts_per_5mer_chr6.txt",
        result_5mer_chr7 = "{folder}{subfolder_name}5mCartograph_REPEATS_out2_chr7/5mCartograph_REPEATS_out2_counts_per_5mer_chr7.txt",
        result_5mer_chr8 = "{folder}{subfolder_name}5mCartograph_REPEATS_out2_chr8/5mCartograph_REPEATS_out2_counts_per_5mer_chr8.txt",
        result_5mer_chr9 = "{folder}{subfolder_name}5mCartograph_REPEATS_out2_chr9/5mCartograph_REPEATS_out2_counts_per_5mer_chr9.txt",
        result_5mer_tigs = "{folder}{subfolder_name}5mCartograph_REPEATS_out2_tigs/5mCartograph_REPEATS_out2_counts_per_5mer_tigs.txt",
        result_chr_chr1 = "{folder}{subfolder_name}5mCartograph_REPEATS_out2_chr1/5mCartograph_REPEATS_out2_counts_per_chromosome_chr1.txt",
        result_chr_chr2 = "{folder}{subfolder_name}5mCartograph_REPEATS_out2_chr2/5mCartograph_REPEATS_out2_counts_per_chromosome_chr2.txt",
        result_chr_chr3 = "{folder}{subfolder_name}5mCartograph_REPEATS_out2_chr3/5mCartograph_REPEATS_out2_counts_per_chromosome_chr3.txt",
        result_chr_chr4 = "{folder}{subfolder_name}5mCartograph_REPEATS_out2_chr4/5mCartograph_REPEATS_out2_counts_per_chromosome_chr4.txt",
        result_chr_chr5 = "{folder}{subfolder_name}5mCartograph_REPEATS_out2_chr5/5mCartograph_REPEATS_out2_counts_per_chromosome_chr5.txt",
        result_chr_chr6 = "{folder}{subfolder_name}5mCartograph_REPEATS_out2_chr6/5mCartograph_REPEATS_out2_counts_per_chromosome_chr6.txt",
        result_chr_chr7 = "{folder}{subfolder_name}5mCartograph_REPEATS_out2_chr7/5mCartograph_REPEATS_out2_counts_per_chromosome_chr7.txt",
        result_chr_chr8 = "{folder}{subfolder_name}5mCartograph_REPEATS_out2_chr8/5mCartograph_REPEATS_out2_counts_per_chromosome_chr8.txt",
        result_chr_chr9 = "{folder}{subfolder_name}5mCartograph_REPEATS_out2_chr9/5mCartograph_REPEATS_out2_counts_per_chromosome_chr9.txt",
        result_chr_tigs = "{folder}{subfolder_name}5mCartograph_REPEATS_out2_tigs/5mCartograph_REPEATS_out2_counts_per_chromosome_tigs.txt"
    output:
        dir_3 = directory("{folder}{subfolder_name}5mCartograph_REPEATS_out3/"),
        result_all = "{folder}{subfolder_name}5mCartograph_REPEATS_out3/5mCartograph_REPEATS_smk_out3_log_and_results.txt"        
    shell:
        """
        mkdir -p {output.dir_3}
        python3 /path/to/scripts/5mCartograph_REPEATS/5mCartograph_REPEATS_smk_3_summarize_calculate_and_output.py {input.result_5mer_chr1} {input.result_5mer_chr2} {input.result_5mer_chr3} {input.result_5mer_chr4} {input.result_5mer_chr5} {input.result_5mer_chr6} {input.result_5mer_chr7} {input.result_5mer_chr8} {input.result_5mer_chr9} {input.result_5mer_tigs} {input.result_chr_chr1} {input.result_chr_chr2} {input.result_chr_chr3} {input.result_chr_chr4} {input.result_chr_chr5} {input.result_chr_chr6} {input.result_chr_chr7} {input.result_chr_chr8} {input.result_chr_chr9} {input.result_chr_tigs}
        # python3 /path/to/5mCartograph_repeats/scripts/5mCartograph_REPEATS_smk_3_summarize_calculate_and_output.py {input.result_5mer_chr1} {input.result_5mer_chr2} {input.result_5mer_chr3} {input.result_5mer_chr4} {input.result_5mer_chr5} {input.result_5mer_chr6} {input.result_5mer_chr7} {input.result_5mer_chr8} {input.result_5mer_chr9} {input.result_5mer_tigs} {input.result_chr_chr1} {input.result_chr_chr2} {input.result_chr_chr3} {input.result_chr_chr4} {input.result_chr_chr5} {input.result_chr_chr6} {input.result_chr_chr7} {input.result_chr_chr8} {input.result_chr_chr9} {input.result_chr_tigs}
        """

############################################################################################################################
############################################################################################################################

rule plot:
    input:
        plot_data = "{folder}{subfolder_name}5mCartograph_REPEATS_out3/5mCartograph_REPEATS_smk_out3_log_and_results.txt"
    output:
        png = "{folder}{subfolder_name}5mCartograph_REPEATS_out3/5mCartograph_plot_results_grouped_repeat_gry.png",
        log = "{folder}{subfolder_name}5mCartograph_REPEATS_out3/5mCartograph_plot_results_grouped_repeat_gry_log.txt"
    params:
        title = expand("{title}", title=plot_title)
    shell:
        """
        python3 /path/to/scripts/5mCartograph_REPEATS/5mCartograph_REPEATS_smk_4_plot_results_grouped_repeats_gry.py {input.plot_data} {params.title}
        # python3 /path/to/5mCartograph_repeats/scripts/5mCartograph_REPEATS_smk_4_plot_results_grouped_repeats_gry.py {input.plot_data} {params.title}
        """

############################################################################################################################
############################################################################################################################

rule all:
    input:
        # file that shall be created with the workflow
        result_all = expand("{folder}{subfolder_name}5mCartograph_REPEATS_out3/5mCartograph_REPEATS_smk_out3_log_and_results.txt", folder=folder_name, subfolder_name=subfolder_name),
        png = expand("{folder}{subfolder_name}5mCartograph_REPEATS_out3/5mCartograph_plot_results_grouped_repeat_gry.png", folder=folder_name, subfolder_name=subfolder_name),
        log = expand("{folder}{subfolder_name}5mCartograph_REPEATS_out3/5mCartograph_plot_results_grouped_repeat_gry_log.txt", folder=folder_name, subfolder_name=subfolder_name)
    run:
        print(datetime.now() - startTime)

