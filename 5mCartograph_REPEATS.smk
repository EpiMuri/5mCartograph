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
# Command: snakemake all --cores 10 --snakefile 5mCartograph_REPEATS.smk --configfile path/to/config.yaml

############################################################################################################################
############################################################################################################################

# harbours specific information about the data sets used
# configfile: '/path/to/scripts/config_file'

# variables
tool_folder_name = config['tool_folder_name']
result_folder_name = config['result_folder_name']
subfolder_name = config['subfolder_name']
tsv_path = config['tsv_path']
anno_path = config['anno_path']
plot_title = config['plot_title']

# for efficiency check
startTime = datetime.now()

############################################################################################################################
############################################################################################################################

checkpoint split:
    input:
        DSP = expand("{tsv}", tsv=tsv_path),
        REPEATS_annotation = expand("{anno}", anno=anno_path)   
    output:
        dir_1 = directory("{folder}{subfolder_name}5mCartograph_REPEATS_out1/")
    params:
        subfolder_name = config['subfolder_name']
    shell:
        """
        mkdir -p {output.dir_1}
        python3 {tool_folder_name}5mCartograph_REPEATS_smk_1_read_and_split_DSP_files.py {input.DSP} {input.REPEATS_annotation} {result_folder_name} {params.subfolder_name}
        """

############################################################################################################################
############################################################################################################################

def get_all_pairs(wildcards):
    ckpt = checkpoints.split.get(**wildcards)
    import glob, os
    tsv_files = sorted(glob.glob(os.path.join(ckpt.output[0], "5mCartograph_REPEATS_out1_frequency_*.tsv")))
    gff_files = sorted(glob.glob(os.path.join(ckpt.output[0], "5mCartograph_REPEATS_out1_annotation_*.gff")))

    pairs = []
    names = []
    for ftsv in tsv_files:
        idx = os.path.basename(ftsv).split("_")[4].split(".")[0]
        fgff = os.path.join(ckpt.output[0], f"5mCartograph_REPEATS_out1_annotation_{idx}.gff")
        if fgff in gff_files:
            pairs.append((ftsv, fgff))
            names.append(idx)
    return pairs, names

rule calculate_chr:
    input:
        lambda wildcards: get_all_pairs(wildcards)[0][get_all_pairs(wildcards)[1].index(wildcards.i)][0],
        lambda wildcards: get_all_pairs(wildcards)[0][get_all_pairs(wildcards)[1].index(wildcards.i)][1]
    output:
        dir_2 = directory("{folder}{subfolder_name}5mCartograph_REPEATS_out2_{i}/"),
        result_5mer = "{folder}{subfolder_name}5mCartograph_REPEATS_out2_{i}/5mCartograph_REPEATS_out2_counts_per_5mer_{i}.txt",
        result_chr = "{folder}{subfolder_name}5mCartograph_REPEATS_out2_{i}/5mCartograph_REPEATS_out2_counts_per_chromosome_{i}.txt"
    shell:
        """
        mkdir -p {output.dir_2}
        python3 {tool_folder_name}5mCartograph_REPEATS_smk_2_chr_specific_counts.py {input[0]} {input[1]} {result_folder_name} {subfolder_name}
        """

############################################################################################################################
############################################################################################################################

# get all files generated
def get_subset_result_files(wildcards):
    pairs, names = get_all_pairs(wildcards)
    # print(pairs)
    outs_1 = []
    outs_2 = []
    for ftsv, fgff in pairs:
        idy = os.path.basename(ftsv).split("_")[4].split(".")[0]
        outs_1.append(f"{result_folder_name}{subfolder_name}5mCartograph_REPEATS_out2_" + idy + "/5mCartograph_REPEATS_out2_counts_per_5mer_" + idy + ".txt")
        outs_2.append(f"{result_folder_name}{subfolder_name}5mCartograph_REPEATS_out2_" + idy + "/5mCartograph_REPEATS_out2_counts_per_chromosome_" + idy + ".txt")
    # return lists
    return outs_1 + outs_2


rule result_summary:
    input:
        get_subset_result_files
    output:
        dir_3 = directory("{folder}{subfolder_name}5mCartograph_REPEATS_out3/"),
        result_all = "{folder}{subfolder_name}5mCartograph_REPEATS_out3/5mCartograph_REPEATS_smk_out3_log_and_results.txt"        
    shell:
        """
        mkdir -p {output.dir_3}
        python3 {tool_folder_name}5mCartograph_REPEATS_smk_3_summarize_calculate_and_output.py {input} {result_folder_name} {subfolder_name}
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
        python3 {tool_folder_name}5mCartograph_REPEATS_smk_4_plot_results_grouped_repeats_gry.py {input.plot_data} {params.title}
        """

############################################################################################################################
############################################################################################################################

rule all:
    input:
        # file that shall be created with the workflow
        result_all = expand("{folder}{subfolder_name}5mCartograph_REPEATS_out3/5mCartograph_REPEATS_smk_out3_log_and_results.txt", folder=result_folder_name, subfolder_name=subfolder_name),
        png = expand("{folder}{subfolder_name}5mCartograph_REPEATS_out3/5mCartograph_plot_results_grouped_repeat_gry.png", folder=result_folder_name, subfolder_name=subfolder_name),
        log = expand("{folder}{subfolder_name}5mCartograph_REPEATS_out3/5mCartograph_plot_results_grouped_repeat_gry_log.txt", folder=result_folder_name, subfolder_name=subfolder_name)
    run:
        print(datetime.now() - startTime)

