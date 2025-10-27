# snakefile for the analysis of 5mC in the whole genome of sugar beet
# the original basic_5mC_analysis python script was split into three subscripts

# FILES NEEDED:
# 5mCartograph.smk -> parallelization of subscript 2 and coordination of subscripts
# Config_5mCargotraph_....yaml -> configuration of sample specific files (individual for each sample)
# 5mCartograph_smk_1_read_and_split_DSP_files.py -> split and save the DSP output and annotation file
# 5mCartograph_smk_2_chr_specific_counts.py -> calculate 5mer and chromosome specific data (the chromosomes are processed individually, while the contigs run in the 10th slot)
# 5mCartograph_smk_3_summarize_calculate_and_output.py -> summarizes data of subsets and calculates data of interest
# 5mCartograph_smk_4_plot_results.py -> plots results with matplotlib

# tsv needs to be placed in destination directory
# snakefile needs to be used in 'snakemake' environment 
# Command: snakemake all --cores 10 --snakefile 5mCartograph.smk --configfile path/to/config.yaml

############################################################################################################################
############################################################################################################################

# harbours specific information about the data sets used
# configfile: '/path/to/config_file.py'

# variables
tool_folder_name = config['tool_folder_name']
result_folder_name = config['result_folder_name']
tsv_path = config['tsv_path']
anno_path = config['anno_path']
plot_title_1 = config['plot_title_1']
plot_title_2 = config['plot_title_2']

############################################################################################################################
############################################################################################################################

checkpoint split:
    input:
        DSP = expand("{tsv}", tsv=tsv_path),
        annotation = expand("{anno}", anno=anno_path)
    output:
        dir_1 = directory("{folder}5mCartograph_out1/")
    shell:
        """
        mkdir -p {output.dir_1}
        python3 {tool_folder_name}5mCartograph_smk_1_read_and_split_DSP_files.py {input.DSP} {input.annotation} {result_folder_name}
        """

############################################################################################################################
############################################################################################################################

def get_all_pairs(wildcards):
    ckpt = checkpoints.split.get(**wildcards)
    import glob, os
    tsv_files = sorted(glob.glob(os.path.join(ckpt.output[0], "5mCartograph_out1_frequency_*.tsv")))
    gff_files = sorted(glob.glob(os.path.join(ckpt.output[0], "5mCartograph_out1_annotation_*.gff")))

    pairs = []
    names = []
    for ftsv in tsv_files:
        idx = os.path.basename(ftsv).split("_")[3].split(".")[0]
        fgff = os.path.join(ckpt.output[0], f"5mCartograph_out1_annotation_{idx}.gff")
        if fgff in gff_files:
            pairs.append((ftsv, fgff))
            names.append(idx)
    return pairs, names

rule calculate_chr:
    input:
        lambda wildcards: get_all_pairs(wildcards)[0][get_all_pairs(wildcards)[1].index(wildcards.i)][0],
        lambda wildcards: get_all_pairs(wildcards)[0][get_all_pairs(wildcards)[1].index(wildcards.i)][1]
    output:
        dir_2 = directory("{folder}5mCartograph_out2_{i}/"),
        result_5mer = "{folder}5mCartograph_out2_{i}/5mCartograph_out2_counts_per_5mer_{i}.txt",
        result_chr = "{folder}5mCartograph_out2_{i}/5mCartograph_out2_counts_per_chromosome_{i}.txt"
    shell:
        """
        mkdir -p {output.dir_2}
        python3 {tool_folder_name}5mCartograph_smk_2_chr_specific_counts.py {input[0]} {input[1]} {result_folder_name}
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
        idy = os.path.basename(ftsv).split("_")[3].split(".")[0]
        outs_1.append(f"{result_folder_name}5mCartograph_out2_" + idy + "/5mCartograph_out2_counts_per_5mer_" + idy + ".txt")
        outs_2.append(f"{result_folder_name}5mCartograph_out2_" + idy + "/5mCartograph_out2_counts_per_chromosome_" + idy + ".txt")
    # return lists
    return outs_1 + outs_2

rule result_summary:
    input:
        get_subset_result_files
    output:
        dir_3 = directory("{folder}5mCartograph_out3/"),
        result_all = "{folder}5mCartograph_out3/5mCartograph_out3_log_and_results.txt"        
    shell:
        """
        mkdir -p {output.dir_3}
        python3 {tool_folder_name}5mCartograph_smk_3_summarize_calculate_and_output.py {input} {result_folder_name}
        """

############################################################################################################################
############################################################################################################################

rule plot_genome_wide:
    input:
        plot_data = "{folder}5mCartograph_out3/5mCartograph_out3_log_and_results.txt"
    output:
        png = "{folder}5mCartograph_out3/5mCartograph_plot_results_grouped_genome_wide_gry_bar_hight.png",
        log = "{folder}5mCartograph_out3/5mCartograph_plot_results_grouped_genome_wide_gry_log.txt"
    params:
        title = expand("{title}", title=plot_title_1)
    shell:
        """
        python3 {tool_folder_name}5mCartograph_smk_4_plot_results_grouped_genome_wide_gry.py {input.plot_data} {params.title}
        """

############################################################################################################################
############################################################################################################################

rule plot_all_genes:
    input:
        plot_data = "{folder}5mCartograph_out3/5mCartograph_out3_log_and_results.txt"
    output:
        png = "{folder}5mCartograph_out3/5mCartograph_plot_results_grouped_all_genes_gry.png",
        log = "{folder}5mCartograph_out3/5mCartograph_plot_results_grouped_all_genes_gry_log.txt"
    params:
        title = expand("{title}", title=plot_title_2)
    shell:
        """
        python3 {tool_folder_name}5mCartograph_smk_4_plot_results_grouped_all_genes_gry.py {input.plot_data} {params.title}
        """

############################################################################################################################
############################################################################################################################

rule all:
    input:
        # file that shall be created with the workflow
        result_all = expand("{folder}5mCartograph_out3/5mCartograph_out3_log_and_results.txt", folder=result_folder_name),
        png_1 = expand("{folder}5mCartograph_out3/5mCartograph_plot_results_grouped_genome_wide_gry_bar_hight.png", folder=result_folder_name),
        log_1 = expand("{folder}5mCartograph_out3/5mCartograph_plot_results_grouped_genome_wide_gry_log.txt", folder=result_folder_name),
        png_2 = expand("{folder}5mCartograph_out3/5mCartograph_plot_results_grouped_all_genes_gry.png", folder=result_folder_name),
        log_2 = expand("{folder}5mCartograph_out3/5mCartograph_plot_results_grouped_all_genes_gry_log.txt", folder=result_folder_name)

