import subprocess

workdir: "."

INPUT_BED = "initial/gwas_data.bed"
INPUT_BIM = "initial/gwas_data.bim"
INPUT_FAM = "initial/gwas_data.fam"


GROUPS = ["group_1", "group_2", "group_3", "group_4", "group_5"]
THRESHOLDS = ["0.01", "0.05", "0.1"]

INITIAL_DIR = "initial"
QC_DIR = "qc"
FINAL_DIR = "qc_final"
PCA_DIR = "pca"
PLOT_DIR = "plots"
GWAS_DIR = "gwas"
PGS_DIR = "pgs"

rule all:
    input:
        f"{FINAL_DIR}/QC_output_final.bed",
        f"{FINAL_DIR}/QC_output_final.bim",
        f"{FINAL_DIR}/QC_output_final.fam",
        f"{PLOT_DIR}/qc_visualization.png",
        f"{PLOT_DIR}/pca_plot.png",
        f"{PLOT_DIR}/linear_gwaslab_mqq.png",
        f"{PLOT_DIR}/linear_gwaslab_manhattan.png",
        f"{PLOT_DIR}/mlma_gwaslab_mqq.png",
        f"{PLOT_DIR}/mlma_gwaslab_manhattan.png",
        f"{PLOT_DIR}/pgs_evaluation.png"

#
# -----QC-----
#

# I. Sample quality control
# 1. sex check
rule sex_check:
    input:
        bed = INPUT_BED,
        bim = INPUT_BIM,
        fam = INPUT_FAM
    output:
        sexcheck = f"{QC_DIR}/QC_output_sexcheck.sexcheck"
    shell:
        """
        plink --bfile initial/gwas_data --check-sex --allow-no-sex --out {QC_DIR}/QC_output_sexcheck
        """

rule identify_sex_problems:
    input:
        sexcheck = f"{QC_DIR}/QC_output_sexcheck.sexcheck"
    output:
        fails = f"{QC_DIR}/fail-sexcheck-qc.txt"
    shell:
        """
        awk '$4 == 0 || $5 == "PROBLEM" {{print $1, $2}}' {input.sexcheck} > {output.fails}
        """

rule filter_sex_problems:
    input:
        bed = INPUT_BED,
        bim = INPUT_BIM,
        fam = INPUT_FAM,
        fails = f"{QC_DIR}/fail-sexcheck-qc.txt"
    output:
        bed = f"{QC_DIR}/QC_output_sex_filtered.bed",
        bim = f"{QC_DIR}/QC_output_sex_filtered.bim",
        fam = f"{QC_DIR}/QC_output_sex_filtered.fam"
    shell:
        """
        plink --bfile initial/gwas_data --remove {input.fails} --make-bed --out {QC_DIR}/QC_output_sex_filtered
        """

# 2. missingness analysis
rule missingness_analysis:
    input:
        bed = f"{QC_DIR}/QC_output_sex_filtered.bed",
        bim = f"{QC_DIR}/QC_output_sex_filtered.bim",
        fam = f"{QC_DIR}/QC_output_sex_filtered.fam"
    output:
        imiss = f"{QC_DIR}/QC_output_missingness.imiss",
        lmiss = f"{QC_DIR}/QC_output_missingness.lmiss"
    shell:
        """
        plink --bfile {QC_DIR}/QC_output_sex_filtered --missing --out {QC_DIR}/QC_output_missingness
        """

# 3. heterozygosity analysis
rule heterozygosity_analysis:
    input:
        bed = f"{QC_DIR}/QC_output_sex_filtered.bed",
        bim = f"{QC_DIR}/QC_output_sex_filtered.bim",
        fam = f"{QC_DIR}/QC_output_sex_filtered.fam"
    output:
        het = f"{QC_DIR}/QC_output_het.het"
    shell:
        """
        plink --bfile {QC_DIR}/QC_output_sex_filtered --het --out {QC_DIR}/QC_output_het
        """

rule analyze_het_in_r:
    input:
        het = f"{QC_DIR}/QC_output_het.het",
        imiss = f"{QC_DIR}/QC_output_missingness.imiss"
    output:
        outliers = f"{QC_DIR}/wrong_het_missing_values.txt"
    script:
        "scripts/analyze_heterozygosity.R"

rule filter_heterozygosity:
    input:
        bed = f"{QC_DIR}/QC_output_sex_filtered.bed",
        bim = f"{QC_DIR}/QC_output_sex_filtered.bim",
        fam = f"{QC_DIR}/QC_output_sex_filtered.fam",
        outliers = f"{QC_DIR}/wrong_het_missing_values.txt"
    output:
        bed = f"{QC_DIR}/QC_output_het_filtered.bed",
        bim = f"{QC_DIR}/QC_output_het_filtered.bim",
        fam = f"{QC_DIR}/QC_output_het_filtered.fam"
    shell:
        """
        plink --bfile {QC_DIR}/QC_output_sex_filtered --remove {input.outliers} --make-bed --out {QC_DIR}/QC_output_het_filtered
        """

# 4. IBD check
rule pruning_for_ibd:
    input:
        bed = f"{QC_DIR}/QC_output_het_filtered.bed",
        bim = f"{QC_DIR}/QC_output_het_filtered.bim",
        fam = f"{QC_DIR}/QC_output_het_filtered.fam"
    output:
        prune_in = f"{QC_DIR}/QC_output_pruning.prune.in",
        prune_out = f"{QC_DIR}/QC_output_pruning.prune.out"
    shell:
        """
        plink --bfile {QC_DIR}/QC_output_het_filtered --indep-pairwise 500kb 5 0.2 --out {QC_DIR}/QC_output_pruning
        """

rule calculate_ibd:
    input:
        bed = f"{QC_DIR}/QC_output_het_filtered.bed",
        bim = f"{QC_DIR}/QC_output_het_filtered.bim",
        fam = f"{QC_DIR}/QC_output_het_filtered.fam",
        prune_in = f"{QC_DIR}/QC_output_pruning.prune.in"
    output:
        genome = f"{QC_DIR}/QC_output_ibd.genome"
    shell:
        """
        plink --bfile {QC_DIR}/QC_output_het_filtered --extract {input.prune_in} --genome --min 0.185 --out {QC_DIR}/QC_output_ibd
        """

rule analyze_ibd_in_r:
    input:
        genome = f"{QC_DIR}/QC_output_ibd.genome"
    output:
        ibd_fails = f"{QC_DIR}/wrong_ibd.txt"
    script:
        "scripts/analyze_relatedness.R"

rule filter_relatedness:
    input:
        bed = f"{QC_DIR}/QC_output_het_filtered.bed",
        bim = f"{QC_DIR}/QC_output_het_filtered.bim",
        fam = f"{QC_DIR}/QC_output_het_filtered.fam",
        ibd_fails = f"{QC_DIR}/wrong_ibd.txt"
    output:
        bed = f"{QC_DIR}/QC_output_sample_final.bed",
        bim = f"{QC_DIR}/QC_output_sample_final.bim",
        fam = f"{QC_DIR}/QC_output_sample_final.fam"
    shell:
        """
        plink --bfile {QC_DIR}/QC_output_het_filtered --remove {input.ibd_fails} --make-bed --out {QC_DIR}/QC_output_sample_final
        """

# II. SNP quality control
# SNP quality control by chips
rule create_group_files:
    input:
        metadata = f"{INITIAL_DIR}/metadata.txt"
    output:
        group_files = expand(f"{INITIAL_DIR}/group_{{n}}.txt", n=[1, 2, 3, 4, 5])
    script:
        "scripts/divide_dataset_by_chip.R"

rule divide_by_group:
    input:
        bed = f"{QC_DIR}/QC_output_sample_final.bed",
        bim = f"{QC_DIR}/QC_output_sample_final.bim",
        fam = f"{QC_DIR}/QC_output_sample_final.fam",
        group_file = f"{INITIAL_DIR}/{{group}}.txt"
    output:
        bed = f"{QC_DIR}/QC_output_{{group}}.bed",
        bim = f"{QC_DIR}/QC_output_{{group}}.bim",
        fam = f"{QC_DIR}/QC_output_{{group}}.fam"
    wildcard_constraints:
        group = "|".join(GROUPS)
    shell:
        """
        plink --bfile {QC_DIR}/QC_output_sample_final --keep {input.group_file} --make-bed --out {QC_DIR}/QC_output_{wildcards.group}
        """

rule filter_snps_by_group:
    input:
        bed = f"{QC_DIR}/QC_output_{{group}}.bed",
        bim = f"{QC_DIR}/QC_output_{{group}}.bim",
        fam = f"{QC_DIR}/QC_output_{{group}}.fam"
    output:
        bed = f"{QC_DIR}/QC_output_geno_filtered_{{group}}.bed",
        bim = f"{QC_DIR}/QC_output_geno_filtered_{{group}}.bim",
        fam = f"{QC_DIR}/QC_output_geno_filtered_{{group}}.fam"
    wildcard_constraints:
        group = "|".join(GROUPS)
    shell:
        """
        plink --bfile {QC_DIR}/QC_output_{wildcards.group} --geno 0.05 --hwe 0.00001 --maf 0.01 --make-bed --out {QC_DIR}/QC_output_geno_filtered_{wildcards.group}
        """

# merge filtered groups
rule create_merge_list:
    input:
        expand(f"{QC_DIR}/QC_output_geno_filtered_{{group}}.bed", group=GROUPS)
    output:
        merge_list = f"{QC_DIR}/merge_list.txt"
    run:
        with open(output.merge_list, 'w') as f:
            for group in GROUPS:
                fam_file = f"{QC_DIR}/QC_output_geno_filtered_{group}.fam"
                try:
                    result = subprocess.run(
                        ["wc", "-l", fam_file],
                        capture_output=True,
                        text=True,
                        check=True
                    )
                    lines = int(result.stdout.strip().split()[0])
                    if lines > 50:
                        f.write(f"{QC_DIR}/QC_output_geno_filtered_{group}\n")
                    else:
                        print(f"Skipping {fam_file} with length of {lines}")
                except Exception:
                    print(f"Error checking {fam_file} length")

rule merge_filtered_groups:
    input:
        merge_list = f"{QC_DIR}/merge_list.txt",
        bedfiles = expand(f"{QC_DIR}/QC_output_geno_filtered_{{group}}.bed", group=GROUPS),
        bimfiles = expand(f"{QC_DIR}/QC_output_geno_filtered_{{group}}.bim", group=GROUPS),
        famfiles = expand(f"{QC_DIR}/QC_output_geno_filtered_{{group}}.fam", group=GROUPS)
    output:
        bed = f"{FINAL_DIR}/QC_output_final.bed",
        bim = f"{FINAL_DIR}/QC_output_final.bim",
        fam = f"{FINAL_DIR}/QC_output_final.fam"
    shell:
        """
        plink --merge-list {input.merge_list} --make-bed --out {FINAL_DIR}/QC_output_final
        """

rule create_qc_visualizations:
    input:
        final_bed = f"{FINAL_DIR}/QC_output_final.bed",
        final_bim = f"{FINAL_DIR}/QC_output_final.bim",
        final_fam = f"{FINAL_DIR}/QC_output_final.fam"
    output:
        plots = f"{PLOT_DIR}/qc_visualization.png"
    script:
        "scripts/visualize_qc.R"

#
# -----PCA-----
#

rule prune_snps_for_pca:
    input:
        bed = f"{FINAL_DIR}/QC_output_final.bed",
        bim = f"{FINAL_DIR}/QC_output_final.bim",
        fam = f"{FINAL_DIR}/QC_output_final.fam"
    output:
        prune_in = f"{PCA_DIR}/gwas.prune.in",
        prune_out = f"{PCA_DIR}/gwas.prune.out"
    shell:
        """
        plink --bfile {FINAL_DIR}/QC_output_final --indep-pairwise 100kb 5 0.2 --out {PCA_DIR}/gwas
        """

rule run_pca:
    input:
        bed = f"{FINAL_DIR}/QC_output_final.bed",
        bim = f"{FINAL_DIR}/QC_output_final.bim",
        fam = f"{FINAL_DIR}/QC_output_final.fam",
        prune_in = f"{PCA_DIR}/gwas.prune.in"
    output:
        eigenvec = f"{PCA_DIR}/gwas.eigenvec",
        eigenval = f"{PCA_DIR}/gwas.eigenval"
    shell:
        """
        plink --bfile {FINAL_DIR}/QC_output_final --extract {input.prune_in} --pca 20 --out {PCA_DIR}/gwas
        """

rule plot_pca:
    input:
        eigenvec = f"{PCA_DIR}/gwas.eigenvec",
        eigenval = f"{PCA_DIR}/gwas.eigenval"
    output:
        pca_plot = f"{PLOT_DIR}/pca_plot.png",
        scree_plot = f"{PLOT_DIR}/scree_plot.png"
    script:
        "scripts/plot_pca.R"

#
# -----GWAS-----
#

rule prepare_covariates:
    input:
        eigenvec = f"{PCA_DIR}/gwas.eigenvec"
    output:
        covariates = f"{GWAS_DIR}/pca_covariates.txt"
    shell:
        """
        awk '{{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12}}' {input.eigenvec} > {output.covariates}
        """

rule prepare_phenotype:
    input:
        height = f"{INITIAL_DIR}/height.txt"
    output:
        formatted_height = f"{INITIAL_DIR}/right_format_height.txt"
    shell:
        """
        awk '{{print $1, $1, $2}}' {input.height} > {output.formatted_height}
        """

rule run_gwas_linear:
    input:
        bed = f"{FINAL_DIR}/QC_output_final.bed",
        bim = f"{FINAL_DIR}/QC_output_final.bim",
        fam = f"{FINAL_DIR}/QC_output_final.fam",
        phenotype = f"{INITIAL_DIR}/right_format_height.txt",
        covariates = f"{GWAS_DIR}/pca_covariates.txt"
    output:
        gwas_results = f"{GWAS_DIR}/gwas_height.assoc.linear"
    shell:
        """
        plink --bfile {FINAL_DIR}/QC_output_final --pheno {input.phenotype} --covar {input.covariates} --linear --allow-no-sex --out {GWAS_DIR}/gwas_height
        """

rule extract_add_results:
    input:
        linear_results = f"{GWAS_DIR}/gwas_height.assoc.linear"
    output:
        add_results = f"{GWAS_DIR}/gwas_height.ADD.assoc.linear"
    shell:
        """
        awk 'NR==1 || $5=="ADD"' {input.linear_results} > {output.add_results}
        """

rule calculate_grm:
    input:
        bed = f"{FINAL_DIR}/QC_output_final.bed",
        bim = f"{FINAL_DIR}/QC_output_final.bim",
        fam = f"{FINAL_DIR}/QC_output_final.fam"
    output:
        grm_bin = f"{GWAS_DIR}/gwas_grm.grm.bin",
        grm_n = f"{GWAS_DIR}/gwas_grm.grm.N.bin",
        grm_id = f"{GWAS_DIR}/gwas_grm.grm.id"
    shell:
        """
        plink --bfile {FINAL_DIR}/QC_output_final --make-grm-bin --out {GWAS_DIR}/gwas_grm
        """

rule run_gcta_mlma:
    input:
        bed = f"{FINAL_DIR}/QC_output_final.bed",
        bim = f"{FINAL_DIR}/QC_output_final.bim",
        fam = f"{FINAL_DIR}/QC_output_final.fam",
        grm_bin = f"{GWAS_DIR}/gwas_grm.grm.bin",
        grm_n = f"{GWAS_DIR}/gwas_grm.grm.N.bin",
        grm_id = f"{GWAS_DIR}/gwas_grm.grm.id",
        phenotype = f"{INITIAL_DIR}/right_format_height.txt",
        covariates = f"{GWAS_DIR}/pca_covariates.txt"
    output:
        mlma_results = f"{GWAS_DIR}/gwas_height_gcta_mlma.mlma"
    shell:
        """
        gcta64 --mlma --bfile {FINAL_DIR}/QC_output_final --grm {GWAS_DIR}/gwas_grm --pheno {input.phenotype} --qcovar {input.covariates} --thread-num 4 --out {GWAS_DIR}/gwas_height_gcta_mlma
        """

rule prepare_gwas_data:
    input:
        gwas_results = f"{GWAS_DIR}/gwas_height.assoc.linear",
        mlma_results = f"{GWAS_DIR}/gwas_height_gcta_mlma.mlma"
    output:
        linear_clean = f"{GWAS_DIR}/linear_clean.tsv",
        mlma_clean = f"{GWAS_DIR}/mlma_clean.tsv"
    script:
        "scripts/prepare_gwas_data.R"

rule run_gwaslab_plots:
    input:
        script = "scripts/run_gwaslab.py",
        linear_clean = f"{GWAS_DIR}/linear_clean.tsv",
        mlma_clean = f"{GWAS_DIR}/mlma_clean.tsv"
    output:
        linear_mqq = f"{PLOT_DIR}/linear_gwaslab_mqq.png",
        linear_manhattan = f"{PLOT_DIR}/linear_gwaslab_manhattan.png",
        mlma_mqq = f"{PLOT_DIR}/mlma_gwaslab_mqq.png",
        mlma_manhattan = f"{PLOT_DIR}/mlma_gwaslab_manhattan.png"
    conda:
        "gwaslab_environment.yaml"
    shell:
        """
        python {input.script}
        """
        
#
# -----PGS PREDICTION-----
#

rule split_data_for_pgs:
    input:
        height = f"{INITIAL_DIR}/right_format_height.txt",
        covariates = f"{GWAS_DIR}/pca_covariates.txt"
    output:
        train_samples = f"{PGS_DIR}/train_samples.txt",
        test_samples = f"{PGS_DIR}/test_samples.txt",
        train_heights = f"{PGS_DIR}/train_heights.txt",
        test_heights = f"{PGS_DIR}/test_heights.txt",
        train_covariates = f"{PGS_DIR}/train_covariates.txt",
        test_covariates = f"{PGS_DIR}/test_covariates.txt"
    script:
        "scripts/split_data.R"

rule create_train_dataset:
    input:
        bed = f"{FINAL_DIR}/QC_output_final.bed",
        bim = f"{FINAL_DIR}/QC_output_final.bim",
        fam = f"{FINAL_DIR}/QC_output_final.fam",
        train_samples = f"{PGS_DIR}/train_samples.txt"
    output:
        bed = f"{PGS_DIR}/pgs_train_data.bed",
        bim = f"{PGS_DIR}/pgs_train_data.bim",
        fam = f"{PGS_DIR}/pgs_train_data.fam"
    shell:
        """
        plink --bfile {FINAL_DIR}/QC_output_final --keep {input.train_samples} --make-bed --out {PGS_DIR}/pgs_train_data
        """
        
rule create_test_dataset:
    input:
        bed = f"{FINAL_DIR}/QC_output_final.bed",
        bim = f"{FINAL_DIR}/QC_output_final.bim",
        fam = f"{FINAL_DIR}/QC_output_final.fam",
        test_samples = f"{PGS_DIR}/test_samples.txt"
    output:
        bed = f"{PGS_DIR}/pgs_test_data.bed",
        bim = f"{PGS_DIR}/pgs_test_data.bim",
        fam = f"{PGS_DIR}/pgs_test_data.fam"
    shell:
        """
        plink --bfile {FINAL_DIR}/QC_output_final --keep {input.test_samples} --make-bed --out {PGS_DIR}/pgs_test_data
        """
        
rule gwas_on_training_data:
    input:
        bed = f"{PGS_DIR}/pgs_train_data.bed",
        bim = f"{PGS_DIR}/pgs_train_data.bim",
        fam = f"{PGS_DIR}/pgs_train_data.fam",
        phenotype = f"{PGS_DIR}/train_heights.txt",
        covariates = f"{PGS_DIR}/train_covariates.txt"
    output:
        gwas_results = f"{PGS_DIR}/gwas_train.assoc.linear"
    shell:
        """
        plink --bfile {PGS_DIR}/pgs_train_data --pheno {input.phenotype} --covar {input.covariates} --linear --allow-no-sex --out {PGS_DIR}/gwas_train
        """
        
rule extract_add_results_train:
    input:
        linear_results = f"{PGS_DIR}/gwas_train.assoc.linear"
    output:
        add_results = f"{PGS_DIR}/gwas_train.ADD.assoc.linear"
    shell:
        """
        awk 'NR==1 || $5=="ADD"' {input.linear_results} > {output.add_results}
        """
        
rule ld_clumping:
    input:
        bed = f"{PGS_DIR}/pgs_train_data.bed",
        bim = f"{PGS_DIR}/pgs_train_data.bim",
        fam = f"{PGS_DIR}/pgs_train_data.fam",
        gwas_results = f"{PGS_DIR}/gwas_train.ADD.assoc.linear"
    output:
        clumped = f"{PGS_DIR}/pgs_clumped.clumped"
    shell:
        """
        plink --bfile {PGS_DIR}/pgs_train_data --clump {input.gwas_results} --clump-p1 0.0001 --clump-p2 0.01 --clump-r2 0.50 --clump-kb 500 --out {PGS_DIR}/pgs_clumped
        """
        
rule extract_clumped_snps:
    input:
        clumped = f"{PGS_DIR}/pgs_clumped.clumped"
    output:
        snp_list = f"{PGS_DIR}/pgs_clumped_snps.txt"
    shell:
        """
        awk 'NR>1 {{print $3}}' {input.clumped} > {output.snp_list}
        """
        
rule create_score_files:
    input:
        gwas_results = f"{PGS_DIR}/gwas_train.ADD.assoc.linear"
    output:
        score_file = f"{PGS_DIR}/score_file_p{{threshold}}.txt"
    params:
        threshold = lambda wildcards: wildcards.threshold
    wildcard_constraints:
        threshold = "|".join(THRESHOLDS)
    shell:
        """
        awk '$9 < {params.threshold} {{print $2, $4, $7}}' {input.gwas_results} > {output.score_file}
        """
        
rule calculate_pgs:
    input:
        bed = f"{PGS_DIR}/pgs_test_data.bed",
        bim = f"{PGS_DIR}/pgs_test_data.bim",
        fam = f"{PGS_DIR}/pgs_test_data.fam",
        score_file = f"{PGS_DIR}/score_file_p{{threshold}}.txt"
    output:
        profile = f"{PGS_DIR}/pgs_results_p{{threshold}}.profile"
    wildcard_constraints:
        threshold = "|".join(THRESHOLDS)
    shell:
        """
        plink --bfile {PGS_DIR}/pgs_test_data --score {input.score_file} 1 2 3 --out {PGS_DIR}/pgs_results_p{wildcards.threshold}
        """
        
rule evaluate_pgs:
    input:
        pgs_profiles = expand(f"{PGS_DIR}/pgs_results_p{{threshold}}.profile", threshold=THRESHOLDS),
        test_heights = f"{PGS_DIR}/test_heights.txt",
        test_covariates = f"{PGS_DIR}/test_covariates.txt"
    output:
        evaluation_plot = f"{PLOT_DIR}/pgs_evaluation.png"
    params:
        thresholds = THRESHOLDS
    script:
        "scripts/evaluate_pgs.R"
