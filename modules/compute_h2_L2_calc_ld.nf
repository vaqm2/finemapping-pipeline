#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process compute_h2_L2_calc_ld {
    label 'big_mem'
    input:
        tuple val chr,
        val annotation_files_prefix,
        path annotations_file,
        path ld_file,
        path metrics_file,
        val weight_files_prefix,
        path weights_file,
        val ld_dir,
        path polyfun_script,
        val out,
        path munged_sumstats,
        file ld_files
    output:
        tuple path("${out}.annot_coeff_ridge.${chr}.txt"), 
        path("${out}.${chr}.snpvar_ridge_constrained.gz"),
        path("${out}.${chr}.snpvar_ridge.gz"),
        path("${out}.${chr}.l2.M"),
        path("${out}.${chr}.bins.parquet"),
        path("${out}.${chr}.l2.ldscore.parquet")

    script:
        """
        python ./polyfun.py \
            --compute-h2-L2 \
            --compute-ldscores \
            --output-prefix $out \
            --sumtats $munged_sumstats \
            --ref-ld-chr ${annotation_file_prefix}. \
            --w-ld-chr ${weight_file_prefix}. \
            --ld-wind-kb 1000 \
            --ld-ukb \
            --ld-dir $ld_dir \
            --allow-missing \
            --chr $chr
        """
}