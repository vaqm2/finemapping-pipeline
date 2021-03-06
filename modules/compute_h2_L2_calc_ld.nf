#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process compute_h2_L2_calc_ld {
    label 'big_mem'

    input:
        tuple val(chr),
        val(annotation_files_prefix),
        path(annotations_file),
        path(ld_file),
        path(metrics_file),
        val(weight_files_prefix),
        path(weights_file),
        val(bfile_prefix),
        path(bed),
        path(bim),
        path(fam),
        path(polyfun_script),
        val(out),
        path(munged_sumstats)

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
            --ref-ld-chr ${annotation_files_prefix}. \
            --w-ld-chr ${weight_files_prefix}. \
            --ld-wind-kb 1000 \
            --bfile ${bfile_prefix}${chr} \
            --allow-missing \
            --chr $chr
        """
}