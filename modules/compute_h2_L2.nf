#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process compute_h2_L2 {
    input:
        tuple val(chr),
            val(base_prefix),
            path(ld),
            path(metrics),
            path(sumstats),
            val(out_prefix),
            path(annotations),
            path(weights)
    output:
        tuple path("${out_prefix}.annot_coeff_ridge.*.txt"),
            path("${out_prefix}.*.snpvar_ridge_constrained.*.gz"),
            path("${out_prefix}.*.snpvar_ridge.*.gz"), 
            path("${out_prefix}.*.l2.M"),
            path("${out_prefix}.*.bins.parquet")

    script:
        """
        python polyfun.py \
            --compute-h2-L2 \
            --output-prefix $out_prefix \
            --sumtats $sumstats \
            --ref-ld-chr $annotations \
            --w-ld-chr $weights \
            --ld-wind-kb 1000 \
            --allow-missing
        """
}