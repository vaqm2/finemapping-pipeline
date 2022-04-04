#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process reestimate_snp_h2 {
    label 'mod_mem'
    publishDir launchDir
    input:
        path polyfun_script,
        val out,
        path munged_sumstats,
        val weight_file_prefix,
        path weight_files,
        path step2_outputs
    output:
        path "${out}.*.snpvar_constrained.gz"
    script:
        """
        python ./polyfun.py \
            --compute-h2-bins \
            --output-prefix $out \
            --sumstats $munged_sumstats \
            --w-ld-chr ${weight_file_prefix}.
        """
}