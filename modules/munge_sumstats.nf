#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process munge_sumstats {
    label 'low_mem'
    input:
        tuple path(sumstats),
            val(n),
            val(info),
            val(maf),
            val(out_prefix),
            path(munge_sumstats_script)
    output:
        path("${out_prefix}_munged.parquet")
    script:
        """
        python munge_sumstats.py --sumstats ${sumstats} \
            --n ${n} \
            --min-info ${info} \
            --min-maf ${maf} \
            --out ${out_prefix}_munged.parquet
        """
}