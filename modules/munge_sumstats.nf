#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process munge_sumstats {
    label 'low_mem'

    input:
        tuple path(sumstats),
        val(N),
        val(out),
        path(munge_sumstats_script)
    
    output:
        path("${out}_munged.parquet")
    
    script:
        """
        python ./munge_polyfun_sumstats.py \
            --sumstats $sumstats \
            --n $N \
            --out ${out}_munged.parquet
        """
}