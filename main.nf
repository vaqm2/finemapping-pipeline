#!/usr/bin/env nextflow

nextflow.enable.dsl = 2
params.help = false

include { munge_sumstats } from './modules/munge_sumstats.nf'
include { compute_h2_L2 } from './modules/compute_h2_L2.nf'

def help_msg() {
    log.info """
    Pipeline to perform finemapping on SAIGE mixed linear model association results using PolyFun, SuSiE & FINEMAP
    Author: Vivek Appadurai | vivek.appadurai@regionh.dk
    USAGE: nextflow run main.nf

    Options:
    --assoc <saige_results.assoc> [A file of association test results from the SAIGE mixed linear model analysis]
    --n <10000> [Sample size of the association analysis]
    --help [Prints this message]
    --info <0.6> [Filter for imputation INFO scores]
    --maf <0.001> [Filter for minor allele frequency of SNPs]
    --out <my_out_prefix> [Prefix for output files]
    """
}

if(params.help) {
    help_msg()
    exit 0
}

log.info """
============================================================================
I B P _ F I N E M A P P I N G _ P I P E L I N E _ V 1 . 0 - N E X T F L O W
============================================================================
Association Test                 : $params.assoc
Sample Size                      : $params.n
Output Prefix                    : $params.out
Annotations                      : $params.annotation
Weights                          : $params.weights
LD Cache                         : $params.ld
=============================================================================
"""

workflow {

    // Step 1: Munge sumstats and store in parquet format for PolyFun

    process munge_sumstats {
        input:
            path sumstats from Channel.fromPath(params.assoc)
            val N from Channel.of(params.n)
            val out from Channel.of(params.out)
        output:
            path("${out}_munged.parquet") into sumstats_munged_ch
        script:
            """
            python munge_polyfun_sumstats.py \
                --sumstats $sumstats \
                --n $N \
                --out ${out}_munged.parquet
            """
    }

    // Step 2: Calculate per SNP h2 using L2-regularized S-LDSC

    process compute_h2_L2 {
        input:
            val out from Channel.of(params.out)
            path munged_sumstats from sumstats_munged_ch
            val annotation_files_prefix from Channel.of(params.annotation)
            val weight_files_prefix from Channel.of(params.weights)
            path annotation_files from Channel.fromPath("${params.annotation}.*")
            path weight_files from Channel.fromPath("${params.weights}.*")
        output:
            path("${out}.annot_coeff_ridge.*.txt") into bin_ch
            path("${out}.*.snpvar_ridge_constrained.gz") into bin_ch
            path("${out}.*.snpvar_ridge.gz") into bin_ch
            path("${out}.*.l2.M") into bin_ch
            path("${out}.*.bins.parquet") into bin_ch

        script:
            """
            python polyfun.py \
                --compute-h2-L2 \
                --output-prefix $out \
                --sumtats $munged_sumstats \
                --ref-ld-chr $annotation_file_prefix \
                --w-ld-chr $weight_file_prefix \
                --ld-wind-kb 1000 \
                --allow-missing
            """
    }

    // Step 3: Compute LD Scores for each SNP bin

    process compute_ld_scores {
        input:
            val out from Channel.of(params.out)
            val ld_dir from Channel.of(params.ld)
            path ld_files from Channel.fromPath(params.ld)
        output:
            path("${out}.*.l2.ldscore.paruqet") into ld_scores_ch
        script:
            """
            python polyfun.py \
                --compute-ldscores \
                --output-prefix $out \
                --ld-ukb \
                --ld-dir $ld_dir \
                --chr $chr
            """
    }

    // Step 4: Re-calculate per SNP h2 using S-LDSC to use as priors for functional finemapping

    process reestimate_snp_h2 {
        input:
            val out from Channel.of(params.out)
            path munged_sumstats from sumstats_munged_ch
            val weight_file_prefix from Channel.of(params.weights)
        output:
            path("${out}.*.snpvar_constrained.gz") into functional_priors_ch
        script:
            """
            python polyfun.py \
                --compute-h2-bins \
                --output-prefix $out \
                --sumstats $munged_sumstats \
                --w-ld-chr $weight_file_prefix
            """
    }
}