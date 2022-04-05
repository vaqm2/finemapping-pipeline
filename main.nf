#!/usr/bin/env nextflow

nextflow.enable.dsl = 2
import groovy.json.JsonSlurper

include { munge_sumstats } from './modules/munge_sumstats.nf' 
include { compute_h2_L2_calc_ld } from './modules/compute_h2_L2_calc_ld.nf'
include { reestimate_snp_h2 } from './modules/reestimate_snp_h2.nf'

def help_msg() {
    log.info """
    Pipeline to perform finemapping on SAIGE mixed linear model association results using PolyFun, SuSiE & FINEMAP
    Author: Vivek Appadurai | vivek.appadurai@regionh.dk
    USAGE: nextflow run main.nf

    Options:
    --assoc <saige_results.assoc> [A file of association test results from the SAIGE mixed linear model analysis]
    --n <10000> [Sample size of the association analysis]
    --help [Prints this message]
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
Summary Stats : $params.assoc
Sample Size   : $params.n
Output Prefix : $params.out
Annotations   : $params.annotation
Weights       : $params.weights
LD Cache      : $params.ld
=============================================================================
"""

String annotations   = new File(params.annotation).text
String weights       = new File(params.weights).text
def annotations_dict = new JsonSlurper().parseText(annotations)
def weights_dict     = new JsonSlurper().parseText(weights)
def weights_list     = []

annotations_ch = Channel.of(1..22) | map {
    a -> [file(annotations_dict[a.toString()]."ann").getBaseName(),
    annotations_dict[a.toString()]."ann",
    annotations_dict[a.toString()]."ld",
    annotations_dict[a.toString()]."m"
    ]
}

weights_ch = Channel.of(1..22) | map {
    a -> [file(weights_dict[a.toString()]).getBaseName(),
    weights_dict[a.toString()]
    ]
}

weights_dict.each { key, value -> weights_list.add(value) }

workflow {
    // Step 1: Munge sumstats and store in parquet format for PolyFun

    Channel.of(params.assoc) \
    | combine(Channel.of(params.n)) \
    | combine(Channel.of(params.out)) \
    | combine(Channel.of(params.munge_sumstats_script)) \
    | munge_sumstats \
    | set { sumstats_munged_ch }

    /* Steps 2 & 3: Calculate per SNP h2 using L2-regularized S-LDSC, 
    partition SNPs into bins and compute LD scores per bin */

    Channel.of(1..22) \
    | combine(annotations_ch, by: 0) \
    | combine(weights_ch, by: 0) \
    | combine(Channel.of(params.ld)) \
    | combine(Channel.of(params.polyfun_script)) \
    | combine(Channel.of(params.out)) \
    | combine(sumstats_munged_ch) \
    | combine(Channel.fromPath("${params.ld}/*").collect()) \
    | compute_h2_L2_calc_ld \
    | set { per_snp_h2_bin_ld_ch }

    /* Step 4: Re-calculate per SNP h2 using S-LDSC to use as priors for functional finemapping 
    Write SNP prior weights for finemapping to launch directory */

    Channel.of(params.polyfun_script) \
    | combine(Channel.of(params.out)) \
    | combine(sumstats_munged_ch) \
    | combine(Channel.of(params.weights)) \
    | combine(Channel.fromList(weights_list).collect()) \
    | combine(per_snp_h2_bin_ld_ch.collect()) \
    | reestimate_snp_h2
}