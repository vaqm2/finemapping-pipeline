#!/usr/bin/env nextflow

nextflow.enable.dsl = 2
import groovy.json.JsonSlurper

include { munge_sumstats } from './modules/munge_sumstats.nf' 

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
=======================================================================================================================
I B P _ F I N E M A P P I N G _ P I P E L I N E _ V 1 . 0 - N E X T F L O W
=======================================================================================================================
Association Test                 : $params.assoc
Sample Size                      : $params.n
Info Score Threshold             : $params.info
Minor Allele Frequency Threshold : $params.maf
Output Prefix                    : $params.out
=======================================================================================================================
"""

String baseline_ld_ukbb_anno   = new File(params.annotations).text
def baseline_ld_ukbb_anno_dict = new JsonSlurper().parseText(baseline_ld_ukbb_anno)
baseline_ld_ukbb_anno_ch = Channel.of(1..22)
    | map {a -> [
        a,
        file(baseline_ld_ukbb_anno_dict[a.to_string()]."ld").getBaseName(),
        baseline_ld_ukbb_anno_dict[a.to_string()]."ld",
        baseline_ld_ukbb_anno_dict[a.to_string()]."wt",
        baseline_ld_ukbb_anno_dict[a.to_string()]."ann",
        baseline_ld_ukbb_anno_dict[a.to_string()]."m"
    ]

    }

workflow {

    // Step 1: Munge sumstats and store in parquet format for PolyFun

    Channel.of(params.assoc) \
    | combine(Channel.of(params.n)) \
    | combine(Channel.of(params.info)) \
    | combine(Channel.of(params.maf)) \
    | combine(Channel.of(params.out)) \
    | combine(Channel.of(params.munge_sumstats_path)) \
    | munge_sumstats \
    | set(sumstats_munged_ch)

    // Step 2: Run PolyFun with L2-regularized S-LDSC to estimate per SNP heritabilities

    Channel.of(1..22) \
    | combine(baseline_ld_ukbb_anno_ch, by: 0) \
    | combine(sumstats_munged_ch) \
    | combine(Channel.of(params.out)) \
    | compute_l2_h2
    | set(snp_var_files_ch)
}