#!/usr/bin/env python
# coding: utf-8
"""
A script to merge VEP with SpliceAI info and generate annotations according to BRaVa's guidelines.

The required input for the annotations by VEP and SpliceAI is a table format (tab-delimited), e.g. by
`bcftools +split-vep -l $vep_chr | cut -f 2 | tr '\n' '\t' | sed -E 's/(.*)\t/\1\n/' > $table_vep`
`bcftools +split-vep -f '%CHROM %POS %REF %ALT %CSQ\n' -d -A tab $vep_chr >> $table_vep`
`bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/SpliceAI\n' $SpliceAI_chr > $table_spliceAI`

The output follows BRaVa's format which is (2 lines per gene):
ENSG00000187634 var chr1:943315:T:C chr1:962890:T:A ...
ENSG00000187634 anno damaging_missense non_coding ...

Run as
`python brava_create_annot.py -c 22 -v path/to/vep_table.gz -s path/to/SpliceAI.gz -w ./`

####################
GK - July 12th, 2023
"""

import numpy as np
import pandas as pd
import argparse

parser = argparse.ArgumentParser(description="Merge VEP with SpliceAI info and generate annotations according to BRaVa's guidelines.")
parser.add_argument("--chr", "-c", help="which chromosome to precess", required=True, type=str)
parser.add_argument("--vep", "-v", help="file with VEP annotation (tab-delimited)", required=True, type=str)
parser.add_argument("--spliceai", "-s", help="VCF file with SpliceAI annotations", required=True, type=str)
parser.add_argument("--work_dir", "-w", help="working directory (optional)", type=str, default='./')

# Columns to read (VEP)
parser.add_argument("--vep_snp_id_col", help="SNPID (chr:pos:ref:alt) column in VEP table", required=True, type=str)
parser.add_argument("--vep_gene_col", help="GENEID column in VEP table", required=True, type=str)
parser.add_argument("--vep_lof_col", help="LoF column in VEP table", required=True, type=str)
parser.add_argument("--vep_max_pop_col", help="gnomAD_maxAF column in VEP table", required=True, type=str)
parser.add_argument("--vep_mane_select_col", help="MANE SELECT column in VEP table", required=True, type=str)
parser.add_argument("--vep_revel_col", help="REVEL column in VEP table", required=True, type=str)
parser.add_argument("--vep_cadd_phred_col", help="CADD_PHRED column in VEP table", required=True, type=str)
parser.add_argument("--vep_consequence_col", help="Consequence column in VEP table", required=True, type=str)

args = parser.parse_args()

PLOF_CSQS = ["transcript_ablation", "splice_acceptor_variant", "splice_donor_variant", "stop_gained", "frameshift_variant"]

MISSENSE_CSQS = ["stop_lost", "start_lost", "transcript_amplification", "inframe_insertion", "inframe_deletion", "missense_variant"]

SYNONYMOUS_CSQS = ["stop_retained_variant", "synonymous_variant"]

OTHER_CSQS = ["mature_miRNA_variant", "5_prime_UTR_variant", "3_prime_UTR_variant", "non_coding_transcript_exon_variant", "intron_variant",
              "NMD_transcript_variant", "non_coding_transcript_variant", "upstream_gene_variant",
              "downstream_gene_variant", "TFBS_ablation", "TFBS_amplification", "TF_binding_site_variant",
              "regulatory_region_ablation", "regulatory_region_amplification", "feature_elongation",
              "regulatory_region_variant", "feature_truncation", "intergenic_variant"]

def get_annotation(variant):
    """
    A function to assign variants to a class according to BRAVA's guidelines
    """
    cadd_cutoff = 28.1
    revel_cutoff = 0.773
    spliceai_cutoff = 0.20


    if variant[args.vep_lof_col] == 'HC':
        annot = 'pLoF'
    elif variant[args.vep_consequence_col] in MISSENSE_CSQS and (variant[args.vep_revel_col] >= revel_cutoff or variant[args.vep_cadd_phred_col] >= cadd_cutoff):
        annot = 'damaging_missense'
    elif variant['max_DS'] >= spliceai_cutoff:
        annot = 'damaging_missense'
    elif variant[args.vep_lof_col] == 'LC':
        annot = 'damaging_missense'
    elif variant[args.vep_consequence_col] in MISSENSE_CSQS or variant[args.vep_consequence_col] in {'inframe_deletion', 'inframe_insertion'}:
        annot = 'other_missense'
    elif variant[args.vep_consequence_col] in OTHER_CSQS:
        annot = 'non_coding'
    elif (variant[args.vep_consequence_col] in SYNONYMOUS_CSQS) and (variant['max_DS'] < spliceai_cutoff or pd.isna(variant['max_DS'])):
        annot = 'synonymous'
    else:
        annot = ''

    return annot

def get_spliceAI_DS( value ):
    """
    Extract the DS_max from the SpliceAI field, which looks like
    T|OR11H1---ENSG0123---ENST0123---yes---protein_coding---NM_00123.1|0.00|0.01|0.00|0.03|22|0|-2|8
    Note that we need to consider cases with missing info, for which '.' is returned.
    IMPORTANT, we must also consider the case of multiple spliceAI scores per variant such as:
    SpliceAI=C|ENSG00000276612---ENSG00000276612.3---ENST00000623476.2---yes---protein_coding---|0.00|0.00|0.00|0.00|15|34|47|7,C|GATD3B---ENSG00000280071.4---ENST00000620528.5---yes---protein_coding---NM_001363758.2|0.00|0.00|0.00|0.00|5|-9|47|7
    """

    # Check if the value is missing
    if pd.isna(value):
        return np.nan

    scores = []

    # Split spliceAI scores
    spliceai_split = value.split('|')
    ncols_per_annotation=9
    n_genes = np.floor(len(spliceai_split) / ncols_per_annotation).astype(int)

    for gene_id in range(n_genes):
        
        # Return the maximum of the relevant components
        start_idx = gene_id * 9
        DS_max = np.max(np.array(spliceai_split[start_idx+2:start_idx+6], dtype=float))

        split_spliceai_info = spliceai_split[start_idx+1].split('---')

        gene_symbol, gene_id, transcript_id = split_spliceai_info[0], split_spliceai_info[1], split_spliceai_info[2]
        scores.append((gene_id, DS_max))

    return scores

def compare_with_hail(vep_df, hail_file):
    hail_vep_df = pd.read_csv(hail_file, sep="\t", na_values="NA")
    hail_vep_df["consequence_category"] = hail_vep_df["consequence_category"].replace(np.nan, '')
    print(f"Pandas VEP df size: {vep_df.shape[0]}")
    print(f"Hail VEP df size: {hail_vep_df.shape[0]}")
    vep_df = vep_df.rename(columns={'rsid':'varid', 'worst_csq_by_gene_canonical.gene_id': 'gene_id'})
    print(vep_df, hail_vep_df)
    merged_df = pd.merge(vep_df, hail_vep_df, on=['varid', 'gene_id'])
    print(f"Combined df size: {merged_df.shape[0]}")
    merged_df['consequence_comparison'] = merged_df.apply(
        lambda row: (row['annotation'] == row['consequence_category']),
        axis=1
    )

    print(merged_df[["annotation", "consequence_category"]][merged_df['consequence_comparison'] != True])

    print(merged_df['consequence_comparison'].value_counts())


    #print(f"Combined VEP df size: {}")


def write_saige_file(vep_df):

    saige_file = f"{args.work_dir}brava_SAIGE_group_chr{args.chr}.txt"

    # write the SAIGE group file to disk
    print(f"Saving the annotation to {saige_file}")
    with open(saige_file, 'w') as fout:
        for gene_id, variants in vep_df.groupby(args.vep_gene_col):
            annotated = ~variants.annotation.isna()

            snps  = ' '.join([snp for snp in variants[args.vep_snp_id_col][annotated].values])
            annos = ' '.join([anno for anno in variants.annotation[annotated].values])

            fout.write(f"{gene_id} var {snps}\n")
            fout.write(f"{gene_id} anno {annos}\n")
    
if __name__=='__main__':
    print("Will make the annotation file for chrom-" + args.chr)

    vep_cols_to_read = [args.vep_snp_id_col, args.vep_gene_col, args.vep_lof_col, args.vep_max_pop_col, args.vep_mane_select_col,
                        args.vep_revel_col, args.vep_cadd_phred_col, args.vep_lof_col, args.vep_consequence_col]

    vep_df = pd.read_csv(args.vep, sep='\t', usecols=vep_cols_to_read, na_values='.')

    # gnomAD maxAF FILTER
    num_total_rows = vep_df.shape[0]
    vep_df = vep_df[(vep_df[args.vep_max_pop_col] <= 0.01) | (vep_df[args.vep_max_pop_col].isna())]
    num_kept_rows = vep_df.shape[0]
    print(f"gnomAD FILTER: {num_kept_rows} out of {num_total_rows} remaining")

    spliceai_df = pd.read_csv(args.spliceai, sep='\t', comment='#', header=None, names=["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"], na_values='.', index_col='ID')

    print("Unique SNP-IDs found with SpliceAI:", np.unique(spliceai_df.index).shape[0])

    spliceai_df["spliceai_pairs"] = spliceai_df["INFO"].map(get_spliceAI_DS)
    spliceai_df = spliceai_df.explode("spliceai_pairs")
    spliceai_df[["GENE", "max_DS"]] = pd.DataFrame(spliceai_df['spliceai_pairs'].values.tolist(), index=spliceai_df.index)
    spliceai_df = spliceai_df.drop(["INFO", "spliceai_pairs"], axis='columns')

    # We don't want the ENSEMBL gene ID version - drop the .VERSION
    spliceai_df['GENE'] = spliceai_df['GENE'].astype(str).apply(lambda x: x.split('.')[0])

    # nan max DS_score FILTER
    num_total_rows = spliceai_df.shape[0]
    spliceai_df = spliceai_df[~spliceai_df["max_DS"].isna()]
    num_kept_rows = spliceai_df.shape[0]
    print(f"DS_Score FILTER: {num_kept_rows} out of {num_total_rows} remaining")

    ## Now generate the annotations
    vep_df = vep_df.merge(spliceai_df, left_on=[args.vep_snp_id_col, args.vep_gene_col], right_on=["ID", "GENE"], how='left')
    vep_df["annotation"] = vep_df.apply(get_annotation, axis=1)

    # no annotation FILTER
    num_total_rows = vep_df.shape[0]
    vep_df = vep_df[~(vep_df['annotation'] == '')]
    num_kept_rows = vep_df.shape[0]
    print(f"empty annotation FILTER: {num_kept_rows} out of {num_total_rows} remaining") 

    print(vep_df)
    vep_df["annotation"].value_counts()
    write_saige_file(vep_df)
    #compare_with_hail(vep_df, f"/well/lindgren/barney/brava_annotation/data/vep/annotated/v5/ukb_wes_450k.july.qced.brava.v5.3.chr{args.chr}.worst_csq_by_gene_canonical.txt.gz")
