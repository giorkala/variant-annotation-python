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
import warnings 
warnings.filterwarnings(action='ignore', message='All-NaN axis encountered')

parser = argparse.ArgumentParser(description="Merge VEP with SpliceAI info and generate annotations according to BRaVa's guidelines.")
parser.add_argument("--chr", "-c", help="which chromosome to precess", required=True, type=str)
parser.add_argument("--vep", "-v", help="file with VEP annotation (tab-delimited)", required=True, type=str)
parser.add_argument("--spliceai", "-s", help="file with SpliceAI annotation (tab-delimited; no header)", required=True, type=str)
parser.add_argument("--work_dir", "-w", help="working directory (optional)", type=str, default='./')
args = parser.parse_args()


PLOF_CSQS = ["transcript_ablation", "splice_acceptor_variant", "splice_donor_variant", "stop_gained", "frameshift_variant"]

MISSENSE_CSQS = ["stop_lost", "start_lost", "transcript_amplification", "inframe_insertion", "inframe_deletion", "missense_variant"]

SYNONYMOUS_CSQS = ["stop_retained_variant", "synonymous_variant"]

OTHER_CSQS = ["mature_miRNA_variant", "5_prime_UTR_variant", "3_prime_UTR_variant", "non_coding_transcript_exon_variant", "intron_variant",
              "NMD_transcript_variant", "non_coding_transcript_variant", "upstream_gene_variant",
              "downstream_gene_variant", "TFBS_ablation", "TFBS_amplification", "TF_binding_site_variant",
              "regulatory_region_ablation", "regulatory_region_amplification", "feature_elongation",
              "regulatory_region_variant", "feature_truncation", "intergenic_variant"]

def get_annot( variant ):
    """
    A function to assign variants to a class according to BRAVA's guidelines
    """
    cadd_cutoff = 28.1
    revel_cutoff = 0.773
    spliceai_cutoff = 0.20
    if variant.LoF == 'HC':
        annot = 'pLoF'
    elif variant.SpliceAI_DS_max >= spliceai_cutoff or variant.LoF == 'LC':
        # any variant with SpliceAI ≥ 0.2 + LOFTEE LC
        annot = 'damaging_missense'
    elif variant.Consequence in MISSENSE_CSQS and ( variant.REVEL >= revel_cutoff or variant.CADD_PHRED >= cadd_cutoff ):  
        # missense / start-loss / stop-loss with REVEL ≥ 0.773 AND/OR CADD ≥ 28.1 +
        annot = 'damaging_missense'
    elif variant.Consequence in MISSENSE_CSQS or variant.Consequence in {'inframe_deletion', 'inframe_insertion'}:
        # Other missense/protein-altering: missense / start-loss / stop-loss not categorised in (2) + in frame indels
        annot = 'other_missense'
    elif variant.Consequence in OTHER_CSQS:
        annot = 'non_coding'
    elif variant.Consequence in SYNONYMOUS_CSQS and variant.SpliceAI_DS_max < spliceai_cutoff:
        annot = 'synonymous'
    else:
        annot = 'na'
        
    return annot

def get_spliceAI_DS( value, sep='|' ):
    """
    Extract the last part from the SpliceAI field, which looks like
    T|OR11H1---ENSG0123---ENST0123---yes---protein_coding---NM_00123.1|0.00|0.01|0.00|0.03|22|0|-2|8
    Note that we need to consider cases with missing info, for which '.' is returned.
    """
    tmp = str.split(value, sep='|')
    if len(tmp)<2:
        return np.nan
    else:
        return np.max(np.array( tmp[2:6], dtype=float ))
    
if __name__=='__main__':
    print("Will make the annotation file for chrom-" + args.chr)

    cols_to_read = ["Allele","CADD_PHRED","CADD_RAW","CANONICAL", "Consequence","Gene",
                "gnomAD_AF","gnomAD_AFR_AF","gnomAD_AMR_AF","gnomAD_ASJ_AF","gnomAD_EAS_AF","gnomAD_FIN_AF","gnomAD_NFE_AF","gnomAD_OTH_AF","gnomAD_SAS_AF",
                "IMPACT","INTRON","LoF","LoF_filter","LoF_flags","LoF_info","MANE_PLUS_CLINICAL","MANE_SELECT","MAX_AF","MAX_AF_POPS",
                "PolyPhen","REVEL","SIFT","SIFT_pred","SYMBOL","SYMBOL_SOURCE","VARIANT_CLASS","VEP_canonical"]

    df1 = pd.read_csv( args.vep, sep='\t', low_memory=True, usecols=cols_to_read, compression='gzip') # encoding='cp1252',
    # generate SNP IDs as chr:pos:a1:a2 
    df1['SNPID'] = [':'.join(str.split(x)[:-1]) for x in df1.Allele.values]
    print("Unique SNP-IDs found in the VEP table:", np.unique(df1.SNPID).shape[0])

    indx = df1.MANE_SELECT != '.'
    print( "Rows to keep as MANE-Select: {0} out of {1}".format(sum(indx), len(indx) ))
    df1 = df1[indx]

    temp = df1[ ['gnomAD_AF','gnomAD_AFR_AF','gnomAD_AMR_AF','gnomAD_ASJ_AF','gnomAD_EAS_AF','gnomAD_FIN_AF','gnomAD_NFE_AF','gnomAD_OTH_AF','gnomAD_SAS_AF']].copy()
    temp.replace('.', np.nan, inplace=True)
    temp = temp.astype(float)

    df1['gnomAD_maxAF'] = [np.nanmax(x[1]) for x in temp.iterrows()]

    indx = np.isnan( df1.gnomAD_maxAF )
    print( "Rows to discard due to missing gnomAD AF info (nans): {0} out of {1}".format(sum(indx), len(indx) ))
    df1 = df1[~indx]

    indx = df1['gnomAD_maxAF'] <= 0.01 
    print( "Rows to discard due to AF filter: {0} out of {1}.".format(len(indx)-sum(indx), len(indx) ))
    df1 = df1[indx]

    # read input from SpliceAI
    df2 = pd.read_csv( args.spliceai, sep='\t', header=None)
    df2[1] = df2[1].astype(str)
    df2['SNPID'] = df2[[0,1,3,4]].agg(':'.join, axis=1)
    df2.set_index('SNPID', inplace=True)
    print("Unique SNP-IDs found with SpliceAI:", np.unique(df2.index).shape[0])

    # update the first dataframe with the DS_score
    temp = df2.loc[ df1.SNPID ][5].values
    df1['SpliceAI_DS_max'] = [get_spliceAI_DS(x) for x in temp]
    indx = np.isnan( df1.SpliceAI_DS_max )
    print( "Rows to discard due to missing SpliceAI info: {0} out of {1}".format(sum(indx), len(indx) ))
    df1 = df1[~indx]

    ## Now generate the annotations
    df3 = df1[['SNPID','Gene','LoF','Consequence','REVEL','CADD_PHRED','SpliceAI_DS_max']].copy()
    df3.replace('.', np.nan, inplace=True)
    df3 = df3.astype({'REVEL':float,'CADD_PHRED':float})

    annot_vars = {}
    for x in df3[:].iterrows():
        tmp = get_annot(x[1])
        annot_vars[ x[1].SNPID ] = ( tmp, x[1].Gene )

    df_annot = pd.DataFrame.from_dict( annot_vars, orient='index' )
    df_annot.columns = ['Class','Gene']

    total = 0
    for c in df_annot.Class.unique():
        x = df_annot.query( 'Class==@c' ).index
        total += len(x)
        print( c, len(x) )
        
    print("Genes found {0} out of {1} originally.".format( len(df_annot.Gene.unique()), len(df1.Gene.unique())  ))

    # write the resulting annotiation to disk
    print("Saving the annotation to", args.work_dir+"annotation_"+args.chr+".txt")
    with open( args.work_dir+"annotation_"+args.chr+".txt", 'w' ) as fout:
        for gene in df_annot.Gene.unique():
            temp = df_annot.query( 'Gene==@gene' )
            fout.write( gene+' var') # , end=''
            for snpid in temp.index:
                fout.write( ' ' + snpid)  
                
            fout.write( '\n' + gene + ' anno')
            for snpid in temp.index:
                fout.write( ' ' + temp.loc[snpid].Class ) 
            
            fout.write('\n')

# end-of-script
