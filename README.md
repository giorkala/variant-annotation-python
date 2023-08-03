# BRAVA SAIGE Annotation

This script merges VEP with SpliceAI info and generates SAIGE annotations according to BRaVa's guidelines. The required inputs for the annotations by VEP and SpliceAI are table format files (tab-delimited). 

## Requirements

This script requires `numpy` and `pandas` to be installed in your Python environment.

## Usage

```bash
python brava_create_annot.py -c <chromosome> -v <path_to_vep_table> -s <path_to_SpliceAI_table> -w <working_directory> --vep_snp_id_col <VEP SNP ID column> --vep_gene_col <VEP GENE ID column> --vep_lof_col <VEP LoF column> --vep_max_pop_col <VEP gnomAD_maxAF column> --vep_mane_select_col <VEP MANE SELECT column> --vep_revel_col <VEP REVEL column> --vep_cadd_phred_col <VEP CADD_PHRED column> --vep_consequence_col <VEP Consequence column>
```

## Arguments

* `-c/--chr`: Specify the chromosome to process. This is a required argument.
* `-v/--vep`: Specify the file path for the VEP annotation (tab-delimited). This is a required argument.
* `-s/--spliceai`: Specify the file path for the SpliceAI annotations. This is a required argument.
* `-w/--work_dir`: Specify the working directory. This is an optional argument, with the default being the current directory ('./').
* `--vep_snp_id_col`: Specify the SNP ID (chr:pos:ref:alt) column in the VEP table. This is a required argument.
* `--vep_gene_col`: Specify the GENE ID column in the VEP table. This is a required argument.
* `--vep_lof_col`: Specify the LoF column in the VEP table. This is a required argument.
* `--vep_max_pop_col`: Specify the gnomAD_maxAF column in the VEP table. This is a required argument.
* `--vep_revel_col`: Specify the REVEL column in the VEP table. This is a required argument.
* `--vep_cadd_phred_col`: Specify the CADD_PHRED column in the VEP table. This is a required argument.
* `--vep_consequence_col`: Specify the Consequence column in the VEP table. This is a required argument.

## Outputs

The script will generate the SAIGE annotation file for the given chromosome in the working directory. The output file name will follow the pattern `brava_SAIGE_group_chr<chromosome>.txt`.

