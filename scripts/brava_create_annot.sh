#!/usr/bin/env bash
#
# @description annotate variants using Hail
# @depends quality controlled MatrixTables with variants.
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=annotate
#SBATCH --chdir=/well/lindgren/barney/brava_annotation
#SBATCH --output=logs/annotate.log
#SBATCH --error=logs/annotate.errors.log
#SBATCH --partition=short
#SBATCH --cpus-per-task 1
#SBATCH --array=1-23

#set -o errexit
#set -o nounset

#readonly array_idx=$( get_array_task_id )
#readonly chr=$( get_chr ${array_idx} )

chr=21

readonly vep_dir="data/vep"
readonly vep="${vep_dir}/ukb_wes_450k.chr${chr}.tsv.gz"

readonly spliceai_dir="/well/lindgren/barney/spliceai/out"
readonly spliceai_path="${spliceai_dir}/ukb_wes_450k.qced.v6.sites_only.${chr}.vcf"

readonly out_dir="data/out"
readonly out_prefix="${out_dir}/ukb_wes_450k.july.qced.brava.v5.2.chr${chr}.tsv"

readonly vep_snp_id_col="rsid"
readonly vep_gene_col="worst_csq_by_gene_canonical.gene_id"
readonly vep_lof_col="worst_csq_by_gene_canonical.lof"
readonly vep_max_pop_col="gnomAD_exomes_POPMAX_AF"
readonly vep_mane_select_col="worst_csq_by_gene_canonical.mane_select"
readonly vep_revel_col="worst_csq_by_gene_canonical.revel_score"
readonly vep_cadd_phred_col="worst_csq_by_gene_canonical.cadd_phred"
readonly vep_consequence_col="worst_csq_by_gene_canonical.most_severe_consequence"

mkdir -p ${out_dir}

set_up_pythonpath_legacy  

python3 scripts/brava_create_annot.py \
    --chr=$chr \
    --vep=$vep \
    --spliceai=$spliceai_path \
    --work_dir=$work_dir \
    --vep_snp_id_col=$vep_snp_id_col \
    --vep_gene_col=$vep_gene_col \
    --vep_lof_col=$vep_lof_col \
    --vep_max_pop_col=$vep_max_pop_col \
    --vep_mane_select_col=$vep_mane_select_col \
    --vep_revel_col=$vep_revel_col \
    --vep_cadd_phred_col=$vep_cadd_phred_col \
    --vep_consequence_col=$vep_consequence_col
