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

readonly array_idx=$( get_array_task_id )
readonly chr=$( get_chr ${array_idx} )

readonly vep_dir="data/vep/hail-vep105"
readonly vep="${vep_dir}/ukb_wes_450k.july.qced.chr${chr}.vep_with_gnomad.ht"

readonly spliceai_dir="data/spliceai"
readonly spliceai_path="${spliceai_dir}/ukb_wes_450k.spliceai.chr${chr}.ht"

readonly out_dir="data/vep/annotated/v5"
readonly out_prefix="${out_dir}/ukb_wes_450k.july.qced.brava.v5.chr${chr}"

mkdir -p ${out_dir}

set_up_pythonpath_legacy  

python3 scripts/brava_create_annot.py -c ${chr} -v path/to/vep_table.gz -s path/to/SpliceAI.gz -w ./

