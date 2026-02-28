#!/bin/bash
# Usage: bash 02_joint_genotype.sh <Assembly_Name> <Project_Root_Dir>

ASSEMBLY=$1
ROOT=$2
CHRO_LIST="${ROOT}/config/${ASSEMBLY}.chro_list.txt"
SAMPLE_LIST="${ROOT}/config/${ASSEMBLY}.sample_list.txt"
RESULT_DIR="${ROOT}/results/${ASSEMBLY}"

mkdir -p ${RESULT_DIR}/GLnexus_Merge && cd ${RESULT_DIR}/GLnexus_Merge

for chro in $(cat $CHRO_LIST); do
    echo "Joint calling for chromosome: $chro"
    
    # Collect all gVCF paths for this chromosome
    input_vcfs=""
    while IFS=$'\t' read -r pop name path mode; do
        sample_id="${pop}-${name}"
        vcf_path="${RESULT_DIR}/${sample_id}/${sample_id}.${chro}.g.vcf.gz"
        if [ -f "$vcf_path" ]; then
            input_vcfs="$input_vcfs $vcf_path"
        fi
    done < "$SAMPLE_LIST"

    # Run GLnexus
    mkdir -p $chro && cd $chro
    singularity run -B /share:/share \
        ${ROOT}/bin/glnexus_v1.4.1.sif glnexus_cli \
        --config DeepVariantWGS $input_vcfs \
        | bcftools view | bgzip -@ 10 -c > ${ASSEMBLY}.${chro}.merged.vcf.gz
    cd ..
done