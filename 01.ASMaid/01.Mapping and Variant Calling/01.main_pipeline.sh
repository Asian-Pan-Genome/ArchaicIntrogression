#!/bin/bash

# Usage: bash 00_main_pipeline.sh <Assembly_Name> <Ref_Fasta_Path> <Project_Root_Dir> <Data_Root_Dir>
# Description: Main entry point for the mapping pipeline (Optimized for Nature Publication).
# It handles reference indexing, sample list generation, and SLURM job submission.

if [ $# -ne 4 ]; then
    echo "Error: Missing arguments."
    echo "Usage: $0 <Assembly_Name> <Ref_Fasta_Path> <Project_Root_Dir> <Data_Root_Dir>"
    exit 1
fi

ASSEMBLY=$1    # Assembly Name (e.g., C115-CMG02-Mat)
REF_PATH=$2    # Full path to the reference fasta
ROOT_DIR=$3    # Project working directory
DATA_ROOT=$4   # Root directory where raw fastq and lists are stored

# Define internal directory structure
RESULT_DIR="${ROOT_DIR}/results/${ASSEMBLY}"
CODE_DIR="${ROOT_DIR}/scripts"
CONFIG_DIR="${ROOT_DIR}/config"

mkdir -p "$RESULT_DIR"
mkdir -p "$CONFIG_DIR"

# Function: Generate Sample List
# This function creates a metadata TSV file used for batch processing.
generate_sample_config() {
    local hgdp_list_path="${DATA_ROOT}/HGDP.list/"
    local output_list="${CONFIG_DIR}/${ASSEMBLY}.sample_list.txt"
    
    echo "Generating sample configuration list at: $output_list"
    
    # Check if data directory exists before proceeding
    if [ ! -d "$hgdp_list_path" ]; then
        echo "Error: HGDP list directory not found at $hgdp_list_path"
        exit 1
    fi

    # Initialize/Clear the list file
    > "$output_list" 

    # Format: Population [tab] Sample_ID [tab] Data_Path [tab] Mapping_Mode
    # Mode 'mem': paired-end modern DNA; Mode 'aln': archaic DNA
    
    # 1. Randomly select modern HGDP samples for comparative analysis
    # 10 East African samples
    cat "${hgdp_list_path}Mbuti.list" | sort -R | head -3 | awk -v root="$DATA_ROOT" '{print $1"\t"$2"\t"root"/"$3"\tmem"}' >> "$output_list"
    cat "${hgdp_list_path}BantuKenya.list" | sort -R | head -3 | awk -v root="$DATA_ROOT" '{print $1"\t"$2"\t"root"/"$3"\tmem"}' >> "$output_list"
    cat "${hgdp_list_path}San.list" | awk -v root="$DATA_ROOT" '{print $1"\t"$2"\t"root"/"$3"\tmem"}' >> "$output_list"
    cat "${hgdp_list_path}BantuSouthAfrica.list" | sort -R | head -2 | awk -v root="$DATA_ROOT" '{print $1"\t"$2"\t"root"/"$3"\tmem"}' >> "$output_list"
    
    # 3 European and 3 Oceanian samples (Optional)
    cat "${hgdp_list_path}Sardinian.list" | sort -R | head -3 | awk -v root="$DATA_ROOT" '{print $1"\t"$2"\t"root"/"$3"\tmem"}' >> "$output_list"
    cat "${hgdp_list_path}Papuan.list" | sort -R | head -3 | awk -v root="$DATA_ROOT" '{print $1"\t"$2"\t"root"/"$3"\tmem"}' >> "$output_list"

    # 2. Add Archaic samples with 'aln' mode for damaged DNA
    echo -e "Nean\tAltai\t${DATA_ROOT}/Data/Altai\taln" >> "$output_list"
    echo -e "Nean\tCha\t${DATA_ROOT}/Data/Cha\taln" >> "$output_list"
    echo -e "Den\tDen\t${DATA_ROOT}/Data/Den\taln" >> "$output_list"
    
    echo "Configuration generated successfully."
}

# --- Step 1: Prepare Reference Genome ---
echo "Indexing reference genome..."
cd "$RESULT_DIR" || exit
cp "$REF_PATH" ./
REF_FILE=$(basename "$REF_PATH")
samtools faidx "$REF_FILE"

# Generate chromosome list (autosomes 1-22 + X)
chro_list="${CONFIG_DIR}/${ASSEMBLY}.chro_list.txt"
cat "${REF_FILE}.fai" | cut -f 1 | head -23 > "$chro_list"

# --- Step 2: Generate Config File (sample list: HGDP + Archaic) ---
generate_sample_config

# --- Step 3: Batch Job Submission via SLURM ---
while IFS=$'\t' read -r pop name data_path mode; do
    sample_id="${pop}-${name}"
    echo "Submitting SLURM job for: $sample_id (Mode: $mode)"
    
    mkdir -p "$sample_id" && cd "$sample_id" || exit
    
    # Submit the unified mapping template
    sbatch "${CODE_DIR}/00.sub_worker.mapping.variant_calling.sh" \
        --sample "$sample_id" \
        --ref "../$REF_FILE" \
        --data "$data_path" \
        --mode "$mode" \
        --assembly "$ASSEMBLY" \
        --root "$ROOT_DIR"
        
    cd ..
done < "${CONFIG_DIR}/${ASSEMBLY}.sample_list.txt"

echo "All jobs submitted. Check squeue for status."