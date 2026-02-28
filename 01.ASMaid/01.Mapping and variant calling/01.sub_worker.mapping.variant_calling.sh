#!/bin/bash
#SBATCH --job-name=Proc_Sample
#SBATCH --cpus-per-task=10
#SBATCH --mem=100g
#SBATCH -t 300:00:00

# Parse arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --sample) SAMPLE="$2"; shift ;;
        --ref) REF="$2"; shift ;;
        --data_dir) DATA_DIR="$2"; shift ;;
        --mode) MODE="$2"; shift ;; # 'mem' or 'aln'
        --assembly) ASSEMBLY="$2"; shift ;;
        --root) ROOT="$2"; shift ;;
    esac
    shift
done

THREADS=$SLURM_CPUS_PER_TASK
RG="@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tLB:LIB\tPL:ILLUMINA"
CHRO_LIST="${ROOT}/config/${ASSEMBLY}.chro_list.txt"

# --- Step 1: Mapping & Merging ---
# Handle multiple fastq files in the directory
if [ "$MODE" == "mem" ]; then
    # Modern HGDP: Paired-end
    for fq1 in $(ls ${DATA_DIR}/*_1.clean.fastq.gz); do
        name=$(basename $fq1 | sed 's/_1.clean.fastq.gz//')
        fq2=$(echo $fq1 | sed 's/_1.clean.fastq.gz/_2.clean.fastq.gz/')
        bwa mem -t $THREADS -R "$RG" $REF $fq1 $fq2 | samtools view -@ $THREADS -bS > ${name}.bam
    done
else
    # Archaic: Single-end
    for fq in $(ls ${DATA_DIR}/*.fastq.gz); do
        name=$(basename $fq | sed 's/.fastq.gz//')
        bwa aln -t $THREADS -n 0.01 -l 16500 -o 2 $REF $fq > ${name}.sai
        bwa samse -r "$RG" $REF ${name}.sai $fq | samtools view -@ $THREADS -bS > ${name}.bam
        rm ${name}.sai
    done
fi

# Post-alignment processing for each bam
for b in *.bam; do
    prefix=${b%.bam}
    samtools view -@ $THREADS -bh -q 30 -F 2308 $b -o ${prefix}.filt.bam
    sambamba markdup -r -t $THREADS ${prefix}.filt.bam ${prefix}.rmdup.bam
    samtools sort -@ $THREADS ${prefix}.rmdup.bam -o ${prefix}.sorted.bam
    rm $b ${prefix}.filt.bam ${prefix}.rmdup.bam
done

# Merge all bams for this sample
samtools merge -@ $THREADS ${SAMPLE}.merge.bam *.sorted.bam
samtools sort -@ $THREADS ${SAMPLE}.merge.bam -o ${SAMPLE}.final.bam
samtools index -@ $THREADS ${SAMPLE}.final.bam
rm *.sorted.bam ${SAMPLE}.merge.bam

# --- Step 2: Variant Calling Branch ---
if [ "$MODE" == "mem" ]; then
    # Modern Samples -> DeepVariant
    for chro in $(cat $CHRO_LIST); do
        singularity run -B /share:/share \
            ${ROOT}/bin/deepvariant.simg /opt/deepvariant/bin/run_deepvariant \
            --model_type WGS --ref $REF --reads ${SAMPLE}.final.bam \
            --regions $chro --output_gvcf ${SAMPLE}.${chro}.g.vcf.gz --num_shards $THREADS
    done
    touch FinishDV
else
    # Archaic Samples -> GATK HaplotypeCaller
    # Ensure .dict exists
    DICT=$(echo $REF | sed 's/.fasta/.dict/')
    [ ! -f $DICT ] && samtools dict -o $DICT $REF
    
    for chro in $(cat $CHRO_LIST); do
        gatk HaplotypeCaller -R $REF -I ${SAMPLE}.final.bam \
            -L $chro -ERC GVCF -O ${SAMPLE}.${chro}.g.vcf.gz
    done
    touch FinishGATK
fi