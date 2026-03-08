#!/usr/bin/sh

if [ $# -ne 3 ];then
	echo "Usage: bash $0 <gff> <fasta> <output_prefix>"
	exit 1 
fi


GFF=$1
FASTA=$2
OUTPUT_PREFIX=$3

if [ ! -f $FASTA.fai ];then
	echo "no fai for fasta input, ln -s this fasta ..."
	ln -s $FASTA
	FASTA=$(basename $FASTA)
fi


PFAM_DB="/share/home/zhanglab/user/suomingyu/src/PublicDataDownload/Pfam/Pfam-A.hmm"
GENE_NAME="MGAM"

TAG="MANE_Select"
#TAG="Ensembl_canonical"

# 1. Extract MANE Transcript coordinate ..
echo "Extract MANE Transcript coordinate for $GENE_NAME ..."
MAP_INFO=$(grep "$GENE_NAME;" "$GFF" | grep "$TAG" | awk '$3=="transcript" {print $1"\t"$4"\t"$5}' | head -n 1)

if [ -z "$MAP_INFO" ]; then
    echo "No MANE Transcript coordinate found for $GENE_NAME"
    exit 1
fi

echo $MAP_INFO
CHROM=$(echo "$MAP_INFO" | awk '{print $1}' )
START=$(echo "$MAP_INFO" | awk '{print $2}')
END=$(echo "$MAP_INFO" | awk '{print $3}')


# 2. Statistics of the N base (Gap)
echo "Statistics of the Gap ..."
samtools faidx $FASTA $CHROM":"$START"-"$END > $OUTPUT_PREFIX.$GENE_NAME.fa
## tr -c: complement, -d : delete
N_COUNT=$(grep -v ">" $OUTPUT_PREFIX.$GENE_NAME.fa | tr -cd 'Nn' | wc -c)
HAS_GAP=$([ "$N_COUNT" -gt 0 ] && echo "Yes($N_COUNT)" || echo "No")
echo "Gap?: $HAS_GAP"

# 3. Extract protein sequence (use gffread)
echo "Extract protein sequence ..."
##  -y : write a protein fasta file with the translation of CDS for each record
cat $GFF | grep "$GENE_NAME;" | grep $TAG > ${OUTPUT_PREFIX}.tmp.gff
gffread  ${OUTPUT_PREFIX}.tmp.gff  -g "$FASTA" -y "${OUTPUT_PREFIX}.${GENE_NAME}_prot.fa"  --adj-stop

if [ ! -s "${OUTPUT_PREFIX}.${GENE_NAME}_prot.fa" ]; then
    echo "Error: Protein extraction failed (File is empty)."
    exit 1
fi


# 4. Use hmmscan to predict domains
## --domblout : save parseable table of per-domain hits to file
## --noali : do not output alignments, to reduce output size
## -E : report models <= this E-value threshold in output  [10.0]
echo "Running hmmscan for Domain prediction ..."
hmmscan --domtblout "$OUTPUT_PREFIX.$GENE_NAME.domtab" --noali -E 1e-5 "$PFAM_DB" "${OUTPUT_PREFIX}.${GENE_NAME}_prot.fa" > /dev/null

# 5. Statistics of the domains
## Extract domain name 
DOMAIN_LIST=$(grep -v "^#" "$OUTPUT_PREFIX.$GENE_NAME.domtab" | awk '{print $1}')
TOTAL_COUNT=$(echo "$DOMAIN_LIST" | wc -w)
DISTINCT_COUNT=$(echo "$DOMAIN_LIST" | tr ' ' '\n' | sort -u | wc -l)
DETAIL_INFO=$(echo "$DOMAIN_LIST" | tr ' ' '\n' | awk '{count[$1]++} END {for (i in count) printf i":"count[i]", "}' | sed 's/, $//')

if [ "$TOTAL_COUNT" -eq 0 ]; then
    DISTINCT_COUNT=0
    DETAIL_INFO="None"
fi

echo -e "Chromosome\tStart\tEnd\tHas_Gap\tDistinct_Domain_Num\tTotal_Domain_Num\tDomain_Detail_List" > $OUTPUT_PREFIX.$GENE_NAME.stat.result.txt
echo -e "$CHROM\t$START\t$END\t$HAS_GAP\t$DISTINCT_COUNT\t$TOTAL_COUNT\t$DETAIL_INFO" >> $OUTPUT_PREFIX.$GENE_NAME.stat.result.txt

echo "Done!"
