#!/bin/bash
log() {
    local level=$1
    local message=$2
    local timestamp=$(date "+%Y-%m-%d %H:%M:%S,%3N")

    echo -e "${timestamp} - ${level} - ${message}"
}

file=$1  ## sv file
hap1_sv_file=$2  ## hap1 sv bed
hap2_sv_file=$3  ## hap2 sv bed
output_bed=$4    ## output bed

bcftools query -f '%CHROM\t%POS\t%ID\t%SVTYPE\t[%INFO/SVLEN]\t[%GT]\t%INFO/TIG_REGION\n' "$file" | awk 'BEGIN{OFS="\t"} {
    print $1, $2-1, $3, $4, $5,$6, $7, $8
}' > results/sv1.bed 
sed -i '1i#CHROM\tPOS\tID\tSVTYPE\tSVLEN\tGT\tQUERY' results/sv1.bed

grep -E 'SDR|DUP|HighDup|TRANS|COMPLEX' "$hap1_sv_file" | awk 'BEGIN{OFS="\t"} { print $1, $2, $4, $5,$6, $8, $9}' > results/hap1.sv2.bed
grep -E 'SDR|DUP|HighDup|TRANS|COMPLEX' "$hap2_sv_file" | awk 'BEGIN{OFS="\t"} { print $1, $2, $4, $5,$6, $8, $9}' > results/hap2.sv2.bed
cat results/sv1.bed results/hap1.sv2.bed results/hap2.sv2.bed > "$output_bed"
rm results/sv1.bed results/hap1.sv2.bed results/hap2.sv2.bed
rm -r half/ align_hap1.flt.paf align_hap2.flt.paf align_hap1.final.paf align_hap2.final.paf
