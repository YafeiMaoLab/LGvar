#!/bin/bash

folder=$1       
input_vcf=$2   
variants=$3    

if [ -z "$folder" ] || [ -z "$input_vcf" ] || [ ! -f "$input_vcf" ]; then
  echo "Usage: $0 <output_folder> <input_vcf> <variants_list>"
  exit 1
fi

mkdir -p "${folder}"

process_sv=false
process_indel=false
process_snv=false

if [ "$variants" = "all" ] || [ -z "$variants" ]; then
  process_sv=true
  process_indel=true
  process_snv=true
else
  IFS=',' read -r -a variant_arr <<< "$variants"
  for v in "${variant_arr[@]}"; do
    case "$v" in
      inv) process_sv=true ;;
      ins|del) process_indel=true ;;
      snv) process_snv=true ;;
    esac
  done
fi

# 1. SV
if $process_sv; then
    bcftools view -i 'abs(SVLEN)>=50 || SVTYPE="INV"' "$input_vcf" | bcftools sort -Oz -o ${folder}/sortSV.vcf.gz 2>/dev/null 
    tabix ${folder}/sortSV.vcf.gz
fi

# 2. INDEL
if $process_indel; then
    bcftools view -i 'abs(SVLEN)<50 && abs(SVLEN)>0' "$input_vcf" | bcftools sort -Oz -o ${folder}/sortindel.vcf.gz 2>/dev/null
    tabix ${folder}/sortindel.vcf.gz
fi

# 3. SNV
if $process_snv; then
    bcftools view -i 'SVTYPE="SNV"' "$input_vcf" | bcftools sort -Oz -o ${folder}/sortsnv.vcf.gz 2>/dev/null
    tabix ${folder}/sortsnv.vcf.gz
fi

wait
