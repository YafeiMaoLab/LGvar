#!/bin/bash

folder=$1  ## output folder
input_vcf=$2  ## input vcf file
#variants="all"  ## variants type
variants=$3

mkdir -p "${folder}"

process_inv=false
process_indel=false
process_snv=false

if [ "$variants" = "all" ]; then
  process_inv=true
  process_indel=true
  process_snv=true
else
  IFS=',' read -r -a variant_arr <<< "$variants"
  for v in "${variant_arr[@]}"; do
    case "$v" in
      inv) process_inv=true ;;
      ins) process_indel=true;;
      del) process_indel=true;;
      snv) process_snv=true ;;
    esac
  done
fi

awk -F'\t' '{
    if($0 ~ /^#/){  
        print $0
    } else{
        OFS="\t";
        $11=(length($4)-length($5) >0 ? length($4)-length($5) : length($5)-length($4))
        print $1, $2, $3, $4, $5, $6, $7, $8, $9,$10, $11
    }
}' "$input_vcf" > "${folder}/count.vcf"

cd "$folder" || exit 1

# inversion
if $process_inv; then
  awk -F'\t' '{
      if($0 ~ /^#/){ 
          print $0
      }else{
          OFS="\t";
          if($3 ~/INV/ && length($4) >= 50){
              print $1, $2, $3, $4, $5, $6, $7, $8, $9,$10
          }else if($11>=50){
              print $1, $2, $3, $4, $5, $6, $7, $8, $9,$10
          }
      }
  }' count.vcf > SV.vcf

  if [ -f "SV.vcf.gz" ]; then rm "SV.vcf.gz"; fi
  bgzip SV.vcf > /dev/null 2>&1
  bcftools sort SV.vcf.gz -o sortSV.vcf.gz > /dev/null 2>&1
  bcftools index -t sortSV.vcf.gz > /dev/null 2>&1
fi

# indel
if $process_indel; then
  awk -F'\t' '{
      if($0 ~ /^#/){ 
          print $0
      }else{
          OFS="\t";
          if(($3 ~/INS/ || $3 ~ /DEL/)){
              print $1, $2, $3, $4, $5, $6, $7, $8, $9,$10
          }
      }
  }' count.vcf > indel.vcf

  if [ -f "indel.vcf.gz" ]; then rm "indel.vcf.gz"; fi
  bgzip indel.vcf > /dev/null 2>&1
  bcftools sort indel.vcf.gz -o sortindel.vcf.gz > /dev/null 2>&1
  bcftools index -t sortindel.vcf.gz > /dev/null 2>&1
fi

# snv
if $process_snv; then
  awk -F'\t' '{
      if($0 ~ /^#/){ 
          print $0
      }else{
          OFS="\t";
          if($3 ~/SNV/){
              print $1, $2, $3, $4, $5, $6, $7, $8, $9,$10
          }
      }
  }' count.vcf > snv.vcf

  if [ -f "snv.vcf.gz" ]; then rm "snv.vcf.gz"; fi
  bgzip snv.vcf > /dev/null 2>&1
  bcftools sort snv.vcf.gz -o sortsnv.vcf.gz > /dev/null 2>&1
  bcftools index -t sortsnv.vcf.gz > /dev/null 2>&1
fi
