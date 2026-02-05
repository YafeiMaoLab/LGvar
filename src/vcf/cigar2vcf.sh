#!/bin/bash
log() {
    local level=$1
    local message=$2
    local timestamp=$(date "+%Y-%m-%d %H:%M:%S,%3N")

    echo -e "${timestamp} - ${level} - ${message}"
}

if [ ! -d "temp" ]; then
    mkdir temp
fi

# check params
if [ $# -lt 8 ]; then
  echo "Usage: $0 <pacigar.txt> <output_prefix> <input_vcf> <reference_fasta> <query_fasta> <SDR_output> <python_script> <bed_output> [variant_type]"
  echo "  variant_type: Optional - specify a single variant type to process (snv, ins, del, inv, trans, sdr, dup, highdup, complex)"
  echo "  If no variant_type is specified, all types will be processed."
  exit 1
fi

# get params
filein=$1       # cigar
fileout=$2      # file prefix
input_vcf=$3    # vcf file
reference_fasta=$4 # reference genome
query_fasta=$5    # query genome
SDR_output=$6     # sdr results
python_script=$7  # scripts
bed_output=$8     # bed file
#variant_type="all"  # merge all variants
variant_type=$9

# VCF header
rm -f temp/vcf_header.txt
touch temp/vcf_header.txt
echo -e "##fileformat=VCFv4.2" >>temp/vcf_header.txt
time=$(date +'%Y%m%d')
echo -e "##fileDate=${time}" >>temp/vcf_header.txt
echo -e "##source=LGvar" >>temp/vcf_header.txt
echo -e "##reference=file:$reference_fasta" >>temp/vcf_header.txt
while IFS= read -r line; do
    contigname=$(echo "$line" |cut -f1)
    contiglen=$(echo "$line" |cut -f2)
    echo -e "##contig=<ID=${contigname},length=${contiglen}>" >>temp/vcf_header.txt
done < "${reference_fasta}.fai"
echo -e "##INFO=<ID=ID,Number=1,Type=String,Description=\"Variant ID\">" >>temp/vcf_header.txt
echo -e "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Variant type\">" >>temp/vcf_header.txt
echo -e "##INFO=<ID=SVLEN,Number=.,Type=Integer,Description=\"Variant length\">" >>temp/vcf_header.txt
echo -e "##INFO=<ID=TIG_REGION,Number=.,Type=String,Description=\"Contig region where variant was found (one per alt with h1 before h2 for homozygous calls)\">" >>temp/vcf_header.txt
echo -e "##INFO=<ID=QUERY_STRAND,Number=.,Type=String,Description=\"Strand of variant in the contig relative to the reference (order follows TIG_REGION)\">" >>temp/vcf_header.txt
echo -e "##INFO=<ID=INNER_REF,Number=.,Type=String,Description=\"Inversion inner breakpoint in reference coordinates (order follows TIG_REGION)\">" >>temp/vcf_header.txt
echo -e "##INFO=<ID=INNER_TIG,Number=.,Type=String,Description=\"Inversion inner breakpoint in contig coordinates (order follows TIG_REGION)\">" >>temp/vcf_header.txt
echo -e "##INFO=<ID=HOM_REF,Number=.,Type=String,Description=\"Perfect breakpoint homology (SV sequence vs reference). Format 'X,Y' where X homology upstream, and Y is homology downstream. Homology vs reference is often better for DEL.\">" >>temp/vcf_header.txt
echo -e "##INFO=<ID=HOM_TIG,Number=.,Type=String,Description=\"Perfect breakpoint homology (SV sequence vs contig). Format 'X,Y' where X homology upstream, and Y is homology downstream.  Homology vs contig is often better for INS.\">" >>temp/vcf_header.txt
echo -e "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">" >>temp/vcf_header.txt
echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample" >>temp/vcf_header.txt

# Preprocess
log "INFO" "Preprocessing input file..."

awk '{
    split($0, parts, "_")
    ref_part = parts[1]
    rest = ""
    for (i = 2; i <= length(parts); i++) {
        if (i == 2) {
            rest = parts[i]
        } else {
            rest = rest "_" parts[i]
        }
    }
    split(ref_part, ref, /:|-/)
    ref_chrom = ref[1]
    ref_start = ref[2]
    ref_end = ref[3]
    last_colon = length(rest)
    while (substr(rest, last_colon, 1) != ":" && last_colon > 0) {
        last_colon--
    }
    if (last_colon > 0) {
        query_chrom = substr(rest, 1, last_colon - 1)
        coord_part = substr(rest, last_colon + 1)
        split(coord_part, query, "-")
        query_start = query[1]
        query_end = query[2]
    } else {
        query_chrom = "N/A"
        query_start = "N/A"
        query_end = "N/A"
    }
    print $0 "\t" ref_chrom "\t" ref_start "\t" query_chrom "\t" query_start
}' "$filein" > temp/addout.txt

awk -F'\t' 'NR>=2 {
    if ($5 == "SNP_DEL") {
        printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", $9, $10 + $2 - 1, $10 + $2 + $4 - 1, $11, $12 + $3, $12 + $3, $1, $4, $5, $6, $7, $8
    } else if ($5 == "SNP_INS") {
        printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", $9, $10 + $2 - 1, $10 + $2 - 1, $11, $12 + $3 - 1, $12 + $3 + $4 -1, $1, $4, $5, $6, $7, $8
    } else if ($5 == "INDEL_DEL") {
        printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", $9, $10 + $2 - 1, $10 + $2 + $4 - 1, $11, $12 + $3, $12 + $3, $1, $4, $5, $6, $7, $8
    } else if ($5 == "INDEL_INS") {
        printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", $9, $10 + $2 - 1, $10 + $2 - 1, $11, $12 + $3, $12 + $3 + $4, $1, $4, $5, $6, $7, $8
    } else if ($5 == "SNP") {
        printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", $9, $10 + $2, $10 + $2 + $4 - 1, $11, $12 + $3, $12 + $3 + $4 - 1, $1, $4, $5, $6, $7, $8
    }
}' temp/addout.txt > temp/CIGARend.txt
sed -i 's/SNP_DEL/DEL/g' temp/CIGARend.txt
sed -i 's/INDEL_DEL/DEL/g' temp/CIGARend.txt
sed -i 's/SNP_INS/INS/g' temp/CIGARend.txt
sed -i 's/INDEL_INS/INS/g' temp/CIGARend.txt

# Specific variants process
if [ "$variant_type" = "all" ]; then
  process_snv=true
  process_ins=true
  process_del=true
  process_inv=true
  process_trans=true
  process_sdr=true
  process_dup=true
  process_highdup=true
  process_complex=true
else
  IFS=',' read -ra types <<< "$variant_type"
  
  process_snv=false
  process_ins=false
  process_del=false
  process_inv=false
  process_trans=false
  process_sdr=false
  process_dup=false
  process_highdup=false
  process_complex=false

  for type in "${types[@]}"; do
    case "$type" in
      snv)      process_snv=true ;;
      ins)      process_ins=true ;;
      del)      process_del=true ;;
      inv)      process_inv=true ;;
      trans)    process_trans=true ;;
      sdr)      process_sdr=true ;;
      dup)      process_dup=true ;;
      highdup)  process_highdup=true ;;
      complex)  process_complex=true ;;
      *)
        exit 1
        ;;
    esac
  done
fi

# SNV
if $process_snv; then
  log "INFO" "Processing SNVs..."
  awk -F'\t' '$9=="SNP"{print$0}' temp/CIGARend.txt > temp/oursnv.txt
  file="temp/oursnv.txt"
  awk -F'\t' -v OFS='\t' '{
    id = "ID=" $1 "-" $2 "-" $9 "-" $10 "-" $11 ";"
    sv_info = "SVTYPE=" $9 ";"
    region = "TIG_REGION=" $4 ":" $5 "-" $6 "," $4 ":" $5 "-" $6 ";"
    strand = "QUERY_STRAND=" $12 "," $12
    
    print $0, $1 "-" $2 "-" $9 "-" $10 "-" $11, id sv_info region strand, "GT", "1/0"
  }' ${file} > temp/testbe.txt
  awk -F'\t' '{printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", $1, $2, $13, $10, $11,".",".",$14,$15,$16}' temp/testbe.txt > temp/testchr1snvend.txt
  cat temp/testchr1snvend.txt > temp/endcigar.vcf
  sed 's/SNP/SNV/g' temp/endcigar.vcf > temp/end2cigar.vcf
  rm temp/testchr1snvend.txt
fi

# DEL, INS, INV
if $process_del || $process_ins || $process_inv; then
  log "INFO" "Processing SVs..."
  awk '$7 ~ "DEL" || $7 ~"INS" ||$7 ~"INV"{print $0}' $input_vcf | awk -F'\t' '($9!=0 || $10!=0){print $0}' > temp/SDRend.txt
  sed -i 's\no\+\g' temp/SDRend.txt
  sed -i 's/SDR_INV/INV/g' temp/SDRend.txt
  sed -i 's/SV_INV/INV/g' temp/SDRend.txt
  sed -i 's/SDR_INS/INS/g' temp/SDRend.txt
  sed -i 's/SV_INS/INS/g' temp/SDRend.txt
  sed -i 's/SDR_DEL/DEL/g' temp/SDRend.txt
  sed -i 's/SV_DEL/DEL/g' temp/SDRend.txt

  python $python_script --r $reference_fasta --q $query_fasta --i temp/SDRend.txt --o $SDR_output
fi

# INS
if $process_ins; then
  log "INFO" "Processing insertions..."
  awk -F'\t' '$9=="INS"{print$0}' temp/CIGARend.txt > temp/ourins.txt
  file="temp/ourins.txt"
  awk -F'\t' -v OFS='\t' '{
    id = $1 "-" $2 "-" $9 "-" $8;
    info = "ID=" id ";SVTYPE=" $9 ";SVLEN=" $8 ";TIG_REGION=" $4 ":" $5 "-" $6 "," $4 ":" $5 "-" $6 ";QUERY_STRAND=" $12 "," $12;
    print $0, id, info, "GT", "1/0"
  }' ${file} > temp/testbe.txt
  awk -F'\t' '{printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", $1, $2, $13, $10, $11,".",".",$14,$15,$16}' temp/testbe.txt > temp/testchr1insend.txt

  rm temp/ourins.txt
  rm temp/testbe.txt
  cat <(awk 'index($3, "INS"){print$0}' "$SDR_output") temp/testchr1insend.txt > temp/ourinsend.txt
  rm temp/testchr1insend.txt
fi

# DEL
if $process_del; then
  log "INFO" "Processing deletions..."
  awk -F'\t' '$9=="DEL"{print$0}' temp/CIGARend.txt > temp/ourdel.txt
  file="temp/ourdel.txt"
  awk -F'\t' -v OFS='\t' '{
    id = $1 "-" $2 "-" $9 "-" $8;
    info = "ID=" id ";SVTYPE=" $9 ";SVLEN=" "-" $8 ";TIG_REGION=" $4 ":" $5 "-" $6 "," $4 ":" $5 "-" $6 ";QUERY_STRAND=" $12 "," $12;
    print $0, id, info, "GT", "1/0"
  }' ${file} > temp/testbe.txt
  awk -F'\t' '{printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", $1, $2, $13, $10, $11,".",".",$14,$15,$16}' temp/testbe.txt > temp/testchr1delend.txt
  rm temp/ourdel.txt
  rm temp/testbe.txt
  cat <(awk 'index($3, "DEL"){print$0}' "$SDR_output") temp/testchr1delend.txt > temp/ourdelend.txt
  rm temp/testchr1delend.txt
fi

# INV
if $process_inv; then
  log "INFO" "Processing inversions..."
  cat <(awk 'index($3, "INV"){print$0}' "$SDR_output") > temp/ourinvend.txt
fi

# Merge
log "INFO" "Merge results..."
> temp/variants.vcf
if $process_snv; then
  cat temp/end2cigar.vcf >> temp/variants.vcf
fi
if $process_ins; then
  cat temp/ourinsend.txt >> temp/variants.vcf
fi
if $process_del; then
  cat temp/ourdelend.txt >> temp/variants.vcf
fi
if $process_inv; then
  cat temp/ourinvend.txt >> temp/variants.vcf
fi

cat temp/vcf_header.txt temp/variants.vcf > "$fileout"
hap=$(basename "$filein" | cut -c 1-4)

# Bed results
log "INFO" "Creating BED files..."
> temp/LGvar.bed
echo -e "#CHROM\tPOS\tEND\tID\tSVTYPE\tSVLEN\tHAP\tGT\tQUERY" > temp/LGvar.bed

if $process_snv; then
  paste <(grep 'SNV' temp/variants.vcf | awk -v hap=${hap} 'OFS="\t"{print $1,$2,$3,$2+1,"SNV",1,hap,"1|."}') <(grep 'SNV' temp/variants.vcf |awk -F'TIG_REGION=' '{print $2}' | awk -F',' '{print $1}') > temp/SNV.bed
  cat temp/SNV.bed >> temp/LGvar.bed
fi

if $process_ins; then
  paste <(grep 'INS' temp/variants.vcf |awk -v hap=${hap} 'OFS="\t"{print $1,$2,$3}') <(grep 'INS' temp/variants.vcf |awk -F'[-\t]' -v hap=${hap} 'OFS="\t"{print $2+1,"INS",$6,hap,"1|."}') <(grep 'INS' temp/variants.vcf|awk -F'TIG_REGION=' '{print $2}' |awk -F',' '{print $1}') > temp/INS.bed
  cat temp/INS.bed >> temp/LGvar.bed
fi

if $process_del; then
  paste <(grep 'DEL' temp/variants.vcf |awk -v hap=${hap} 'OFS="\t"{print $1,$2,$3}') <(grep 'DEL' temp/variants.vcf |awk -F'[-\t]' -v hap=${hap} 'OFS="\t"{print $2+$6,"DEL",$6,hap,"1|."}') <(grep 'DEL' temp/variants.vcf |awk -F'TIG_REGION=' '{print $2}' |awk -F',' '{print $1}') > temp/DEL.bed
  cat temp/DEL.bed >> temp/LGvar.bed
fi

if $process_trans; then
  grep 'TRANS' "$input_vcf" |awk -v hap=${hap} 'OFS="\t"{print $1,$2,$1"-"$2"-TRANS-"$9,$3,"TRANS",$9,hap,"1|.",$4":"$5"-"$6}' > temp/TRANS.bed
  cat temp/TRANS.bed >> temp/LGvar.bed
fi

if $process_sdr; then
  grep 'SDR_NM' "$input_vcf" |awk -v hap=${hap} 'OFS="\t"{print $1,$2,$1"-"$2"-SDR-"$9,$3,"SDR",$9,hap,"1|.",$4":"$5"-"$6}' > temp/SDR.bed
  grep 'SV_NM' "$input_vcf" |awk -v hap=${hap} 'OFS="\t"{print $1,$2,$1"-"$2"-SDR-"$9,$3,"SDR",$9,hap,"1|.",$4":"$5"-"$6}' > temp/SV_NM.bed
  #grep 'SDR_NM(INV)' "$input_vcf" |awk -v hap=${hap} 'OFS="\t"{print $1,$2,$1"-"$2"-SDR-"$9,$3,"SDR(INV)",$9,hap,"1|.",$4":"$5"-"$6}' > SDR_NMINV.bed
  #grep 'SV_NM(INV)' "$input_vcf" |awk -v hap=${hap} 'OFS="\t"{print $1,$2,$1"-"$2"-SDR-"$9,$3,"SDR(INV)",$9,hap,"1|.",$4":"$5"-"$6}' > SV_NMINV.bed
  grep 'INV-INV' "$input_vcf" |awk -v hap=${hap} 'OFS="\t"{print $1,$2,$1"-"$2"-SDR-"$9,$3,"INV-INV",$9,hap,"1|.",$4":"$5"-"$6}' > temp/INV_INV.bed
  #cat SDR.bed SV_NM.bed SDR_NMINV.bed SV_NMINV.bed INV_INV.bed >> LGvar.bed
  cat temp/SDR.bed temp/SV_NM.bed temp/INV_INV.bed >> temp/LGvar.bed
fi

if $process_dup; then
  grep 'DUP' "$input_vcf" |awk -v hap=${hap} 'OFS="\t"{print $1,$2,$1"-"$2"-DUP-"$9,$3,"DUP",$9,hap,"1|.",$4":"$5"-"$6}' > temp/DUP.bed
  cat temp/DUP.bed >> temp/LGvar.bed
fi

if $process_highdup; then
  grep 'high-dup' "$input_vcf" |awk -v hap=${hap} 'OFS="\t"{print $1,$2,$1"-"$2"-HighDup",$3,"HighDup",".",hap,"1|.","."}' > temp/Highdup.bed
  cat temp/Highdup.bed >> temp/LGvar.bed
fi

if $process_inv; then
  grep 'SDR_INV' "$input_vcf" |awk -v hap=${hap} 'OFS="\t"{print $1,$2,$1"-"$2"-INV-"$9,$3,"INV",$9,hap,"1|.",$4":"$5"-"$6}' > temp/SDR_INV.bed
  grep 'SV_INV' "$input_vcf" |awk -v hap=${hap} 'OFS="\t"{print $1,$2,$1"-"$2"-INV-"$9,$3,"INV",$9,hap,"1|.",$4":"$5"-"$6}' > temp/SV_INV.bed
  cat temp/SDR_INV.bed temp/SV_INV.bed >> temp/LGvar.bed
fi

if $process_complex; then
  grep 'COMPLEX' "$input_vcf" |awk -v hap=${hap} 'OFS="\t"{print $1,$2,$1"-"$2"-SDR_COMPLEX-"$9,$3,"SDR_COMPLEX",$9,hap,"1|.",$4":"$5"-"$6}' > temp/complex.bed
  cat temp/complex.bed >> temp/LGvar.bed
fi

awk 'OFS="\t"{print $1,$2,$4,$3,$5,$6,$7,$8,$9}' temp/LGvar.bed > $bed_output

if [ -f "results/LGvarhap2.bed" ]; then
    awk 'BEGIN{OFS="\t"} 
        /^#/ {print; next}  
        NF >= 8 {$8 = ".|1"}  
        {print}              
    ' results/LGvarhap2.bed > temp/temp.bed

    mv temp/temp.bed results/LGvarhap2.bed
fi

rm -r temp

log "INFO" "Done! Created bed file for $variant_type"
