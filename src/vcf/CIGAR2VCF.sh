#!/bin/bash
set -euo pipefail

# ---------- logging ----------
log() {
    local level=$1
    local message=$2
    echo -e "$(date '+%Y-%m-%d %H:%M:%S,%3N') - ${level} - ${message}"
}

# ---------- usage ----------
if [ $# -lt 8 ]; then
    echo "Usage: $0 <pacigar.txt> <output_vcf> <input_vcf> <reference_fasta> <query_fasta> <sdr_output> <python_script> <bed_output> [variant_type]"
    echo "  variant_type: all (default) | comma-separated list of: snv,ins,del,inv,trans,sdr,dup,highdup,complex"
    exit 1
fi

# ---------- parameters ----------
readonly cigar_in=$1          # input CIGAR table
readonly vcf_out=$2           # output VCF file
readonly input_vcf=$3         # upstream VCF (for SV/SDR/DUP etc.)
readonly ref_fasta=$4         # reference FASTA
readonly query_fasta=$5       # query FASTA
readonly sdr_output=$6        # SDR python script output
readonly python_script=$7     # python processing script
readonly bed_output=$8        # output BED file
readonly variant_type=${9:-all}

# ---------- temp directory ----------
TMPDIR=$(mktemp -d -t cigar2vcf_XXXXXX)
trap 'rm -rf "$TMPDIR"' EXIT

# ---------- parse variant types ----------
declare -A PROCESS
for t in snv ins del inv trans sdr dup highdup complex; do
    PROCESS[$t]=false
done

if [ "$variant_type" = "all" ]; then
    for t in "${!PROCESS[@]}"; do PROCESS[$t]=true; done
else
    IFS=',' read -ra requested <<< "$variant_type"
    for t in "${requested[@]}"; do
        if [[ -v PROCESS[$t] ]]; then
            PROCESS[$t]=true
        else
            log "ERROR" "Unknown variant type: '$t'"
            exit 1
        fi
    done
fi

hap=$(basename "$cigar_in" | cut -c1-4)

# ==========================================================
# 1. Build VCF header
# ==========================================================
log "INFO" "Building VCF header..."
{
    echo "##fileformat=VCFv4.2"
    echo "##fileDate=$(date +'%Y%m%d')"
    echo "##source=LGvar"
    echo "##reference=file:${ref_fasta}"
    awk '{printf "##contig=<ID=%s,length=%s>\n", $1, $2}' "${ref_fasta}.fai"
    echo '##INFO=<ID=ID,Number=1,Type=String,Description="Variant ID">'
    echo '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Variant type">'
    echo '##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Variant length">'
    echo '##INFO=<ID=TIG_REGION,Number=.,Type=String,Description="Contig region where variant was found (one per alt with h1 before h2 for homozygous calls)">'
    echo '##INFO=<ID=QUERY_STRAND,Number=.,Type=String,Description="Strand of variant in the contig relative to the reference (order follows TIG_REGION)">'
    echo '##INFO=<ID=INNER_REF,Number=.,Type=String,Description="Inversion inner breakpoint in reference coordinates (order follows TIG_REGION)">'
    echo '##INFO=<ID=INNER_TIG,Number=.,Type=String,Description="Inversion inner breakpoint in contig coordinates (order follows TIG_REGION)">'
    echo '##INFO=<ID=HOM_REF,Number=.,Type=String,Description="Perfect breakpoint homology (SV sequence vs reference). Format '"'"'X,Y'"'"' where X homology upstream, and Y is homology downstream. Homology vs reference is often better for DEL.">'
    echo '##INFO=<ID=HOM_TIG,Number=.,Type=String,Description="Perfect breakpoint homology (SV sequence vs contig). Format '"'"'X,Y'"'"' where X homology upstream, and Y is homology downstream.  Homology vs contig is often better for INS.">'
    echo '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">'
    echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample"
} > "$TMPDIR/vcf_header.txt"

# ==========================================================
# 2. Preprocess CIGAR table → resolve chrom / coord columns
#    then map type names and format into standardised columns
#
#    Output columns (cigar_parsed.txt):
#    1:ref_chrom 2:ref_pos_start 3:ref_pos_end
#    4:qry_chrom 5:qry_pos_start 6:qry_pos_end
#    7:orig_id   8:svlen  9:svtype  10:ref_seq  11:qry_seq  12:strand
# ==========================================================
log "INFO" "Preprocessing CIGAR input..."

awk -F'\t' '
BEGIN { OFS="\t"
    rename["SNP_DEL"]  = "DEL"
    rename["INDEL_DEL"]= "DEL"
    rename["SNP_INS"]  = "INS"
    rename["INDEL_INS"]= "INS"
    rename["SNP"]      = "SNP"
}
NR == 1 { next }
{
    id = $1
    if (!match(id, /:[0-9]+-[0-9]+_/)) {
        next
    }

    split_pos = RSTART + RLENGTH - 1
    ref_part   = substr(id, 1, split_pos - 1)
    query_part = substr(id, split_pos + 1)

    if (match(ref_part, /:[0-9]+-[0-9]+$/)) {
        ref_chrom = substr(ref_part, 1, RSTART - 1)
        split(substr(ref_part, RSTART + 1), ra, "-")
    }

    if (match(query_part, /:[0-9]+-[0-9]+$/)) {
        qry_chrom = substr(query_part, 1, RSTART - 1)
        split(substr(query_part, RSTART + 1), qa, "-")
    }

    ref_s = ra[1] + 0
    qry_s = qa[1] + 0
    
    typ = $5
    len = $4 + 0

    if (typ == "SNP_DEL" || typ == "INDEL_DEL") {
        r1 = ref_s + $2 - 1;  r2 = ref_s + $2 + len - 1
        q1 = qry_s + $3;       q2 = q1
    } else if (typ == "SNP_INS" || typ == "INDEL_INS") {
        r1 = ref_s + $2 - 1;  r2 = r1
        q1 = qry_s + $3 - 1;  q2 = qry_s + $3 + len - 1
    } else if (typ == "SNP") {
        r1 = ref_s + $2;       r2 = ref_s + $2 + len - 1
        q1 = qry_s + $3;       q2 = qry_s + $3 + len - 1
    } else {
        next
    }

    svtype = (typ in rename) ? rename[typ] : typ
    print ref_chrom, r1, r2, qry_chrom, q1, q2, $1, len, svtype, $6, $7, $8
}
' "$cigar_in" > "$TMPDIR/cigar_parsed.txt"

# ==========================================================
# 3. Process SNV  (SNP → SNV)
# ==========================================================
vcf_snv="$TMPDIR/snv.vcf"
if ${PROCESS[snv]}; then
    log "INFO" "Processing SNVs..."
    awk -F'\t' '
    $9 == "SNP" {
        id   = $1 "-" $2 "-SNV-" $10 "-" $11
        info = "ID=" id ";SVTYPE=SNV;TIG_REGION=" $4 ":" $5 "-" $6 ";QUERY_STRAND=" $12
        printf "%s\t%s\t%s\t%s\t%s\t.\t.\t%s\tGT\t1/0\n", $1, $2, id, $10, $11, info
    }
    ' "$TMPDIR/cigar_parsed.txt" > "$vcf_snv"
fi

# ==========================================================
# 4. Run python script for DEL / INS / INV  (once, if needed)
# ==========================================================
if ${PROCESS[del]} || ${PROCESS[ins]} || ${PROCESS[inv]}; then
    awk '$7 ~ /DEL|INS|INV/' "$input_vcf" \
    | sed 's/SDR_INV\|SV_INV/INV/g;
           s/SDR_INS\|SV_INS/INS/g;
           s/SDR_DEL\|SV_DEL/DEL/g' \
    > "$TMPDIR/sv_input.txt"

    if [ -s "$TMPDIR/sv_input.txt" ]; then
        python "$python_script" \
            --r "$ref_fasta" --q "$query_fasta" \
            --i "$TMPDIR/sv_input.txt" --o "$sdr_output"
    else
        cp "$cigar_in" "$sdr_output"
    fi
fi

# ---------- helper: format CIGAR INS/DEL rows into VCF ----------
# Usage: format_cigar_sv <svtype> <cigar_parsed> <outfile>
format_cigar_sv() {
    local svtype=$1 infile=$2 outfile=$3
    local svlen_prefix=""
    [ "$svtype" = "DEL" ] && svlen_prefix="-"
    awk -F'\t' -v T="$svtype" -v P="$svlen_prefix" '
    $9 == T {
        id   = $1 "-" $2 "-" T "-" $8
        info = "ID=" id ";SVTYPE=" T ";SVLEN=" P $8 ";TIG_REGION=" $4 ":" $5 "-" $6 ";QUERY_STRAND=" $12
        printf "%s\t%s\t%s\t%s\t%s\t.\t.\t%s\tGT\t1/0\n", $1, $2, id, $10, $11, info
    }
    ' "$infile" > "$outfile"
}

# ==========================================================
# 5. INS
# ==========================================================
vcf_ins="$TMPDIR/ins.vcf"
if ${PROCESS[ins]}; then
    log "INFO" "Processing insertions..."
    format_cigar_sv INS "$TMPDIR/cigar_parsed.txt" "$TMPDIR/cigar_ins.vcf"
    {
        awk 'index($3,"INS")' "$sdr_output"
        cat "$TMPDIR/cigar_ins.vcf"
    } > "$vcf_ins"
fi

# ==========================================================
# 6. DEL
# ==========================================================
vcf_del="$TMPDIR/del.vcf"
if ${PROCESS[del]}; then
    log "INFO" "Processing deletions..."
    format_cigar_sv DEL "$TMPDIR/cigar_parsed.txt" "$TMPDIR/cigar_del.vcf"
    {
        awk 'index($3,"DEL")' "$sdr_output"
        cat "$TMPDIR/cigar_del.vcf"
    } > "$vcf_del"
fi

# ==========================================================
# 7. INV
# ==========================================================
vcf_inv="$TMPDIR/inv.vcf"
if ${PROCESS[inv]}; then
    log "INFO" "Processing inversions..."
    awk 'index($3,"INV")' "$sdr_output" > "$vcf_inv"
fi

# ==========================================================
# 8. Merge VCF and prepend header
# ==========================================================
log "INFO" "Merging VCF records..."
{
    cat "$TMPDIR/vcf_header.txt"
    {
        ${PROCESS[snv]} && cat "$vcf_snv"  || true
        ${PROCESS[ins]} && cat "$vcf_ins"  || true
        ${PROCESS[del]} && cat "$vcf_del"  || true
        ${PROCESS[inv]} && cat "$vcf_inv"  || true
    } | awk '
    {
        if (!seen[$0]++) {
            print $0
        }
    }'
} > "$vcf_out"

# ==========================================================
# 9. Build BED file  (scan variants.vcf once per type with awk)
# ==========================================================
log "INFO" "Creating BED file..."
{
    echo -e "#CHROM\tPOS\tID\tSVTYPE\tSVLEN\tHAP\tGT\tQUERY"

    if ${PROCESS[snv]}; then
        awk -v hap="$hap" -F'\t' '
        /SNV/ {
            match($8, /TIG_REGION=([^;,]+)/, m)
            print $1, $2, $3, "SNV", 1, hap, "1|.", m[1]
        }' OFS='\t' "$vcf_out"
    fi

    if ${PROCESS[ins]}; then
        awk -v hap="$hap" -F'\t' '
        /\tINS\t|SVTYPE=INS/ {
            match($8, /SVLEN=([0-9]+)/, sl)
            match($8, /TIG_REGION=([^;,]+)/, m)
            print $1, $2, $3, "INS", sl[1], hap, "1|.", m[1]
        }' OFS='\t' "$vcf_out"
    fi

    if ${PROCESS[del]}; then
        awk -v hap="$hap" -F'\t' '
        /\tDEL\t|SVTYPE=DEL/ {
            match($8, /SVLEN=-?([0-9]+)/, sl)
            match($8, /TIG_REGION=([^;,]+)/, m)
            print $1, $2, $3, "DEL", sl[1], hap, "1|.", m[1]
        }' OFS='\t' "$vcf_out"
    fi

    if ${PROCESS[trans]}; then
        awk -v hap="$hap" 'OFS="\t"
        /TRANS/ {print $1,$2,$1"-"$2"-TRANS-"$9,"TRANS",$9,hap,"1|.",$4":"$5"-"$6}
        ' "$input_vcf"
    fi

    if ${PROCESS[sdr]}; then
        awk -v hap="$hap" 'OFS="\t"
        /SDR_NM|SV_NM/ {print $1,$2,$1"-"$2"-SDR-"$9,"SDR",$9,hap,"1|.",$4":"$5"-"$6}
        /INV-INV/      {print $1,$2,$1"-"$2"-SDR-"$9,"INV-INV",$9,hap,"1|.",$4":"$5"-"$6}
        ' "$input_vcf"
    fi

    if ${PROCESS[dup]}; then
        awk -v hap="$hap" 'OFS="\t"
        /DUP/ {print $1,$2,$1"-"$2"-DUP-"$9,"DUP",$9,hap,"1|.",$4":"$5"-"$6}
        ' "$input_vcf"
    fi

    if ${PROCESS[highdup]}; then
        awk -v hap="$hap" 'OFS="\t"
        /high-dup/ {print $1,$2,$1"-"$2"-HighDup","HighDup",".",hap,"1|.","."}
        ' "$input_vcf"
    fi

    if ${PROCESS[inv]}; then
        awk -v hap="$hap" 'OFS="\t"
        /SDR_INV|SV_INV/ {print $1,$2,$1"-"$2"-INV-"$9,"INV",$9,hap,"1|.",$4":"$5"-"$6}
        ' "$input_vcf"
    fi

    if ${PROCESS[complex]}; then
        awk -v hap="$hap" 'OFS="\t"
        /COMPLEX/ {print $1,$2,$1"-"$2"-SDR_COMPLEX-"$9,"SDR_COMPLEX",$9,hap,"1|.",$4":"$5"-"$6}
        ' "$input_vcf"
    fi

} | awk 'OFS="\t" { print $1,$2,$3,$4,$5,$6,$7,$8 }' \
  | sort | uniq | awk '$4=="SNV" || $4=="INS" || $4=="DEL" || $4=="INV" || $4=="DUP" || $4=="INV-INV" || $4=="HighDup" || $4=="TRANS" || $4=="SDR" || $4=="SDR_COMPLEX" || $4=="NA"'> "$bed_output"

# ==========================================================
# 10. Fix hap2 GT if present
# ==========================================================
hap2_bed="results/LGvarhap2.bed"
if [ -f "$hap2_bed" ]; then
    awk 'BEGIN{OFS="\t"} /^#/{print;next} NF>=8{$8=".|1"} {print}' \
        "$hap2_bed" > "$TMPDIR/hap2_fixed.bed"
    mv "$TMPDIR/hap2_fixed.bed" "$hap2_bed"
fi

log "INFO" "Done!"
