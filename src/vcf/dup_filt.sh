#!/bin/bash
paf=$1
file=$2
out=$3
script=$4

sort -k8,8n "$paf" |awk '{OFS=FS="\t"}{print $6,$8,$9}' > temp/refpaf.region
awk '{count[$0]++} END {for (line in count) if (count[line] > 1) print line}' temp/refpaf.region > temp/duplicates.txt
bedtools intersect -a temp/refpaf.region -b temp/refpaf.region | awk 'NR==FNR{a[$0];next} !($0 in a)' temp/refpaf.region - | sort -k1,1V -k2,2n -k3,3n | uniq | bedtools merge -i - > temp/intersect.txt
cat temp/duplicates.txt temp/intersect.txt | sort -k1,1V -k2,2n -k3,3n > temp/mergedup.txt

awk -F'[:-]' 'BEGIN {OFS="\t"}{print $0,$1,$2}' $file | awk '{OFS=FS="\t";print $0,$10+$2}' | nl > temp/cigartestin.txt

awk -v OFS='\t' 'NR>1 {print $10, $11, $12, $1, $2}' temp/cigartestin.txt > temp/cigartestin.bed

bedtools intersect -a temp/mergedup.txt -b temp/cigartestin.bed -wa -wb > temp/overlaps.txt
awk 'BEGIN {OFS="\t"} {if ($6 >= $2 && $6 <= $3) {print}}' temp/overlaps.txt > temp/overlaps_filt.txt

if [[ -s temp/overlaps_filt.txt ]]; then
    python "$script" temp/overlaps_filt.txt fp
    cat temp/*_del.txt > temp/del.txt
    
    if [ -s temp/del.txt ]; then
        awk 'NR==FNR{d[$1]; next} !($1 in d)' temp/del.txt temp/cigartestin.txt | 
        cut -f 2-9 > "$out"
    else
        cut -f 2-9 temp/cigartestin.txt > "$out"
    fi
else
    cut -f 2-9 temp/cigartestin.txt > "$out"
fi

find . -maxdepth 1 -type f \( -name "refpaf.region" \
                             -o -name "mergedup.txt" \
                             -o -name "cigartestin.txt" \
                             -o -name "cigartestin.bed" \
                             -o -name "overlaps.txt" \
                             -o -name "overlaps_filt.txt" \
                             -o -name "del.txt" \
                             -o -name "duplicates.txt" \
                             -o -name "intersect.txt" \
                             -o -name "*_del.txt" \
                             -o -name "*_counts.txt" \) -delete
