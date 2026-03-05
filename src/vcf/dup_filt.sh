#!/bin/bash

if [ $# -ne 4 ]; then
    echo "Usage: $0 <paf_file> <input_file> <output_file> <script.py>"
    exit 1
fi

paf="$1"
file="$2"
out="$3"
script="$4"

mkdir -p temp

for f in "$paf" "$file" "$script"; do
    if [ ! -f "$f" ]; then
        echo "Error: File $f not found." >&2
        exit 1
    fi
done

sort -k8,8n "$paf" | awk '{OFS=FS="\t"}{print $6,$8,$9}' > temp/refpaf.region

awk '{count[$0]++} END {for (line in count) if (count[line] > 1) print line}' temp/refpaf.region > temp/duplicates.txt

if [ -s temp/refpaf.region ]; then
    bedtools intersect -a temp/refpaf.region -b temp/refpaf.region 2>/dev/null | \
        awk 'NR==FNR{a[$0];next} !($0 in a)' temp/refpaf.region - | \
        sort -k1,1V -k2,2n -k3,3n | uniq | bedtools merge -i - > temp/intersect.txt
else
    > temp/intersect.txt
fi

touch temp/duplicates.txt temp/intersect.txt
cat temp/duplicates.txt temp/intersect.txt | sort -k1,1V -k2,2n -k3,3n > temp/mergedup.txt

awk -F'[:-]' 'BEGIN {OFS="\t"}{print $0,$1,$2}' "$file" | \
    awk '{OFS=FS="\t";print $0,$10+$2}' | nl > temp/cigartestin.txt

awk -v OFS='\t' 'NR>1 {print $10, $11, $12, $1, $2}' temp/cigartestin.txt > temp/cigartestin.bed

if [ -s temp/mergedup.txt ] && [ -s temp/cigartestin.bed ]; then
    bedtools intersect -a temp/mergedup.txt -b temp/cigartestin.bed -wa -wb > temp/overlaps.txt
else
    > temp/overlaps.txt
fi

if [ -s temp/overlaps.txt ]; then
    awk 'BEGIN {OFS="\t"} {if ($6 >= $2 && $6 <= $3) print}' temp/overlaps.txt > temp/overlaps_filt.txt
else
    > temp/overlaps_filt.txt
fi

if [ -s temp/overlaps_filt.txt ]; then
    python "$script" temp/overlaps_filt.txt fp
    shopt -s nullglob
    del_files=(temp/*_del.txt)
    shopt -u nullglob
    if [ ${#del_files[@]} -gt 0 ]; then
        cat "${del_files[@]}" > temp/del.txt
    else
        > temp/del.txt
    fi

    if [ -s temp/del.txt ]; then
        awk 'NR==FNR{d[$1]; next} !($1 in d)' temp/del.txt temp/cigartestin.txt | cut -f 2-9 > "$out"
    else
        cut -f 2-9 temp/cigartestin.txt > "$out"
    fi
else
    cut -f 2-9 temp/cigartestin.txt > "$out"
fi

find temp/ -type f \( -name "refpaf.region" \
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
