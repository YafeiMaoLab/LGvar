#!/bin/bash
usage() {
    echo "Usage: $0 -i input_file -r reference_fasta -q query_fasta -o output_file -t threads -c chunk_size -p paf_file -f fraction"
    echo "Options:"
    echo "  -i  Input TSV file"
    echo "  -r  Reference genome FASTA"
    echo "  -q  Query haplotype FASTA"
    echo "  -o  Output file"
    echo "  -t  Threads"
    echo "  -c  Chunk_size"
    echo "  -p  paf_file"
    echo "  -f  fraction"
    exit 1
}

while getopts "i:r:q:o:t:c:p:f:" opt; do
    case $opt in
        i) input="$OPTARG" ;;
        r) ref_path="$OPTARG" ;;
        q) hap_path="$OPTARG" ;;
        o) output="$OPTARG" ;;
        t) threads="$OPTARG" ;;
        c) chunk_size="$OPTARG" ;;
        p) paf_file="$OPTARG" ;;
        f) fraction="$OPTARG" ;;
        *) usage ;;
    esac
done

if [[ -z "$input" || -z "$ref_path" || -z "$hap_path" || -z "$output" || -z "$threads" || -z "$chunk_size" || -z "$paf_file" || -z "$fraction" ]]; then
    echo "Error: Missing required parameters!"
    usage
fi

tmpdir=$(mktemp -d -p . -t tmp.XXXXXXXXXX)
max_threads=$threads
chunk_size=$chunk_size
fraction=$fraction

mkdir -p half

process_paf_inversions() {
    local paf_file="$1"
    local output_file="$2"
    local new_rows_to_add=()
    
    while IFS=$'\t' read -ra paf_cols; do
        if [[ "${paf_cols[4]}" == "-" ]]; then
            local match_length="${paf_cols[9]}"
            
            if (( match_length < 100000 )); then
                local query_name="${paf_cols[0]}"
                local query_chr="${query_name%%:*}"
                local target_name="${paf_cols[5]}"
                local target_chr="${target_name%%:*}"
                
                local query_start="${paf_cols[2]}"
                local query_end="${paf_cols[3]}"
                local target_start="${paf_cols[7]}"
                local target_end="${paf_cols[8]}"
                
                local query_name_coords="${query_name#*:}"
                local query_name_start="${query_name_coords%-*}"
                local target_name_coords="${target_name#*:}"
                local target_name_start="${target_name_coords%-*}"
                
                local actual_target_start=$((target_name_start + target_start))
                local actual_target_end=$((target_name_start + target_end))
                local actual_query_start=$((query_name_start + query_start))
                local actual_query_end=$((query_name_start + query_end))
                
                if ((match_length < 10000)); then
                    anno="SV_INV"
                else
                    anno="SDR_INV"
                fi
                local new_row=(
                    "$target_chr"
                    "$actual_target_start"
                    "$actual_target_end"
                    "$query_chr"
                    "$actual_query_start"
                    "$actual_query_end"
                    "$anno"
                    "${paf_cols[4]}"
                    "${paf_cols[10]}"
                    "${paf_cols[9]}"
                )
                
                new_rows_to_add+=("$(IFS=$'\t'; echo "${new_row[*]}")")
            fi
        fi
    done < "$paf_file"
    
    if [[ ${#new_rows_to_add[@]} -gt 0 ]]; then
        for new_row in "${new_rows_to_add[@]}"; do
            printf "%s\n" "$new_row" >> "$output_file"
        done
    fi
}

export -f process_paf_inversions

process_chunk() {
    local chunk_file="$1"
    local chunk_tmpdir=$(mktemp -d -p "$tmpdir")
    local chunk_output="$chunk_tmpdir/result.tsv"
    local fraction=$fraction
    
    local local_paf_file="$chunk_tmpdir/aln.paf"
    local local_minimap2_log="$chunk_tmpdir/minimap2.log"
    
    while IFS=$'\t' read -r -a cols; do
        newline=("${cols[@]}")

        if [[ "${cols[0]}" == "ref_chr" ]]; then
            printf "%s\n" "$(IFS=$'\t'; echo "${newline[*]}")" >> "$chunk_output"
            continue
        fi

        # SDR_NM/SV_NM
        if [[ "${cols[6]}" == "SDR_NM" || "${cols[6]}" == "SV_NM" ]]; then
            ref_chr="${cols[0]}"
            ref_start="${cols[1]}"
            ref_end="${cols[2]}"
            q_chr="${cols[3]}"
            q_start="${cols[4]}"
            q_end="${cols[5]}"
            ref_len="${cols[8]}"
            q_len="${cols[9]}"

            len=$(( ref_len < q_len ? ref_len : q_len ))
            max_len=$(( ref_len > q_len ? ref_len : q_len))

            if (( max_len >= 1000000 )); then
                printf "%s\n" "$(IFS=$'\t'; echo "${newline[*]}")" >> "$chunk_output"
                continue
            fi
            # inversion re-identification
            ref_fa="$chunk_tmpdir/ref_${ref_chr}_${ref_start}-${ref_end}.fa"
            query_fa="$chunk_tmpdir/query_${q_chr}_${q_start}-${q_end}.fa"
            
            if samtools faidx "$ref_path" "${ref_chr}:${ref_start}-${ref_end}" > "$ref_fa" 2>/dev/null &&
               samtools faidx "$hap_path" "${q_chr}:${q_start}-${q_end}" > "$query_fa" 2>/dev/null
            then
                minimap2 -t 24 -cx asm20 --eqx --secondary=no "$ref_fa" "$query_fa" > "$local_paf_file" 2>> "$local_minimap2_log"
                
                if [[ -s "$local_paf_file" ]]; then
                    total_matches=0
                    has_negative=false
                    has_positive=false
                    all_negative=true

                    while IFS=$'\t' read -ra paf_cols; do
                        if [[ "${paf_cols[4]}" == "-" ]]; then
                            has_negative=true
                        else
                            has_positive=true
                            all_negative=false
                        fi
                        ((total_matches += paf_cols[9]))  
                    done < "$local_paf_file"

                    if [[ "$has_negative" == "true" ]] && [[ "$all_negative" == "false" ]] && (( total_matches < max_len * fraction / 10)) && ((max_len < 500000)); then
                        process_paf_inversions "$local_paf_file" "$chunk_output"
                    fi

                    if [[ "$all_negative" == "true" ]] && (( total_matches > max_len * fraction / 10)); then
                        if [[ "${cols[6]}" == "SDR_NM" ]]; then
                            newline[6]="SDR_INV"
                        elif [[ "${cols[6]}" == "SV_NM" ]]; then
                            newline[6]="SV_INV"
                        fi
                        newline[7]="-"
                        printf "%s\n" "$(IFS=$'\t'; echo "${newline[*]}")" >> "$chunk_output"
                    elif ((max_len < 500000)); then
                        process_paf_inversions "$local_paf_file" "$chunk_output"
                    fi

                    cat "$local_paf_file" >> "$chunk_tmpdir/merged.paf"
                fi
            else
                echo "ERROR: Failed to extract sequences for ${ref_chr}:${ref_start}-${ref_end}" >&2
            fi

        # SDR_INS 
        elif [[ "${cols[6]}" == "SDR_INS" ]]; then
            ref_chr="${cols[0]}"
            ref_start="${cols[1]}"
            local ins_result="$chunk_tmpdir/ins.result"

            awk -v chr="$ref_chr" '$6 == chr' "$paf_file" | cut -f6,8-9 | sort -k2,2n |
            bedtools merge -i - -d 500 |  
            bedtools closest -a <(echo -e "$ref_chr\t$ref_start\t$((ref_start+1))") -b - -D a |
            awk -v OFS="\t" '{
                dist_start = ($3 - $6) < 0 ? ($6 - $3) : ($3 - $6);
                dist_end = ($3 - $7) < 0 ? ($7 - $3) : ($3 - $7);
                min_dist = (dist_start < dist_end) ? dist_start : dist_end;
                
                if (min_dist <= 10) print "YES";
                else print "NO";
            }' > "$ins_result"
            
            if [[ -s "$ins_result" ]] && grep -q "YES" "$ins_result"; then
                  printf "%s\n" "$(IFS=$'\t'; echo "${newline[*]}")" >> "$chunk_output"
                  continue
            else
                  continue
            fi
        fi

        printf "%s\n" "$(IFS=$'\t'; echo "${newline[*]}")" >> "$chunk_output"
        
    done < "$chunk_file"

    cat "$chunk_output"
}

export -f process_chunk
export ref_path hap_path tmpdir paf_file fraction

> half/minimap2.paf

split -l $chunk_size --numeric-suffixes --additional-suffix=".tsv" "$input" "$tmpdir/chunk_"

find "$tmpdir" -name "chunk_*.tsv" | sort | \
    xargs -P $max_threads -I {} bash -c 'process_chunk "{}"' > "$output.tmp"

{
    head -n 1 "$input"
    grep -v "^ref_chr" "$output.tmp"
} > "$output"

find "$tmpdir" -name "merged.paf" -exec cat {} + > half/minimap2.paf 2>/dev/null || true

rm -rf "$tmpdir" "$output.tmp"
