#!/bin/bash

if [ -z "$1" ]; then
  echo "Usage: $0 <paf_file>"
  exit 1
fi

log() {
    local level=$1
    local message=$2
    local timestamp=$(date "+%Y-%m-%d %H:%M:%S,%3N")

    echo -e "${timestamp} - ${level} - ${message}"
}

paf_file="$1"
length="$2"
output="$3"

awk '
{
  query_chr = $1;
  ref_chr = $6;
  
  key = ref_chr ":" query_chr;
  count[key]++;
  
  query_length = $4 - $3;
  ref_length = $9 - $8;
  
  total_query_length[key] += query_length;
  total_ref_length[key] += ref_length;
  #total_length[key] += (query_length + ref_length);
  
  if (!(ref_chr in ref_lengths)) {
    ref_lengths[ref_chr] = $7;  # reflen
  }
  if (!(query_chr in query_lengths)) {
    query_lengths[query_chr] = $2;  # querylen
  }
}
END {
  for (ref_chr in ref_lengths) {
    printf "Reference Chromosome: %s (Length: %d)\n", ref_chr, ref_lengths[ref_chr];
    for (key in count) {
      split(key, arr, ":");
      current_ref_chr = arr[1];
      current_query_chr = arr[2];
      
      if (current_ref_chr == ref_chr) {
        printf "  Query Chromosome: %s (Length: %d)\n", current_query_chr, query_lengths[current_query_chr];
        printf "    Total Alignment Records: %d\n", count[key];
        printf "    Total Query Length: %d\n", total_query_length[key];
        printf "    Total Ref Length: %d\n", total_ref_length[key];
      }
    }
    printf "----------------------------------------\n";
  }
}
' "$paf_file" > "$output".txt

declare -A ref_to_queries
current_ref=""
current_query=""
total_query_length=0
total_ref_length=0

while IFS= read -r line; do
    if [[ $line =~ ^Reference\ Chromosome:\ ([^ ]+) ]]; then
        current_ref=${BASH_REMATCH[1]}
        current_query=""
        total_query_length=0
        total_ref_length=0
    elif [[ $line =~ ^\ *Query\ Chromosome:\ ([^ ]+) ]]; then
        current_query=${BASH_REMATCH[1]}
    elif [[ $line =~ ^\ *Total\ Query\ Length:\ ([0-9]+) ]]; then
        total_query_length=${BASH_REMATCH[1]}
    elif [[ $line =~ ^\ *Total\ Ref\ Length:\ ([0-9]+) ]]; then
        total_ref_length=${BASH_REMATCH[1]}
        
        if [ $total_query_length -gt "$length" ] && [ $total_ref_length -gt "$length" ]; then
            if [ -n "$current_ref" ] && [ -n "$current_query" ]; then
                if [ -z "${ref_to_queries[$current_ref]}" ]; then
                    ref_to_queries[$current_ref]=$current_query
                else
                    if [[ ! "${ref_to_queries[$current_ref]}" =~ $current_query ]]; then
                        ref_to_queries[$current_ref]="${ref_to_queries[$current_ref]},$current_query"
                    fi
                fi
            fi
        fi
    fi
done < "$output".txt

echo -e "ref\tquery" > "$output".tsv
for ref in "${!ref_to_queries[@]}"; do
    echo -e "$ref\t${ref_to_queries[$ref]}" >> "$output".tsv
done
rm "$output".txt
log "INFO" "Done! Pairs in $output.tsv"

