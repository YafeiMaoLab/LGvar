#!/bin/bash

if [ -z "$1" ] || [ -z "$2" ] || [ -z "$3" ]; then
  echo "Usage: $0 <paf_file> <length_threshold> <output_prefix>"
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

log "INFO" "Processing PAF file with length threshold: $length"

awk -v len_threshold="$length" '
{
  query_chr = $1;
  ref_chr = $6;
  
  key = ref_chr ":" query_chr;
  
  query_length = $4 - $3;
  ref_length = $9 - $8;
  
  total_query_length[key] += query_length;
  total_ref_length[key] += ref_length;

  ref_list[ref_chr] = 1;
}
END {
  print "ref\tquery"; 
  
  for (ref in ref_list) {
    out_str = ""; 
    
    for (key in total_query_length) {
      split(key, arr, ":");
      curr_ref = arr[1];
      curr_query = arr[2];
      
      if (curr_ref == ref) {
        if (total_query_length[key] > len_threshold && total_ref_length[key] > len_threshold) {
          if (out_str == "") {
            out_str = curr_query;
          } else {
            out_str = out_str "," curr_query;
          }
        }
      }
    }

    if (out_str != "") {
      print ref "\t" out_str;
    }
  }
}
' "$paf_file" > "${output}.tsv"

log "INFO" "Done! Pairs saved in ${output}.tsv"
