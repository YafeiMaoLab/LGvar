import sys
import os
import re
import subprocess

def ensure_fai(fa_path):
    if not fa_path or not os.path.exists(fa_path):
        print(f"Error: Can't find FASTA: {fa_path}")
        sys.exit(1)
        
    fai_path = fa_path + ".fai"
    is_temporary = False
    
    if not os.path.exists(fai_path):
        try:
            subprocess.run(["samtools", "faidx", fa_path], check=True)
            is_temporary = True
        except (subprocess.CalledProcessError, FileNotFoundError) as e:
            sys.exit(1)
            
    return fai_path, is_temporary

def load_fai(fai_path):
    len_dict = {}
    with open(fai_path, 'r') as f:
        for line in f:
            if line.strip():
                fields = line.strip().split('\t')
                len_dict[fields[0]] = int(fields[1])
    return len_dict

def batch_restore_paf(input_paf, output_paf, query_fa, target_fa):
    q_fai, q_is_tmp = ensure_fai(query_fa)
    t_fai, t_is_tmp = ensure_fai(target_fa)
    
    try:
        q_lengths = load_fai(q_fai)
        t_lengths = load_fai(t_fai)

        with open(input_paf, 'r') as infile, open(output_paf, 'w') as outfile:
            for line in infile:
                if not line.strip():
                    continue
                fields = line.strip().split('\t')
                
                q_name_raw = fields[0]
                q_match = re.match(r"([^:]+):(\d+)-(\d+)", q_name_raw)
                
                if q_match:
                    q_chrom, q_frag_start, _ = q_match.group(1), int(q_match.group(2)), int(q_match.group(3))
                    fields[0] = q_chrom
                    fields[1] = str(q_lengths.get(q_chrom, fields[1]))
                    fields[2] = str(q_frag_start - 1 + int(fields[2]))
                    fields[3] = str(q_frag_start - 1 + int(fields[3]))
                
                t_name_raw = fields[5]
                t_match = re.match(r"([^:]+):(\d+)-(\d+)", t_name_raw)
                
                if t_match:
                    t_chrom, t_frag_start, _ = t_match.group(1), int(t_match.group(2)), int(t_match.group(3))
                    fields[5] = t_chrom
                    fields[6] = str(t_lengths.get(t_chrom, fields[6]))
                    fields[7] = str(t_frag_start - 1 + int(fields[7]))
                    fields[8] = str(t_frag_start - 1 + int(fields[8]))
                
                outfile.write("\t".join(fields) + "\n")
        
    finally:
        if q_is_tmp and os.path.exists(q_fai):
            os.remove(q_fai)
        if t_is_tmp and os.path.exists(t_fai):
            os.remove(t_fai)

if __name__ == "__main__":
    if len(sys.argv) < 5:
        print("Usage: python restore_paf.py <input.paf> <output.paf> <query.fa> <target.fa>")
        sys.exit(1)
        
    batch_restore_paf(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
