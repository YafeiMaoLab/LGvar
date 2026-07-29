import argparse
import os
import pandas as pd
from Bio import SeqIO

parser = argparse.ArgumentParser(description='tsv2vcf')
parser.add_argument('--r', required=True, help='Reference fasta')
parser.add_argument('--q', required=True, help='Query fasta')
parser.add_argument('--i', required=True, help='SDR result file')
parser.add_argument('--o', required=True, help='Output VCF-like file')
args = parser.parse_args()

ref_genome = {record.id: record for record in SeqIO.parse(args.r, 'fasta')}
query_genome = {record.id: record for record in SeqIO.parse(args.q, 'fasta')}

rev_comp_table = str.maketrans("ATCGatcgNn", "TAGCtagcNn")

def get_sequence(genome_dict, chrom, start, end, strand="+"):
    if chrom not in genome_dict:
        return "."
    start_pos = max(0, start)
    seq_obj = genome_dict[chrom].seq[start_pos:end]

    if strand == "-":
        return str(seq_obj.reverse_complement())
    return str(seq_obj)

df = pd.read_csv(args.i, sep='\t', header=None)

with open(args.o, 'w') as file:
    file.write("##fileformat=VCFv4.2\n")
    file.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n")
    
    for index, row in df.iterrows():
        ref_chr, ref_start, ref_end = row[0], int(row[1]), int(row[2])
        que_chr, que_start, que_end = row[3], int(row[4]), int(row[5])
        sv_type = row[6]      
        strand = row[7]      
        
        lenall = int(row[9]) if sv_type == 'INS' else int(row[8])
        
        if ref_start == 0 or que_start == 0:
            anchor_ref_start = 0
            anchor_ref_end = max(ref_end, 1)
            anchor_que_start = 0
            anchor_que_end = max(que_end, 1)
            vcf_pos = 1 # VCF 是 1-base
        else:
            anchor_ref_start = ref_start - 1
            anchor_ref_end = ref_end
            anchor_que_start = que_start - 1
            anchor_que_end = que_end
            vcf_pos = ref_start 
            
        actual_strand = "-" if sv_type == 'INV' else strand
        
        selected_reffa = get_sequence(ref_genome, ref_chr, anchor_ref_start, anchor_ref_end, "+")
        selected_quefa = get_sequence(query_genome, que_chr, anchor_que_start, anchor_que_end, actual_strand)

        var_id = f"{ref_chr}-{vcf_pos}-{sv_type}-{lenall}"
        info_field = f"ID={var_id};SVTYPE={sv_type};SVLEN={lenall};TIG_REGION={que_chr}:{que_start}-{que_end};QUERY_STRAND=+,+"

        file.write(f"{ref_chr}\t{vcf_pos}\t{var_id}\t{selected_reffa}\t{selected_quefa}\t.\t.\t{info_field}\tGT\t1|0\n")
