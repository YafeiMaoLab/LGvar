## --*-- coding=iso-8859-1 --*--
import gzip
import re
import time
import numpy as np
from scipy.cluster.hierarchy import linkage, fcluster
from collections import defaultdict
import edlib
import argparse
import logging
import os
import copy
from collections import Counter
from tqdm import tqdm
import multiprocessing
from functools import partial

# logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

class Variant:
    def __init__(self, chrom, pos, ref, alt, var_type, end=None, svlen=None, svtype=None, info=None):
        self.chrom = chrom
        self.pos = pos
        self.ref = ref
        self.alt = alt
        self.var_type = var_type  # 'SNV', 'INDEL', 'SV'
        self.end = end if end else pos + len(ref) - 1
        self.svlen = svlen
        self.svtype = svtype
        self.info = info if info else {}
        self.haplotype = None
        self.id = None
        self.qual = None
        self.filter = None
        self.format = None
        self.sample = None
        self.genotype = None
    
    def get_key(self):
        return (self.chrom, self.var_type, self.pos)
    
    def __repr__(self):
        return f"{self.chrom}:{self.pos}-{self.end} {self.var_type} {self.ref}/{self.alt}"

def parse_vcf(vcf_file, chrom=None, file_type=None):
    """VCF"""
    variants = []
    opener = gzip.open if vcf_file.endswith('.gz') else open
    mode = 'rt' if vcf_file.endswith('.gz') else 'r'
    
    with opener(vcf_file, mode) as f:
        for line in f:
            if line.startswith('#'):
                continue
            
            fields = line.strip().split('\t')
            if len(fields) < 5:
                continue
            
            chrom_field, pos, var_id, ref, alt = fields[:5]
            
            # specific chromosome
            if chrom and chrom_field != chrom:
                continue
            
            qual = fields[5] if len(fields) > 5 else '.'
            filt = fields[6] if len(fields) > 6 else '.'
            info = fields[7] if len(fields) > 7 else '.'
            fmt = fields[8] if len(fields) > 8 else '.'
            sample = fields[9] if len(fields) > 9 else '.'
            
            info_dict = {}
            for item in info.split(';'):
                if '=' in item:
                    key, val = item.split('=', 1)
                    info_dict[key] = val
                else:
                    info_dict[item] = True
            
            svtype = info_dict.get('SVTYPE', None)

            if file_type:
                var_type = file_type.upper()

            elif svtype:
                if svtype == 'SNV':
                    var_type = 'SNV'
                elif svtype in ['INS', 'DEL', 'INV', 'DUP']:
                    var_type = 'SV'
                else:
                    var_type = 'SV'  

            else:
                if len(ref) == 1 and len(alt) == 1:
                    var_type = 'SNV'
                else:
                    var_type = 'INDEL'       

            end = None
            svlen = None
            if var_type == 'SNV':
                end = int(pos) + len(ref) - 1
                svlen = abs(len(alt) - len(ref))
            else:
                if 'END' in info_dict:
                    end = int(info_dict['END'])
                elif 'SVLEN' in info_dict:
                    svlen_val = info_dict['SVLEN']
                    if svlen_val.startswith('-'):
                        svlen = abs(int(svlen_val[1:]))
                        end = int(pos) + svlen
                    else:
                        svlen = int(svlen_val)
                        end = int(pos) + svlen
                else:
                    end = int(pos)
                    svlen = 0
            
            # genotype
            genotype = None
            if fmt and sample:
                fmt_fields = fmt.split(':')
                sample_fields = sample.split(':')
                if 'GT' in fmt_fields:
                    gt_index = fmt_fields.index('GT')
                    genotype = sample_fields[gt_index]
            
            variant = Variant(
                chrom=chrom_field,
                pos=int(pos),
                ref=ref,
                alt=alt,
                var_type=var_type,
                end=end,
                svlen=svlen,
                svtype=svtype,
                info=info_dict
            )
            
            variant.id = var_id
            variant.qual = qual
            variant.filter = filt
            variant.format = fmt
            variant.sample = sample
            variant.genotype = genotype
            
            variants.append(variant)
    
    return variants

def form_clusters(variants, max_distance, small_distance):
    sorted_variants = sorted(variants, key=lambda v: v.get_key())
    clusters = []
    current_cluster = []
    
    for variant in sorted_variants:
        if not current_cluster:
            current_cluster.append(variant)
            continue
        
        last_variant = current_cluster[-1]
        last_key = last_variant.get_key()
        current_key = variant.get_key()

        distance = max_distance if variant.svtype == "SV" else small_distance
        
        if (last_key[0] != current_key[0] or 
            last_key[1] != current_key[1] or 
            abs(last_key[2] - current_key[2]) > distance):
            clusters.append(current_cluster)
            current_cluster = [variant]
        else:
            current_cluster.append(variant)
    
    if current_cluster:
        clusters.append(current_cluster)
    
    return clusters

def compute_variant_similarity(v1, v2, max_distance=500, small_distance=10):
    """Calculate similarity between two variants (1.0=identical, 0.0=completely different)."""
    if v1.chrom != v2.chrom or v1.var_type != v2.var_type:
        return 0.0
    
    pos_diff = abs(v1.pos - v2.pos)
    
    if v1.var_type in ['SNV', 'INDEL']:
        if pos_diff > small_distance:
            return 0.0
    else:
        if pos_diff > max_distance:
            return 0.0
            
    # seq similarity
    seq_sim = 1.0 
    
    if v1.var_type == 'SNV':
        ref_match = v1.ref == v2.ref
        alt_match = v1.alt == v2.alt
        seq_sim = 1.0 if ref_match and alt_match else 0.0
    
    elif v1.var_type == 'INDEL':
        seq1 = v1.alt if len(v1.ref) == 1 else v1.ref
        seq2 = v2.alt if len(v2.ref) == 1 else v2.ref
        if seq1 and seq2:
            align = edlib.align(seq1, seq2, task='path')
            max_len = max(len(seq1), len(seq2))
            seq_sim = 1 - align['editDistance'] / max_len if max_len > 0 else 0.0
    
    elif v1.var_type == 'SV':
        svtype = getattr(v1, 'svtype', None)
        def kmer_similarity(s1, s2, k=21):
            kmers1 = Counter(s1[i:i+k] for i in range(len(s1)-k+1))
            kmers2 = Counter(s2[i:i+k] for i in range(len(s2)-k+1))
            intersection = sum((kmers1 & kmers2).values())
            union = sum((kmers1 | kmers2).values())
            return intersection / union
        
        if svtype == 'INS':
            if pos_diff > max_distance:
                return 0.0
            else:
                seq1 = re.sub(r'^<[^>]+>', '', v1.alt)  
                seq2 = re.sub(r'^<[^>]+>', '', v2.alt)
                if seq1 and seq2:
                    seq_sim = kmer_similarity(seq1, seq2)
        
        elif svtype == 'DEL':
            if pos_diff > max_distance:
                return 0.0
            else:
                seq1 = re.sub(r'^<[^>]+>', '', v1.ref)  
                seq2 = re.sub(r'^<[^>]+>', '', v2.ref)
                if seq1 and seq2:
                    seq_sim = kmer_similarity(seq1, seq2)
        
        elif svtype == 'INV':
            # reciprocal ovelap >= 80%
            overlap = min(v1.end, v2.end) - max(v1.pos, v2.pos)
            min_len = max(v1.end-v1.pos, v2.end-v2.pos)
            seq_sim = 1.0 if (overlap >= 0.8*min_len) else 0.0
    
    return seq_sim  

def cluster_variants(variants, similarity_threshold=0.8, max_distance=500, small_distance=10):
    """Cluster variants based on similarity threshold."""
    if len(variants) <= 1:
        return [variants]

    def similarity_to_distance(v1, v2):
        return 1 - compute_variant_similarity(v1, v2, max_distance, small_distance)

    n = len(variants)
    dist_matrix = np.zeros((n, n))
    for i in range(n):
        for j in range(i+1, n):
            dist_matrix[i,j] = dist_matrix[j,i] = similarity_to_distance(variants[i], variants[j])
    
    # cluster
    Z = linkage(dist_matrix[np.triu_indices(n, k=1)], method='complete')
    clusters = fcluster(Z, 1 - similarity_threshold, criterion='distance')
    
    cluster_dict = defaultdict(list)
    for idx, cluster_id in enumerate(clusters):
        cluster_dict[cluster_id].append(variants[idx])
    
    return list(cluster_dict.values())

def merge_info_fields(cluster, base_info, fields_to_merge=['TIG_REGION', 'QUERY_STRAND']):
    for field in fields_to_merge:
        combined_values = []
        for v in cluster:
            if field in v.info and v.info[field]:
                parts = str(v.info[field]).split(',')
                combined_values.extend(parts)
        
        if combined_values:
            unique_values = list(dict.fromkeys([x.strip() for x in combined_values if x.strip()]))
            base_info[field] = ",".join(unique_values)
    return base_info

def merge_variants(variants_hap1, variants_hap2, max_distance, small_distance, similarity_threshold):
    for v in variants_hap1:
        v.haplotype = 1
    for v in variants_hap2:
        v.haplotype = 2
    
    all_variants = variants_hap1 + variants_hap2
    partitions = form_clusters(all_variants, max_distance, small_distance)
    merged_variants = []
    
    for partition in partitions:
        clusters = cluster_variants(partition, similarity_threshold, max_distance, small_distance)
        
        for cluster in clusters:
            if len(cluster) == 0:
                continue
            merged_id = cluster[0].id
            base_info = copy.deepcopy(cluster[0].info)

            base_info = merge_info_fields(cluster, base_info)

            min_pos = min(v.pos for v in cluster)
            max_end = max(v.end for v in cluster)

            if len(cluster) == 1:
                variant = cluster[0]
                variant.genotype = '1|0' if variant.haplotype == 1 else '0|1'

                variant.id = merged_id
                variant.info = base_info 
                variant.pos = min_pos
                variant.end = max_end
                if variant.var_type == 'SV' and 'END' in base_info:
                    base_info['END'] = str(max_end)
                merged_variants.append(variant)
            
            elif len(cluster) == 2:
                v1, v2 = cluster
                qual1 = float(v1.qual) if v1.qual != '.' else 0
                qual2 = float(v2.qual) if v2.qual != '.' else 0

                variant = copy.deepcopy(v1 if qual1 >= qual2 else v2)
                
                if v1.haplotype == v2.haplotype:
                    variant.genotype = '1|0' if v1.haplotype == 1 else '0|1'
                else:
                    variant.genotype = '1|1'

                variant.id = merged_id
                variant.info = base_info  
                variant.pos = min_pos
                variant.end = max_end
                if variant.var_type == 'SV' and 'END' in base_info:
                    base_info['END'] = str(max_end)
                merged_variants.append(variant)

            else:
                hap1_vars = [v for v in cluster if v.haplotype == 1]
                hap2_vars = [v for v in cluster if v.haplotype == 2]

                best_hap1 = max(hap1_vars, key=lambda v: float(v.qual) if v.qual != '.' else 0) if hap1_vars else None
                best_hap2 = max(hap2_vars, key=lambda v: float(v.qual) if v.qual != '.' else 0) if hap2_vars else None

                if best_hap1 and best_hap2:
                    variant = copy.deepcopy(best_hap1 if (float(best_hap1.qual) if best_hap1.qual != '.' else 0) >= 
                                           (float(best_hap2.qual) if best_hap2.qual != '.' else 0) else best_hap2)
                    variant.genotype = '1|1'
                else:
                    best_var = best_hap1 if best_hap1 else best_hap2
                    variant = copy.deepcopy(best_var)
                    variant.genotype = '1|0' if best_hap1 else '0|1'

                variant.id = merged_id
                variant.info = base_info  
                variant.pos = min_pos
                variant.end = max_end
                if variant.var_type == 'SV' and 'END' in base_info:
                    base_info['END'] = str(max_end)
                merged_variants.append(variant)
    
    return merged_variants

def get_header_file(args):
    input_files = [
        args.hap1_sv, args.hap2_sv,
        args.hap1_indel, args.hap2_indel,
        args.hap1_snv, args.hap2_snv
    ]
    for file in input_files:
        if file and os.path.exists(file):
            return file
    return None

def parse_source_header(vcf_file):
    """header"""
    headers = []
    opener = gzip.open if vcf_file.endswith('.gz') else open
    with opener(vcf_file, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                headers.append(line.strip())
            else:
                break
    return headers

def write_vcf_header(output_handle, args, sample_name="SAMPLE"):
    output_handle.write("##fileformat=VCFv4.2\n")
    output_handle.write(f"##fileDate={time.strftime('%Y%m%d')}\n")
    output_handle.write("##source=LGvar\n")
    
    header_file = get_header_file(args)
    if header_file:
        vcf_header = parse_source_header(header_file)
        for line in vcf_header:
            if line.startswith(('##reference', '##contig')):
                output_handle.write(line + "\n")
    
    output_handle.write("##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">\n")
    output_handle.write("##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Length of the variant\">\n")
    output_handle.write("##INFO=<ID=TIG_REGION,Number=.,Type=String,Description=\"Contig region where variant was found (one per alt with h1 before h2 for homozygous calls)\">\n")
    output_handle.write("##INFO=<ID=QUERY_STRAND,Number=.,Type=String,Description=\"Strand of variant in the contig relative to the reference (order follows TIG_REGION)\">\n")
    output_handle.write("##INFO=<ID=ID,Number=1,Type=String,Description=\"Original variant ID\">\n")
    output_handle.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
    output_handle.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + sample_name + "\n")

def write_variant(variant, output_handle):
    info_fields = []
    for key, value in variant.info.items():
        if value is True:
            info_fields.append(key)
        else:
            info_fields.append(f"{key}={value}")

    if variant.id and variant.id != '.':
        info_fields.append(f"ID={variant.id}")
    
    info_str = ";".join(info_fields)
    qual = variant.qual if variant.qual else "."
    filt = "PASS"
    fmt = "GT"
    sample_gt = variant.genotype if variant.genotype else ".|."

    fields = [
        variant.chrom,
        str(variant.pos),
        variant.id if variant.id else ".",
        variant.ref,
        variant.alt,
        qual,
        filt,
        info_str,
        fmt,
        sample_gt
    ]
    output_handle.write("\t".join(fields) + "\n")

def get_all_chromosomes(args):
    chromosomes = set()
    input_files = [
        args.hap1_sv, args.hap2_sv,
        args.hap1_indel, args.hap2_indel,
        args.hap1_snv, args.hap2_snv
    ]
    
    for vcf_file in input_files:
        if not vcf_file or not os.path.exists(vcf_file):
            continue
            
        opener = gzip.open if vcf_file.endswith('.gz') else open
        mode = 'rt' if vcf_file.endswith('.gz') else 'r'
        
        with opener(vcf_file, mode) as f:
            for line in f:
                if line.startswith('#'):
                    continue
                chrom = line.split('\t')[0]
                chromosomes.add(chrom)
    
    def chrom_key(c):
        try:
            if c.startswith('chr'):
                num = c[3:]
            else:
                num = c
            return (0, int(num)) if num.isdigit() else (1, c)
        except:
            return (2, c)
    
    return sorted(chromosomes, key=chrom_key)

def process_chromosome(chrom, args):
    """Variants of each chromosome"""
    logging.info(f"Processing chromosome: {chrom}")
    
    snv_hap1 = parse_vcf(args.hap1_snv, chrom, 'snv') if args.hap1_snv else []
    snv_hap2 = parse_vcf(args.hap2_snv, chrom, 'snv') if args.hap2_snv else []
    indel_hap1 = parse_vcf(args.hap1_indel, chrom, 'indel') if args.hap1_indel else []
    indel_hap2 = parse_vcf(args.hap2_indel, chrom, 'indel') if args.hap2_indel else []
    sv_hap1 = parse_vcf(args.hap1_sv, chrom, 'sv') if args.hap1_sv else []
    sv_hap2 = parse_vcf(args.hap2_sv, chrom, 'sv') if args.hap2_sv else []
    
    # concate
    merged_snv = merge_variants(snv_hap1, snv_hap2, 
                               max_distance=args.max_distance,
                               small_distance=args.small_distance,
                               similarity_threshold=args.similarity_threshold) if snv_hap1 or snv_hap2 else []
    
    merged_indel = merge_variants(indel_hap1, indel_hap2, 
                                 max_distance=args.max_distance,
                                 small_distance=args.small_distance, 
                                 similarity_threshold=args.similarity_threshold) if indel_hap1 or indel_hap2 else []
    
    merged_sv = merge_variants(sv_hap1, sv_hap2, 
                              max_distance=args.max_distance,
                              small_distance=args.small_distance, 
                              similarity_threshold=args.similarity_threshold) if sv_hap1 or sv_hap2 else []
    
    all_merged = merged_snv + merged_indel + merged_sv
    all_merged.sort(key=lambda v: v.pos)
    
    logging.info(f"Finished processing chromosome: {chrom}")
    return all_merged

def main(args):
    """Multi thread"""
    chromosomes = get_all_chromosomes(args)
    logging.info(f"Found chromosomes: {', '.join(chromosomes)}")
    
    with open(args.output, 'w') as output_handle:
        write_vcf_header(output_handle, args, args.sample_name)
    
    ctx = multiprocessing.get_context('spawn')
    pool = ctx.Pool(processes=min(len(chromosomes), multiprocessing.cpu_count()))
    
    results = []
    try:
        process_func = partial(process_chromosome, args=args)
        
        # show process
        with tqdm(total=len(chromosomes), desc="Merge variants", unit="variant") as pbar:
            for i, chrom_variants in enumerate(pool.imap_unordered(process_func, chromosomes)):
                results.append((chromosomes[i], chrom_variants))
                pbar.update(1)
    except Exception as e:
        logging.error(f"Error during processing: {e}")
        pool.terminate()
        raise
    finally:
        pool.close()
        pool.join()
    
    chrom_index = {chrom: idx for idx, chrom in enumerate(chromosomes)}
    results.sort(key=lambda x: chrom_index[x[0]])

    with open(args.output, 'a') as output_handle:
        for chrom, variants in results:
            for variant in variants:
                write_variant(variant, output_handle)
            logging.info(f"Written chromosome: {chrom}")
    
    logging.info(f"Final VCF written to: {args.output}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Merge variants of two haplotypes with chromosome-wise multiprocessing')
    parser.add_argument('--hap1_snv', help='SNV VCF of hap1')
    parser.add_argument('--hap2_snv', help='SNV VCF of hap2')
    parser.add_argument('--hap1_indel', help='INDEL VCF of hap1')
    parser.add_argument('--hap2_indel', help='INDEL VCF of hap2')
    parser.add_argument('--hap1_sv', help='SV VCF of hap1')
    parser.add_argument('--hap2_sv', help='SV VCF of hap2')
    parser.add_argument('--output', required=True, help='Output VCF')
    parser.add_argument('--sample_name', default="SAMPLE", help='sample name')

    parser.add_argument('--max_distance', type=int, default=500, 
                       help='Max reference distance for two allele to merge [500bp].')
    parser.add_argument('--small_distance', type=int, default=10, 
                       help='Max reference distance for SSV (small variants) merge [10bp].')
    parser.add_argument('--similarity_threshold', type=float, default=0.8, 
                       help='The similarity of variants, used to merge the variants of two haplotypes [0.8].')
    
    args = parser.parse_args()

    if not any([args.hap1_snv, args.hap1_indel, args.hap1_sv]):
        parser.error("Input at least one vcf")
    if not any([args.hap2_snv, args.hap2_indel, args.hap2_sv]):
        parser.error("Input at least one vcf")
    
    main(args)
