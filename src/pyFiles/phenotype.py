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
_SPECIAL_ORDER = {'X': 23, 'Y': 24, 'MT': 25, 'M': 25}

def chrom_sort_key(c):
    name = c[3:] if c.startswith('chr') else c        
    if name.isdigit():
        return (int(name), '')
    if name.upper() in _SPECIAL_ORDER:
        return (_SPECIAL_ORDER[name.upper()], '')
    return (99, name)                                   


class Variant:
    __slots__ = (
        'chrom', 'pos', 'ref', 'alt', 'var_type', 'end', 'svlen', 'svtype',
        'info', 'haplotype', 'id', 'qual', 'filter', 'format', 'sample', 'genotype'
    )

    def __init__(self, chrom, pos, ref, alt, var_type,
                 end=None, svlen=None, svtype=None, info=None):
        self.chrom = chrom
        self.pos = pos
        self.ref = ref
        self.alt = alt
        self.var_type = var_type          # 'SNV', 'INDEL', 'SV'
        self.end = end if end is not None else pos + len(ref) - 1
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
    variants = []
    opener = gzip.open if vcf_file.endswith('.gz') else open
    mode = 'rt' if vcf_file.endswith('.gz') else 'r'

    with opener(vcf_file, mode) as f:
        for line in f:
            if line.startswith('#'):
                continue

            fields = line.rstrip('\n').split('\t')
            if len(fields) < 5:
                continue

            chrom_field = fields[0]
            if chrom and chrom_field != chrom:
                continue

            pos, var_id, ref, alt = fields[1], fields[2], fields[3], fields[4]
            qual   = fields[5] if len(fields) > 5 else '.'
            filt   = fields[6] if len(fields) > 6 else '.'
            info   = fields[7] if len(fields) > 7 else '.'
            fmt    = fields[8] if len(fields) > 8 else '.'
            sample = fields[9] if len(fields) > 9 else '.'

            info_dict = {}
            for item in info.split(';'):
                if '=' in item:
                    k, v = item.split('=', 1)
                    info_dict[k] = v
                else:
                    info_dict[item] = True

            svtype = info_dict.get('SVTYPE', None)

            if file_type:
                var_type = file_type.upper()
            elif svtype:
                var_type = 'SNV' if svtype == 'SNV' else 'SV'
            else:
                var_type = 'SNV' if (len(ref) == 1 and len(alt) == 1) else 'INDEL'

            ipos = int(pos)
            if var_type == 'SNV':
                end   = ipos + len(ref) - 1
                svlen = abs(len(alt) - len(ref))
            else:
                if 'END' in info_dict:
                    end   = int(info_dict['END'])
                    svlen = end - ipos
                elif 'SVLEN' in info_dict:
                    svlen_val = info_dict['SVLEN']
                    svlen = abs(int(svlen_val))
                    end   = ipos + svlen
                else:
                    end   = ipos
                    svlen = 0

            genotype = None
            if fmt and sample:
                fmt_fields    = fmt.split(':')
                sample_fields = sample.split(':')
                if 'GT' in fmt_fields:
                    genotype = sample_fields[fmt_fields.index('GT')]

            variant = Variant(
                chrom=chrom_field, pos=ipos, ref=ref, alt=alt,
                var_type=var_type, end=end, svlen=svlen, svtype=svtype,
                info=info_dict
            )
            variant.id       = var_id
            variant.qual     = qual
            variant.filter   = filt
            variant.format   = fmt
            variant.sample   = sample
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

        last  = current_cluster[-1]
        dist  = max_distance if variant.var_type == 'SV' else small_distance

        if (last.chrom    != variant.chrom or
                last.var_type != variant.var_type or
                abs(last.pos  - variant.pos) > dist):
            clusters.append(current_cluster)
            current_cluster = [variant]
        else:
            current_cluster.append(variant)

    if current_cluster:
        clusters.append(current_cluster)

    return clusters


def _qual_float(v):
    try:
        return float(v.qual)
    except (TypeError, ValueError):
        return 0.0


# def compute_variant_similarity(v1, v2, max_distance=500, small_distance=10):
#     """Calculate similarity between two variants (1.0=identical, 0.0=completely different)."""
#     if v1.chrom != v2.chrom or v1.var_type != v2.var_type:
#         return 0.0

#     pos_diff = abs(v1.pos - v2.pos)

#     if v1.var_type in ('SNV', 'INDEL'):
#         if pos_diff > small_distance:
#             return 0.0
#     else:
#         if pos_diff > max_distance:
#             return 0.0

#     if v1.var_type == 'SNV':
#         return 1.0 if (v1.ref == v2.ref and v1.alt == v2.alt) else 0.0

#     if v1.var_type == 'INDEL':
#         seq1 = v1.alt if len(v1.ref) == 1 else v1.ref
#         seq2 = v2.alt if len(v2.ref) == 1 else v2.ref
#         if not seq1 or not seq2:
#             return 0.0
#         max_len = max(len(seq1), len(seq2))
#         align = edlib.align(seq1, seq2, task='path')
#         return 1 - align['editDistance'] / max_len

#     # SV
#     svtype = v1.svtype

#     def kmer_similarity(s1, s2, k=21):
#         if len(s1) < k or len(s2) < k:
#             max_len = max(len(s1), len(s2))
#             if max_len == 0:
#                 return 0.0
#             ed = edlib.align(s1, s2, task='distance')['editDistance']
#             return 1 - ed / max_len

#         kmers1 = Counter(s1[i:i + k] for i in range(len(s1) - k + 1))
#         kmers2 = Counter(s2[i:i + k] for i in range(len(s2) - k + 1))
#         intersection = sum((kmers1 & kmers2).values())
#         union        = sum((kmers1 | kmers2).values())
#         return intersection / union if union else 0.0

#     _sym_re = re.compile(r'^<[^>]+>')

#     if svtype == 'INS':
#         seq1 = _sym_re.sub('', v1.alt)
#         seq2 = _sym_re.sub('', v2.alt)
#         return kmer_similarity(seq1, seq2) if (seq1 and seq2) else 0.0

#     if svtype == 'DEL':
#         seq1 = _sym_re.sub('', v1.ref)
#         seq2 = _sym_re.sub('', v2.ref)
#         return kmer_similarity(seq1, seq2) if (seq1 and seq2) else 0.0

#     if svtype == 'INV':
#         overlap = min(v1.end, v2.end) - max(v1.pos, v2.pos)
#         min_len = max(v1.end - v1.pos, v2.end - v2.pos)
#         return 1.0 if (min_len > 0 and overlap >= 0.8 * min_len) else 0.0

#     return 0.0

def compute_variant_similarity(v1, v2, max_distance=500, small_distance=10):
    """优化后的变异相似度计算函数（按长度分流 INS，按重叠率计算 DEL）"""
    if v1.chrom != v2.chrom or v1.var_type != v2.var_type:
        return 0.0

    pos_diff = abs(v1.pos - v2.pos)

    if v1.var_type in ('SNV', 'INDEL'):
        if pos_diff > small_distance:
            return 0.0
    else:
        if pos_diff > max_distance:
            return 0.0

    if v1.var_type == 'SNV':
        return 1.0 if (v1.ref == v2.ref and v1.alt == v2.alt) else 0.0

    if v1.var_type == 'INDEL':
        seq1 = v1.alt if len(v1.ref) == 1 else v1.ref
        seq2 = v2.alt if len(v2.ref) == 1 else v2.ref
        if not seq1 or not seq2:
            return 0.0
        max_len = max(len(seq1), len(seq2))
        align = edlib.align(seq1, seq2, task='path')
        return 1 - align['editDistance'] / max_len

    # ==================== 3. 结构变异 (SV) 核心优化逻辑 ====================
    svtype = v1.svtype
    _sym_re = re.compile(r'^<[^>]+>')

    # ------------------ CASE A: DELETION ------------------
    # if svtype == 'DEL':
    #     # 计算交集区域 (Overlap length)
    #     overlap_start = max(v1.pos, v2.pos)
    #     overlap_end = min(v1.end, v2.end)
    #     overlap_len = max(0, overlap_end - overlap_start + 1)
        
    #     if overlap_len == 0:
    #         return 0.0

    #     len1 = v1.end - v1.pos + 1
    #     len2 = v2.end - v2.pos + 1

    #     overlap_ratio1 = overlap_len / len1
    #     overlap_ratio2 = overlap_len / len2
        
    #     if overlap_ratio1 >= 0.8 and overlap_ratio2 >= 0.8:
    #         return 1.0
    #     else:
    #         return 0.0

    # ------------------ CASE B: INSERTION ------------------
    if svtype == 'INS':
        seq1 = _sym_re.sub('', v1.alt)
        seq2 = _sym_re.sub('', v2.alt)
        if not seq1 or not seq2:
            return 0.0

        max_len = max(len(seq1), len(seq2))

        if max_len < 1000:
            align = edlib.align(seq1, seq2, task='distance')
            return 1 - align['editDistance'] / max_len
        else:
            k = 21
            kmers1 = Counter(seq1[i:i + k] for i in range(len(seq1) - k + 1))
            kmers2 = Counter(seq2[i:i + k] for i in range(len(seq2) - k + 1))
            intersection = sum((kmers1 & kmers2).values())
            union        = sum((kmers1 | kmers2).values())
            return intersection / union if union else 0.0
        
    if svtype == 'DEL':
        seq1 = _sym_re.sub('', v1.ref)
        seq2 = _sym_re.sub('', v2.ref)
        if not seq1 or not seq2:
            return 0.0

        max_len = max(len(seq1), len(seq2))

        if max_len < 1000:
            align = edlib.align(seq1, seq2, task='distance')
            return 1 - align['editDistance'] / max_len
        else:
            k = 21
            kmers1 = Counter(seq1[i:i + k] for i in range(len(seq1) - k + 1))
            kmers2 = Counter(seq2[i:i + k] for i in range(len(seq2) - k + 1))
            intersection = sum((kmers1 & kmers2).values())
            union        = sum((kmers1 | kmers2).values())
            return intersection / union if union else 0.0

    # ------------------ CASE C: INVERSION ------------------
    if svtype == 'INV':
        overlap = min(v1.end, v2.end) - max(v1.pos, v2.pos)
        min_len = max(v1.end - v1.pos, v2.end - v2.pos)
        return 1.0 if (min_len > 0 and overlap >= 0.8 * min_len) else 0.0

    return 0.0

# def cluster_variants(variants, similarity_threshold=0.8, max_distance=500, small_distance=10):
#     """Cluster variants based on similarity threshold."""
#     n = len(variants)
#     if n <= 1:
#         return [variants]

#     if n == 2:
#         sim = compute_variant_similarity(
#             variants[0], variants[1], max_distance, small_distance)
#         if sim >= similarity_threshold:
#             return [variants]
#         return [[variants[0]], [variants[1]]]

#     def sim_to_dist(v1, v2):
#         return 1 - compute_variant_similarity(v1, v2, max_distance, small_distance)

#     condensed = [sim_to_dist(variants[i], variants[j])
#                  for i in range(n) for j in range(i + 1, n)]

#     Z        = linkage(condensed, method='complete')
#     labels   = fcluster(Z, 1 - similarity_threshold, criterion='distance')

#     cluster_dict = defaultdict(list)
#     for idx, cid in enumerate(labels):
#         cluster_dict[cid].append(variants[idx])

#     return list(cluster_dict.values())

def cluster_variants(variants, similarity_threshold=0.8, max_distance=500, small_distance=10):
    """聚类函数：增加了同一 Haplotype 隔离约束"""
    n = len(variants)
    if n <= 1:
        return [variants]

    if n == 2:
        if variants[0].haplotype == variants[1].haplotype:

            return [[variants[0]], [variants[1]]]
            
        sim = compute_variant_similarity(
            variants[0], variants[1], max_distance, small_distance)
        if sim >= similarity_threshold:
            return [variants]
        return [[variants[0]], [variants[1]]]

    def sim_to_dist(v1, v2):
        if v1.haplotype == v2.haplotype:
            return 1.0
            
        return 1 - compute_variant_similarity(v1, v2, max_distance, small_distance)

    condensed = [sim_to_dist(variants[i], variants[j])
                 for i in range(n) for j in range(i + 1, n)]

    # 使用 complete linkage 确保 Cluster 内部任意两者的距离都小于阈值
    Z        = linkage(condensed, method='complete')
    labels   = fcluster(Z, 1 - similarity_threshold, criterion='distance')

    cluster_dict = defaultdict(list)
    for idx, cid in enumerate(labels):
        cluster_dict[cid].append(variants[idx])

    return list(cluster_dict.values())


def merge_info_fields(cluster, base_info,
                      fields_to_merge=('TIG_REGION', 'QUERY_STRAND')):
    for field in fields_to_merge:
        combined = []
        for v in cluster:
            if field in v.info and v.info[field]:
                combined.extend(str(v.info[field]).split(','))
        if combined:
            unique = list(dict.fromkeys(x.strip() for x in combined if x.strip()))
            base_info[field] = ','.join(unique)
    return base_info


def _finalize_variant(variant, cluster, merged_id, base_info, genotype):
    variant.id       = merged_id
    variant.info     = base_info
    variant.genotype = genotype
    variant.pos      = min(v.pos for v in cluster)
    variant.end      = max(v.end for v in cluster)
    if variant.var_type == 'SV' and 'END' in base_info:
        base_info['END'] = str(variant.end)
    return variant


def merge_variants(variants_hap1, variants_hap2,
                   max_distance, small_distance, similarity_threshold):
    for v in variants_hap1:
        v.haplotype = 1
    for v in variants_hap2:
        v.haplotype = 2

    all_variants = variants_hap1 + variants_hap2
    partitions   = form_clusters(all_variants, max_distance, small_distance)
    merged       = []

    for partition in partitions:
        clusters = cluster_variants(
            partition, similarity_threshold, max_distance, small_distance)

        for cluster in clusters:
            if not cluster:
                continue

            merged_id = cluster[0].id
            base_info = copy.deepcopy(cluster[0].info)
            base_info = merge_info_fields(cluster, base_info)

            if len(cluster) == 1:
                v        = cluster[0]
                genotype = '1|0' if v.haplotype == 1 else '0|1'
                merged.append(_finalize_variant(v, cluster, merged_id, base_info, genotype))

            elif len(cluster) == 2:
                v1, v2 = cluster
                best     = v1 if _qual_float(v1) >= _qual_float(v2) else v2
                variant  = copy.deepcopy(best)
                genotype = ('1|0' if v1.haplotype == 1 else '0|1') \
                           if v1.haplotype == v2.haplotype else '1|1'
                merged.append(_finalize_variant(
                    variant, cluster, merged_id, base_info, genotype))

            else:
                hap1_vars = [v for v in cluster if v.haplotype == 1]
                hap2_vars = [v for v in cluster if v.haplotype == 2]

                best_hap1 = max(hap1_vars, key=_qual_float) if hap1_vars else None
                best_hap2 = max(hap2_vars, key=_qual_float) if hap2_vars else None

                if best_hap1 and best_hap2:
                    rep      = best_hap1 if _qual_float(best_hap1) >= _qual_float(best_hap2) \
                               else best_hap2
                    variant  = copy.deepcopy(rep)
                    genotype = '1|1'
                else:
                    best_var = best_hap1 or best_hap2
                    variant  = copy.deepcopy(best_var)
                    genotype = '1|0' if best_hap1 else '0|1'

                merged.append(_finalize_variant(
                    variant, cluster, merged_id, base_info, genotype))

    return merged

def get_header_file(args):
    for f in [args.hap1_sv, args.hap2_sv,
              args.hap1_indel, args.hap2_indel,
              args.hap1_snv,  args.hap2_snv]:
        if f and os.path.exists(f):
            return f
    return None


def parse_source_header(vcf_file):
    headers = []
    opener = gzip.open if vcf_file.endswith('.gz') else open
    with opener(vcf_file, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                headers.append(line.rstrip('\n'))
            else:
                break
    return headers


def write_vcf_header(output_handle, args, sample_name='SAMPLE'):
    output_handle.write('##fileformat=VCFv4.2\n')
    output_handle.write(f'##fileDate={time.strftime("%Y%m%d")}\n')
    output_handle.write('##source=LGvar\n')

    header_file = get_header_file(args)
    if header_file:
        for line in parse_source_header(header_file):
            if line.startswith(('##reference', '##contig')):
                output_handle.write(line + '\n')

    output_handle.write('##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">\n')
    output_handle.write('##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length of the variant">\n')
    output_handle.write('##INFO=<ID=TIG_REGION,Number=.,Type=String,Description="Contig region where variant was found">\n')
    output_handle.write('##INFO=<ID=QUERY_STRAND,Number=.,Type=String,Description="Strand of variant in the contig">\n')
    output_handle.write('##INFO=<ID=ID,Number=1,Type=String,Description="Original variant ID">\n')
    output_handle.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
    output_handle.write(f'#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample_name}\n')


def write_variant(variant, output_handle):
    info_parts = []
    for k, v in variant.info.items():
        info_parts.append(k if v is True else f'{k}={v}')
    if variant.id and variant.id != '.':
        info_parts.append(f'ID={variant.id}')

    output_handle.write('\t'.join([
        variant.chrom,
        str(variant.pos),
        variant.id or '.',
        variant.ref,
        variant.alt,
        variant.qual or '.',
        'PASS',
        ';'.join(info_parts),
        'GT',
        variant.genotype or '.|.'
    ]) + '\n')

def get_all_chromosomes(args):
    chromosomes = set()
    input_files = [
        args.hap1_sv,    args.hap2_sv,
        args.hap1_indel, args.hap2_indel,
        args.hap1_snv,   args.hap2_snv,
    ]

    for vcf_file in input_files:
        if not vcf_file or not os.path.exists(vcf_file):
            continue

        opener = gzip.open if vcf_file.endswith('.gz') else open
        found_contig_in_header = False

        with opener(vcf_file, 'rt') as f:
            for line in f:
                if line.startswith('##contig'):
                    m = re.search(r'<ID=([^,>]+)', line)
                    if m:
                        chromosomes.add(m.group(1))
                        found_contig_in_header = True
                elif line.startswith('#'):
                    continue
                else:
                    if not found_contig_in_header:
                        chromosomes.add(line.split('\t', 1)[0])
                    else:
                        chromosomes.add(line.split('\t', 1)[0])

    return sorted(chromosomes, key=chrom_sort_key)

def process_chromosome(chrom, args):
    """Variants of each chromosome"""
    logging.info(f'Processing chromosome: {chrom}')

    snv_hap1   = parse_vcf(args.hap1_snv,   chrom, 'snv')   if args.hap1_snv   else []
    snv_hap2   = parse_vcf(args.hap2_snv,   chrom, 'snv')   if args.hap2_snv   else []
    indel_hap1 = parse_vcf(args.hap1_indel, chrom, 'indel') if args.hap1_indel else []
    indel_hap2 = parse_vcf(args.hap2_indel, chrom, 'indel') if args.hap2_indel else []
    sv_hap1    = parse_vcf(args.hap1_sv,    chrom, 'sv')    if args.hap1_sv    else []
    sv_hap2    = parse_vcf(args.hap2_sv,    chrom, 'sv')    if args.hap2_sv    else []

    kw = dict(max_distance=args.max_distance,
              small_distance=args.small_distance,
              similarity_threshold=args.similarity_threshold)

    merged_snv   = merge_variants(snv_hap1,   snv_hap2,   **kw) if (snv_hap1   or snv_hap2)   else []
    merged_indel = merge_variants(indel_hap1, indel_hap2, **kw) if (indel_hap1 or indel_hap2) else []
    merged_sv    = merge_variants(sv_hap1,    sv_hap2,    **kw) if (sv_hap1    or sv_hap2)    else []

    all_merged = merged_snv + merged_indel + merged_sv
    all_merged.sort(key=lambda v: v.pos)

    logging.info(f'Finished chromosome: {chrom} ({len(all_merged)} variants)')
    return chrom, all_merged


def main(args):
    chromosomes = get_all_chromosomes(args)
    logging.info(f'Found {len(chromosomes)} chromosomes: {", ".join(chromosomes)}')

    with open(args.output, 'w') as fh:
        write_vcf_header(fh, args, args.sample_name)

    ctx  = multiprocessing.get_context('spawn')
    pool = ctx.Pool(processes=min(10, multiprocessing.cpu_count()))

    results = {}
    process_func = partial(process_chromosome, args=args)

    try:
        with tqdm(total=len(chromosomes), desc='Merge variants', unit='chrom') as pbar:
            for chrom, variants in pool.imap_unordered(process_func, chromosomes):
                results[chrom] = variants
                pbar.update(1)
    except Exception as e:
        logging.error(f'Error during processing: {e}')
        pool.terminate()
        raise
    finally:
        pool.close()
        pool.join()

    with open(args.output, 'a') as fh:
        for chrom in chromosomes:                    
            for variant in results.get(chrom, []):
                write_variant(variant, fh)
            logging.info(f'Written chromosome: {chrom}')

    logging.info(f'Done. Output: {args.output}')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Merge variants of two haplotypes with chromosome-wise multiprocessing')
    parser.add_argument('--hap1_snv',   help='SNV VCF of hap1')
    parser.add_argument('--hap2_snv',   help='SNV VCF of hap2')
    parser.add_argument('--hap1_indel', help='INDEL VCF of hap1')
    parser.add_argument('--hap2_indel', help='INDEL VCF of hap2')
    parser.add_argument('--hap1_sv',    help='SV VCF of hap1')
    parser.add_argument('--hap2_sv',    help='SV VCF of hap2')
    parser.add_argument('--output',     required=True, help='Output VCF')
    parser.add_argument('--sample_name', default='SAMPLE', help='Sample name in VCF header')
    parser.add_argument('--max_distance',       type=int,   default=500,
                        help='Max distance for SV merge [500 bp]')
    parser.add_argument('--small_distance',     type=int,   default=10,
                        help='Max distance for SNV/INDEL merge [10 bp]')
    parser.add_argument('--similarity_threshold', type=float, default=0.8,
                        help='Sequence similarity threshold for merging [0.8]')

    args = parser.parse_args()

    if not any([args.hap1_snv, args.hap1_indel, args.hap1_sv]):
        parser.error('Provide at least one hap1 VCF (--hap1_snv / --hap1_indel / --hap1_sv)')
    if not any([args.hap2_snv, args.hap2_indel, args.hap2_sv]):
        parser.error('Provide at least one hap2 VCF (--hap2_snv / --hap2_indel / --hap2_sv)')

    main(args)
