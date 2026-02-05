## --*-- coding=iso-8859-1 --*--
"""
Created on Wed Dec 24 16:56 2025
Modified from R script by Shilong Zhang
@Author: Feifei Zhou
@Author_email: zhoufeifei@sjtu.edu.cn
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from matplotlib.backends.backend_pdf import PdfPages
import argparse

def _merge_alignments(
    df,
    reference_name_col="reference_name",
    reference_start_col="reference_start",
    reference_end_col="reference_end",
    strand_col="strand",
    query_name_col="query_name",
    query_start_col="query_start",
    query_end_col="query_end",
    min_alignment_length=50000,
    mxm_distance=1000000
):
    # Filter short alignment
    df = df[
        (df[reference_end_col] - df[reference_start_col] >= min_alignment_length) |
        (df[query_end_col] - df[query_start_col] >= min_alignment_length)
    ]
    
    # Sort
    df = df.sort_values([
        reference_name_col, 
        query_name_col, 
        strand_col, 
        reference_start_col, 
        query_start_col
    ]).reset_index(drop=True)
    
    if df.empty:
        return pd.DataFrame()

    result = []
    current_group = df.iloc[0].copy()

    for i in range(1, len(df)):
        row = df.iloc[i]
        same_ref = current_group[reference_name_col] == row[reference_name_col]
        same_query = current_group[query_name_col] == row[query_name_col]
        same_strand = current_group[strand_col] == row[strand_col]
        
        ref_dist = abs(row[reference_start_col] - current_group[reference_end_col])
        query_dist = abs(row[query_start_col] - current_group[query_end_col])

        if all([same_ref, same_query, same_strand, 
               ref_dist <= mxm_distance, query_dist <= mxm_distance]):
            current_group[reference_end_col] = max(current_group[reference_end_col], row[reference_end_col])
            current_group[query_end_col] = max(current_group[query_end_col], row[query_end_col])
        else:
            result.append(current_group)
            current_group = row.copy()

    result.append(current_group)
    return pd.DataFrame(result)

def pairwise_alignment_plot(
    ax,
    df,
    r_order=None,
    q_order=None,
    reference_name_col="reference_name",
    reference_start_col="reference_start",
    reference_end_col="reference_end",
    strand_col="strand",
    query_name_col="query_name",
    query_start_col="query_start",
    query_end_col="query_end",
    merge=True,
    mxm_distance=1000000,
    min_alignment_length=50000,
    chrom_distance=10000000,
    inv_middle_distance=10000,
    y_padding=10000000
):
    sdf = df.copy()
    
    # Merge alignment
    if merge:
        sdf = _merge_alignments(
            sdf,
            reference_name_col=reference_name_col,
            reference_start_col=reference_start_col,
            reference_end_col=reference_end_col,
            strand_col=strand_col,
            query_name_col=query_name_col,
            query_start_col=query_start_col,
            query_end_col=query_end_col,
            min_alignment_length=min_alignment_length,
            mxm_distance=mxm_distance
        )
    
    r_order = r_order if r_order is not None else sdf[reference_name_col].unique()
    #r_order = [f'chr{i}' for i in range(1, 23)] + ['chrX', 'chrY']
    q_order = q_order if q_order is not None else sdf[query_name_col].unique()

    r_sizes = sdf.groupby(reference_name_col)[reference_end_col].max().reset_index()
    r_sizes.rename(columns={reference_name_col: 'chrom', reference_end_col: 'size'}, inplace=True)
    q_sizes = sdf.groupby(query_name_col)[query_end_col].max().reset_index()
    q_sizes.rename(columns={query_name_col: 'chrom', query_end_col: 'size'}, inplace=True)

    r_shift = {}
    current = 0
    ref_total = 0
    for chrom in r_order:
        r_shift[chrom] = current + len(r_shift) * chrom_distance
        size = r_sizes[r_sizes['chrom'] == chrom]['size'].values[0]
        ref_total += size + chrom_distance
        current += size

    q_shift = {}
    current = 0
    query_total = 0
    for chrom in q_order:
        q_shift[chrom] = current + len(q_shift) * chrom_distance
        size = q_sizes[q_sizes['chrom'] == chrom]['size'].values[0]
        query_total += size + chrom_distance
        current += size

    sdf[reference_start_col] = sdf[reference_name_col].map(r_shift) + sdf[reference_start_col]
    sdf[reference_end_col] = sdf[reference_name_col].map(r_shift) + sdf[reference_end_col]
    sdf[query_start_col] = sdf[query_name_col].map(q_shift) + sdf[query_start_col]
    sdf[query_end_col] = sdf[query_name_col].map(q_shift) + sdf[query_end_col]

    sdf['.reference_mid'] = (sdf[reference_start_col] + sdf[reference_end_col]) / 2
    sdf['.query_mid'] = (sdf[query_start_col] + sdf[query_end_col]) / 2
    sdf['.both_mid'] = (sdf['.reference_mid'] + sdf['.query_mid']) / 2

    polygons = []
    colors = []
    
    for _, row in sdf.iterrows():
        if row[strand_col] == '+':
            verts = [
                (0, row[reference_start_col]),
                (2, row[query_start_col]),
                (2, row[query_end_col]),
                (0, row[reference_end_col])
            ]
            color = '#2081f9'
        elif row[strand_col] == '-':
            verts = [
                (0, row[reference_start_col]),
                (2, row[query_start_col] + row[reference_end_col] - row[reference_start_col]),
                (2, row[query_start_col]),
                (0, row[reference_end_col]),
                (0, row[reference_start_col])
            ]
            color = '#f99820'
        else:
            raise ValueError("Invalid strand value")
            
        polygons.append(verts)
        colors.append(color)

    coll = PolyCollection(
        polygons,
        facecolors=colors,
        edgecolors='none',
        alpha=0.6
    )
    ax.add_collection(coll)

    for chrom in r_order:
        shift = r_shift.get(chrom, 0)
        size = r_sizes[r_sizes['chrom'] == chrom]['size'].values[0]
        ax.plot([0, 0], [shift, shift + size], color='#444444', lw=2)
    
    for chrom in q_order:
        shift = q_shift.get(chrom, 0)
        size = q_sizes[q_sizes['chrom'] == chrom]['size'].values[0]
        ax.plot([2, 2], [shift, shift + size], color='#444444', lw=2)

    max_y = max(ref_total, query_total) + y_padding
    ax.set_ylim(0, max_y)
    ax.set_xlim(-0.5, 2.5)
    ax.axis('off')

    ax.text(-0.25, -15000000, ','.join(r_order), ha='left', fontsize = 12)
    ax.text(2.25, -15000000, ','.join(q_order), ha='right', fontsize = 12)
    #ax.set_title(f"{r_order[0]} vs {q_order[0]}", pad=20)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Visualize pairwise alignments from a PAF file.")
    parser.add_argument("-p", "--paf", required=True, help="Input PAF file")
    parser.add_argument("-o", "--output", default="alignments.pdf", help="Output PDF file name (default: alignments.pdf)")
    parser.add_argument("-m", "--minimum", type=int, default=50000, help="Filter alignments shorter than this length (default: 50000)")
    parser.add_argument("-d", "--distance", type=int, default=1000000, help="Params for merge alignments (default: 1000000)")
    parser.add_argument("-f", "--file", required=True, help="Chrom pair file for homo chromosome alignments plot")
    
    args = parser.parse_args()
    
    with open(args.paf, 'r') as f:
        line = f.readline().strip()
        num_columns = len(line.split('\t'))

    predefined_columns = [
        'query_name', 'query_length', 'query_start', 'query_end', 'strand',
        'reference_name', 'reference_length', 'reference_start', 'reference_end',
        'n_matches', 'alignment_block_length', 'mapped_length', 'score'
    ]
    
    #additional_columns_needed = num_columns - len(predefined_columns)
    additional_columns_needed = 11
    additional_columns = [f'column{i+1}' for i in range(additional_columns_needed)]
    
    all_columns = predefined_columns + additional_columns
    
    df = pd.read_csv(args.paf, sep='\t', names=all_columns)
    #df = df[(df['reference_name'] != "chrX")]
    df = df[(df['reference_name'] != "chrY")]
    df = df[(df['reference_name'] != "chrM")]

    r_order = df['reference_name'].unique()
    #print(r_order)

    with PdfPages(args.output) as pdf:
        plots_per_page = 24
        rows_per_page = 3
        cols_per_page = 8
        
        total_pages = (len(r_order) + plots_per_page - 1) // plots_per_page
        
        for page in range(total_pages):
            fig, axes = plt.subplots(rows_per_page, cols_per_page, figsize=(40, 15))
            fig.subplots_adjust(hspace=0.5, wspace=0.3)
            
            axes = axes.flatten()
 
            start_idx = page * plots_per_page
            end_idx = min((page + 1) * plots_per_page, len(r_order))
            
            for i, ax in enumerate(axes):
                idx = start_idx + i
                if idx >= end_idx:
                    ax.axis('off')
                    continue
                
                chromosome_order = [f'chr{i}' for i in range(1, 23)]
                order_dict = {chrom: idx for idx, chrom in enumerate(chromosome_order)}
                sort_indices = np.argsort([order_dict.get(chrom, len(chromosome_order)) for chrom in r_order])

                r_order_sorted = r_order[sort_indices]
                chrom = r_order_sorted[idx]
                chrom_df = df[df['reference_name'] == chrom]
                
                if chrom_df.empty:
                    ax.axis('off')
                    continue

                pairwise_alignment_plot(
                    ax,
                    chrom_df,
                    reference_name_col='reference_name',
                    reference_start_col='reference_start',
                    reference_end_col='reference_end',
                    query_name_col='query_name',
                    query_start_col='query_start',
                    query_end_col='query_end',
                    strand_col='strand',
                    min_alignment_length=args.minimum,
                    merge=True,
                    mxm_distance=args.distance
                )
            
            pdf.savefig(fig, bbox_inches='tight')
            plt.close(fig)

