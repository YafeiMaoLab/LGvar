import subprocess
from Bio import SeqIO
import re
import os
from multiprocessing import Pool
import pandas as pd
import numpy as np
import argparse

parser = argparse.ArgumentParser(description='tsv2vcf')
parser.add_argument('--r', help='Reference fasta')
parser.add_argument('--q', help='query fasta')
parser.add_argument('--i', help='SDR result file')
parser.add_argument('--o', help='output dir')
args = parser.parse_args()
QUERY =args.q
REF =args.r
paf=args.i



df = pd.read_csv(paf, sep='\t',header=None)
with open(args.o, 'w') as file:
	## INV
	inv_rows = df[df.iloc[:, 6] == 'INV']
	inv_rows.to_csv('temp/inv.csv', sep='\t', header=None, index=False)
	subprocess.run(f"awk '{{print $4,$5,$6,\"\",\"\",$8}}' temp/inv.csv|tr ' ' '\t' > temp/query.bed", shell=True)
	subprocess.run(f"bedtools getfasta -fi {QUERY} -bed temp/query.bed -s -fo temp/que.fa", shell=True)
	subprocess.run(f"awk '{{print $1,$2,$3}}' temp/inv.csv|tr ' ' '\t'  > temp/ref.bed", shell=True)  
	subprocess.run(f"bedtools getfasta -fi {REF} -bed temp/ref.bed -fo temp/ref.fa", shell=True)
	refsequences = {record.id: record.seq for record in SeqIO.parse('temp/ref.fa', 'fasta')}
	querysequences = {record.id: record.seq for record in SeqIO.parse('temp/que.fa', 'fasta')}
	for index, row in inv_rows.iterrows():
		ref_chr = row[0]
		ref_start = row[1]
		ref_end = row[2]
		que_chr = row[3]
		que_start = row[4]
		que_end = row[5]
		lenall = row[8]
		minus="-"
		selected_reffa = str(refsequences[ref_chr+":"+str(ref_start)+"-"+str(ref_end)])
		selected_quefa = str(querysequences[que_chr+":"+str(que_start)+"-"+str(que_end)+"("+minus+")"])
		three=ref_chr+"-"+str(ref_start)+"-"+"INV"+"-"+str(lenall)
		emp="."
		alllen="ID="+three+";"+"SVTYPE=INV;SVLEN="+str(lenall)+";TIG_REGION="+que_chr+":"+str(que_start)+"-"+str(que_end)+";"+"QUERY_STRAND=+,+"
		GT="GT"
		phase="1|0"
		file.write(f"{ref_chr}\t{ref_start}\t{three}\t{selected_reffa}\t{selected_quefa}\t{emp}\t{emp}\t{alllen}\t{GT}\t{phase}\n")  
	## INS
	ins_rows = df[df.iloc[:, 6] == 'INS']
	ins_rows.to_csv('temp/ins.csv', sep='\t', header=None, index=False)
	subprocess.run(f"awk '{{print $4,$5-1,$6,\"\",\"\",$8}}' temp/ins.csv|tr ' ' '\t' > temp/query.bed", shell=True)
	subprocess.run(f"bedtools getfasta -fi {QUERY} -bed temp/query.bed -s -fo temp/que.fa", shell=True)
	subprocess.run(f"awk '{{print $1,$2-1,$3}}' temp/ins.csv|tr ' ' '\t'  > temp/ref.bed", shell=True)  
	subprocess.run(f"bedtools getfasta -fi {REF} -bed temp/ref.bed -fo temp/ref.fa", shell=True)
	refsequences = {record.id: record.seq for record in SeqIO.parse('temp/ref.fa', 'fasta')}
	querysequences = {record.id: record.seq for record in SeqIO.parse('temp/que.fa', 'fasta')}
	for index, row in ins_rows.iterrows():
		ref_chr = row[0]
		ref_start = row[1]
		ref_end = row[2]
		que_chr = row[3]
		que_start = row[4]
		que_end = row[5]
		lenall = row[9]
		minus=row[7]
		selected_reffa = str(refsequences[ref_chr+":"+str(ref_start-1)+"-"+str(ref_end)])
		selected_quefa = str(querysequences[que_chr+":"+str(que_start-1)+"-"+str(que_end)+"("+minus+")"])
		three=ref_chr+"-"+str(ref_start)+"-"+"INS"+"-"+str(lenall)
		emp="."
		alllen="ID="+three+";"+"SVTYPE=INS;SVLEN="+str(lenall)+";TIG_REGION="+que_chr+":"+str(que_start)+"-"+str(que_end)+";"+"QUERY_STRAND=+,+"
		GT="GT"
		phase="1|0"
		file.write(f"{ref_chr}\t{ref_start}\t{three}\t{selected_reffa}\t{selected_quefa}\t{emp}\t{emp}\t{alllen}\t{GT}\t{phase}\n")  
	## DEL
	del_rows = df[df.iloc[:, 6] == 'DEL']
	del_rows.to_csv('temp/del.csv', sep='\t', header=None, index=False)
	subprocess.run(f"awk '{{print $4,$5-1,$6,\"\",\"\",$8}}' temp/del.csv|tr ' ' '\t' > temp/query.bed", shell=True)
	subprocess.run(f"bedtools getfasta -fi {QUERY} -bed temp/query.bed -s -fo temp/que.fa", shell=True)
	subprocess.run(f"awk '{{print $1,$2-1,$3}}' temp/del.csv|tr ' ' '\t'  > temp/ref.bed", shell=True)  
	subprocess.run(f"bedtools getfasta -fi {REF} -bed temp/ref.bed -fo temp/ref.fa", shell=True)
	refsequences = {record.id: record.seq for record in SeqIO.parse('temp/ref.fa', 'fasta')}
	querysequences = {record.id: record.seq for record in SeqIO.parse('temp/que.fa', 'fasta')}
	for index, row in del_rows.iterrows():
		ref_chr = row[0]
		ref_start = row[1]
		ref_end = row[2]
		que_chr = row[3]
		que_start = row[4]
		que_end = row[5]
		lenall = row[8]
		minus=row[7]
		selected_reffa = str(refsequences[ref_chr+":"+str(ref_start-1)+"-"+str(ref_end)])
		selected_quefa = str(querysequences[que_chr+":"+str(que_start-1)+"-"+str(que_end)+"("+minus+")"])
		three=ref_chr+"-"+str(ref_start)+"-"+"DEL"+"-"+str(lenall)
		emp="."
		alllen="ID="+three+";"+"SVTYPE=DEL;SVLEN="+str(lenall)+";TIG_REGION="+que_chr+":"+str(que_start)+"-"+str(que_end)+";"+"QUERY_STRAND=+,+"
		GT="GT"
		phase="1|0"
		file.write(f"{ref_chr}\t{ref_start}\t{three}\t{selected_reffa}\t{selected_quefa}\t{emp}\t{emp}\t{alllen}\t{GT}\t{phase}\n")
    
rm_file = ["temp/ref.fa", "temp/que.fa", "temp/ins.csv", "temp/del.csv", "temp/inv.csv"]  
for file in rm_file:
    os.remove(file)
