##module load bedtools
##load packages
import re
import time
import logging
import argparse
import subprocess
from Bio import SeqIO
from multiprocessing import Pool

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
parser = argparse.ArgumentParser(description='CIGAR-CALCULATE')
parser.add_argument('--p', help='number of process')
parser.add_argument('--r', help='Reference fasta')
parser.add_argument('--q', help='query fasta')
parser.add_argument('--paf', help='alignment paf file')
parser.add_argument('--o', help='output dir')
args = parser.parse_args()

QUERY =args.q
REF =args.r
paf=args.paf

##get reference and query genome sequences
subprocess.run(f"awk '{{print $1,$3,$4,\"\",\"\",$5}}' {paf} | tr ' ' '\t'  > temp/query.bed", shell=True)
subprocess.run(f"bedtools getfasta -fi {QUERY} -bed temp/query.bed -s -fo temp/que.fa", shell=True)
subprocess.run(f"awk '{{print $6,$8,$9}}' {paf} | tr ' ' '\t'  > temp/ref.bed", shell=True) ##notice the oriention of reference genome
subprocess.run(f"bedtools getfasta -fi {REF} -bed temp/ref.bed -fo temp/ref.fa", shell=True)

refsequences = {record.id: record.seq for record in SeqIO.parse('temp/ref.fa', 'fasta')}
querysequences = {record.id: record.seq for record in SeqIO.parse('temp/que.fa', 'fasta')}

##call SV from CIGAR
def process_cigar(line):
    columns = line.split('\t')
    minus = columns[4]  ##Determine strand
    ref_chr = columns[5]
    ref_start = columns[7]
    ref_end = columns[8]
    query_chr = columns[0]
    query_start = columns[2]
    query_end = columns[3]
    alllen=int(query_end)-int(query_start)
    cigar_str = next((xx for xx in columns if xx.startswith('cg:Z:')), None)  ##get CIGAR string
    cigar_str = cigar_str[5:]
    selected_reffa = refsequences[ref_chr+":"+ref_start+"-"+ref_end]
    selected_quefa = querysequences[query_chr+":"+query_start+"-"+query_end+"("+minus+")"]
    id=ref_chr+":"+ref_start+"-"+ref_end+'_'+query_chr+":"+query_start+"-"+query_end
    ref_pos = 0
    query_pos = 0
    results = []

    ##find all cigar pairs(numbers+types) and record the variations along the string
    i = 0
    while i < len(cigar_str):
        ##search numbers
        num_match = re.match(r'(\d+)', cigar_str[i:])
        if num_match:
            length = int(num_match.group(1))
            i += len(num_match.group(1)) 

            #search cigar digit（=, X, D, I）
            if i < len(cigar_str) and cigar_str[i] in '=XDI':
                variant_type = cigar_str[i]
                i += 1  ##move to the next pair
                
                ##1.match
                if variant_type == '=':
                    ref_pos += length
                    query_pos += length
                    
                ##2.mismatch(SNP)
                elif variant_type == 'X':
                    if length == 1:
                        ref_strseq=selected_reffa[ref_pos:ref_pos+1] 
                        que_strseq=selected_quefa[query_pos:query_pos+1]
                        if(minus=="-"):
                            results.append(f"{id}\t{ref_pos+1}\t{alllen-query_pos}\t1\tSNP\t{ref_strseq}\t{que_strseq}\t{minus}\n") 
                        else:
                            results.append(f"{id}\t{ref_pos+1}\t{query_pos+1}\t1\tSNP\t{ref_strseq}\t{que_strseq}\t{minus}\n")
                        ref_pos += 1
                        query_pos += 1
                    else:
                        for _ in range(length):  ##split into multi SNPs
                            ref_strseq=selected_reffa[ref_pos:ref_pos+1] 
                            que_strseq=selected_quefa[query_pos:query_pos+1]
                            if(minus=="-"):
                                results.append(f"{id}\t{ref_pos+1}\t{alllen-query_pos}\t1\tSNP\t{ref_strseq}\t{que_strseq}\t{minus}\n")
                            else:
                                results.append(f"{id}\t{ref_pos+1}\t{query_pos+1}\t1\tSNP\t{ref_strseq}\t{que_strseq}\t{minus}\n")		
                            ref_pos += 1
                            query_pos += 1
                            
                ##3.deletion
                elif variant_type == 'D':
                    if length == 1:
                        ref_strseq=selected_reffa[ref_pos-1:ref_pos+1] 
                        que_strseq=selected_quefa[query_pos-1:query_pos]
                        if(minus=="-"):
                            results.append(f"{id}\t{ref_pos+1}\t{alllen-query_pos}\t1\tSNP_DEL\t{ref_strseq}\t{que_strseq}\t{minus}\n") 
                        else:
                            results.append(f"{id}\t{ref_pos+1}\t{query_pos+1}\t1\tSNP_DEL\t{ref_strseq}\t{que_strseq}\t{minus}\n")
                    else:
                        ref_strseq=selected_reffa[ref_pos-1:ref_pos+length] 
                        que_strseq=selected_quefa[query_pos-1:query_pos]
                        if(minus=="-"):
                            results.append(f"{id}\t{ref_pos+1}\t{alllen-query_pos}\t{length}\tINDEL_DEL\t{ref_strseq}\t{que_strseq}\t{minus}\n")
                        else:
                            results.append(f"{id}\t{ref_pos+1}\t{query_pos+1}\t{length}\tINDEL_DEL\t{ref_strseq}\t{que_strseq}\t{minus}\n")
                    ref_pos += length
                
                ##4.insertion
                elif variant_type == 'I':
                    if length == 1:
                        ref_strseq=selected_reffa[ref_pos-1:ref_pos] 
                        que_strseq=selected_quefa[query_pos-1:query_pos+1]
                        if(minus=="-"):
                            results.append(f"{id}\t{ref_pos+1}\t{alllen-query_pos}\t1\tSNP_INS\t{ref_strseq}\t{que_strseq}\t{minus}\n")## ref_posf ref_posf query_posf query_posf+1 
                        else:
                            results.append(f"{id}\t{ref_pos+1}\t{query_pos+1}\t1\tSNP_INS\t{ref_strseq}\t{que_strseq}\t{minus}\n")
                    else:
                        ref_strseq=selected_reffa[ref_pos-1:ref_pos] 
                        que_strseq=selected_quefa[query_pos-1:query_pos+length]
                        if(minus=="-"):
                            results.append(f"{id}\t{ref_pos+1}\t{alllen-query_pos-length+1}\t{length}\tINDEL_INS\t{ref_strseq}\t{que_strseq}\t{minus}\n")
                        else:
                            results.append(f"{id}\t{ref_pos+1}\t{query_pos+1}\t{length}\tINDEL_INS\t{ref_strseq}\t{que_strseq}\t{minus}\n")
                    query_pos += length
            else:
                ##other types:add length to the ref and query positions
                ref_pos += length
                query_pos += length
        else:
            break  
            
    with open("temp/"+ref_chr+":"+str(ref_start)+"-"+str(ref_end)+"_"+query_chr+":"+str(query_start)+"-"+str(query_end)+".cigar", 'w') as outfile:
        for line in results:
            outfile.write(line)
    return results



a=time.time()
with open(paf, 'r') as infile, open(args.o, 'w') as outfile:
	write=outfile.write(f"id\tref_start\tquery_start\tlen\ttype\tref_seq\tquery_seq\tQUERYSTRAND\n")
	##if there has many cigar records, the param p is necessary to accelerate this process, default 20.
	with Pool(processes=20) as pool:
		results = pool.map(process_cigar, infile)
	b=time.time()
logging.info(f"Cigar processing time: {(b-a)/60:.2f} minutes")
 
subprocess.run(f"rm temp/query.bed", shell=True)
subprocess.run(f"rm temp/que.fa", shell=True)
subprocess.run(f"rm temp/ref.bed", shell=True)  
subprocess.run(f"rm temp/ref.fa", shell=True)
