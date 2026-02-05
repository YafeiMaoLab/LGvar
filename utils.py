import os
import glob
import time
import sys
import logging
import subprocess
from pathlib import Path

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def run_minimap(ref_path, hap_path, hap_label, chrom_pair_file):
    """1. Whole genome alignments"""
    cmd = [
        "minimap2",
        "-t", "12",
        "-cx", "asm20",
        "--secondary=no",
        "--eqx",
        "-K", "8G",
        "-s", "1000",
        ref_path,
        hap_path,
        "-o", f"align_{hap_label}.paf"
    ]
    
    logging.info("0.Map")
    logging.info(f"Running minimap2 for {hap_label}...")
    
    result = subprocess.run(cmd, check=True)
    
    if result.returncode == 0:
        logging.info(f"Minimap2 completed for {hap_label}")
    
    logging.info(f"Selecting alignments for homologous chromosomes in {hap_label}...")
    process_paf(f"align_{hap_label}.paf",chrom_pair_file,hap_label)

def process_paf(paf_file, chrom_pair_file, hap_label):
    """2. Filter extra alignments"""
    filter_script = Path(os.path.dirname(os.path.abspath(__file__))) / "src/pyFiles/one2multi_filter.py"
    middle_dir = Path("half")
    middle_dir.mkdir(exist_ok=True)
    p_c_chrlen_file = middle_dir/f"{hap_label}_chrlen.txt"
    output_file = f"align_{hap_label}.flt.paf"
    
    cmd = [
        "python", str(filter_script),
        "-m", chrom_pair_file,  ## create automately
        "-f", paf_file,
        "-1", "6",
        "-2", "1"
    ]
    
    logging.info(f"Filtering PAF for {hap_label}...")
    
    with open(output_file, 'w') as f:
        result = subprocess.run(cmd, stdout=f, check=True)
    
    awk_cmd = [
        "awk", "{print $6,$7,$1,$2}",
        f"{output_file}"
    ]
    
    with open(p_c_chrlen_file, 'w') as f:
        subprocess.run(awk_cmd, stdout=f, check=True)

def process_filtering(hap_label, mode, cluster, dellength, centromere, telomere):
    """3. Filter centromere and telomere alignments"""
    nowdic = Path.cwd()
    middle_dir = Path("half")
    middle_dir.mkdir(exist_ok=True)
    
    saffire_dir = Path(f"saffire{hap_label}")
    saffire_dir.mkdir(exist_ok=True)
    p_c_chrlen_file = middle_dir/f"{hap_label}_chrlen.txt"
    output_file = middle_dir/f"{hap_label}_syntenic.tsv"
    
    # run scripts
    chaos_script = Path(os.path.dirname(os.path.abspath(__file__))) / "src/scripts/chaos_filt.r"
    
    if mode == "ctn":
        cmd = [
            "Rscript", str(chaos_script),
            str(p_c_chrlen_file), 
            f"{str(saffire_dir)}/", 
            f"align_{hap_label}.flt.paf", 
            f"align_{hap_label}.final.paf",
            str(cluster),
            str(dellength),
            str(mode)
        ]
        
        logging.info(f"Running chaos filter for {hap_label} in ctn mode...")
        result = subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        
    elif mode == "cts":
        if centromere is None or telomere is None:
            logging.error("Error: When using 'cts' mode, you must provide both --centromere and --telomere arguments.")
            sys.exit(1)
        cmd = [
            "Rscript", str(chaos_script),
            str(p_c_chrlen_file), 
            f"{str(saffire_dir)}/", 
            f"align_{hap_label}.flt.paf", 
            f"align_{hap_label}.final.paf",
            str(cluster),
            str(dellength),
            str(mode),
            str(centromere), str(telomere)
        ]
        
        logging.info(f"Running chaos filter for {hap_label} in cts mode...")
        result = subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    cmd = f"awk 'BEGIN{{OFS=\"\\t\"}} {{print $6, $8, $9, ($8+$9)/2, $1, $3, $4, ($3+$4)/2, $5}}' align_{hap_label}.final.paf | sort -k1,1 -k2,2n > {output_file}"
    subprocess.run(cmd, shell=True, check=True)

def run_cluster_and_call(ref_path, hap_path, hap_label, invcluster, threads, chunk_size, distance, fraction):
    """4. Cluster and SV calling"""
    
    denSDR_dir = Path(f"denSDR{hap_label}")
    denSDR_dir.mkdir(exist_ok=True)
    #dotplot_dir = Path(f"dotplot{hap_label}")
    #dotplot_dir.mkdir(exist_ok=True)
    sdrall_file = denSDR_dir/"SDRall.txt"
    middle_dir = Path("half")
    middle_dir.mkdir(exist_ok=True)
    
    cluster_script = Path(os.path.dirname(os.path.abspath(__file__))) / "src/scripts/denSDR.r"
    fun_script = Path(os.path.dirname(os.path.abspath(__file__))) / "src/scripts/denSDRfun.r"
    
    cmd = [
        "Rscript", str(cluster_script),
        str(fun_script),
        middle_dir/f"{hap_label}_syntenic.tsv",
        middle_dir/f"{hap_label}_chrlen.txt",
        f"{str(denSDR_dir)}/",
        str(invcluster),
        str(distance)
    ]
    
    logging.info(f"SV calling for {hap_label}...")

    result = subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    logging.info(f"SV calling on {hap_label} finished...")
    
    cmd = [
        "cat", f"{denSDR_dir}/*end.tsv", ">", str(sdrall_file)
    ]
    
    subprocess.run(" ".join(cmd), shell=True, check=True)
    
    # realign to find simple inversions from SDR
    realign_script = Path(os.path.dirname(os.path.abspath(__file__))) / "src/vcf/realign.sh"

    cmd = [
        "bash", str(realign_script),
        "-i", str(sdrall_file),
        "-r", ref_path,
        "-q", hap_path,
        "-o", str(denSDR_dir/"SDRall_final.txt"),
        "-t", str(threads),
        "-c", str(chunk_size),
        "-p", f"align_{hap_label}.final.paf",
        "-f", str(fraction)
        
    ]
    
    subprocess.run(cmd, check=True)
    logging.info(f"Inversion re-identification of {hap_label} finished")
    #print(f"{hap_label} realign finished at {time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())}", flush=True)

def run_cigar_processing(ref_path, hap_path, hap_label):
    """ 5.extract variationas from CIGAR """
    cigar_script = Path(os.path.dirname(os.path.abspath(__file__))) / "src/pyFiles/CIGAR.py"
    temp_dir = Path("temp")
    temp_dir.mkdir(exist_ok=True)
    middle_dir = Path("half")
    middle_dir.mkdir(exist_ok=True)
    
    cmd = [
        "python", str(cigar_script),
        "--r", ref_path,
        "--q", hap_path,
        "--paf", f"align_{hap_label}.final.paf",
        "--o", middle_dir/f"{hap_label}cigar.txt"
    ]
    
    logging.info(f"Small variants identification for {hap_label}...")
    
    subprocess.run(cmd, check=True)
    
    output_file = middle_dir/f"{hap_label}cigarend.txt"
    
    with open(output_file, 'w') as dest:
        with open(os.path.join("half", f"{hap_label}cigar.txt"), 'r') as src:
            for line in src:
                dest.write(line)
        
        for file_path in glob.glob(os.path.join("temp", "*.cigar")):
            with open(file_path, 'r') as src:
                for line in src:
                    dest.write(line)

    cmd = [
        "find", str(temp_dir), "-type", "f", "-name", "*.cigar", "-delete"
    ]
    
    subprocess.run(cmd, check=True)
    logging.info(f"Small variants identification for {hap_label} finished...")
    #print(f"CIGAR generated for {hap_label} finished at {time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())}", flush=True)

def run_dup_filtering(hap_label):
    """6. dup filter"""
    dup_filt_script = Path(os.path.dirname(os.path.abspath(__file__))) / "src/vcf/dup_filt.sh"
    filt_script = Path(os.path.dirname(os.path.abspath(__file__))) / "src/pyFiles/dup_filt.py"
    middle_dir = Path("half")
    middle_dir.mkdir(exist_ok=True)
    paf_file = f"align_{hap_label}.final.paf"
    cigar_end_file = middle_dir/f"{hap_label}cigarend.txt"
    output_file = middle_dir/f"{hap_label}cigarout.txt"

    cmd = [
        "bash", str(dup_filt_script),
        paf_file,
        cigar_end_file,
        output_file,
        filt_script
    ]
    
    logging.info("Filtering duplication...")
    subprocess.run(cmd, check=True)

def generate_vcf(ref_path, hap_path, hap_label, variant_type):
    """7. Generate variation results -- vcf"""
    results_dir = Path(f"results")
    results_dir.mkdir(exist_ok=True)

    middle_dir = Path("half")
    middle_dir.mkdir(exist_ok=True)
    cigar_to_vcf_script = Path(os.path.dirname(os.path.abspath(__file__))) / "src/vcf/cigar2vcf.sh"
    py_script = Path(os.path.dirname(os.path.abspath(__file__))) / "src/pyFiles/SDR_vcf.py"
    cigar_out_file = middle_dir/f"{hap_label}cigarout.txt"
    vcf_file = f"{results_dir}/LGvar{hap_label}.vcf"
    sdrall_final_file = f"denSDR{hap_label}/SDRall_final.txt"
    cigarsdr_txt_file = middle_dir/f"{hap_label}cigarsdr.txt"
    lgvarend_bed_file = f"{results_dir}/LGvar{hap_label}.bed"
 
    cmd = [
        "bash", str(cigar_to_vcf_script),
        cigar_out_file,
        vcf_file,
        sdrall_final_file,
        ref_path,
        hap_path,
        cigarsdr_txt_file,
        str(py_script),
        lgvarend_bed_file,
        str(variant_type)
    ]
    
    
    logging.info(f"Generating VCF of {hap_label}...")
    subprocess.run(cmd, check=True)


def split_vcf(hap_label, variant_type):
    """8. split vcf"""
    split_script = Path(os.path.dirname(os.path.abspath(__file__))) / "src/vcf/splitfile.sh"
    hap_dir = Path(f"{hap_label}")
    hap_dir.mkdir(exist_ok=True)
    
    cmd = [
        "bash", str(split_script),
        f"{str(hap_dir)}/",
        f"results/LGvar{hap_label}.vcf",
        str(variant_type)
    ]

    subprocess.run(cmd, check=True)

def integrate_results(hap1_dir, hap2_dir, sample_name, max_distance, small_distance, similarity_threshold, variant_type):
    """9. merge two haplotypes' variation"""
    phenotype_script = Path(os.path.dirname(os.path.abspath(__file__))) / "src/pyFiles/phenotype.py"
    
    variants = list(set(variant_type.split(',')))
    output_vcf = f"sortLGvar_{'_'.join(sorted(variants))}.vcf"

    cmd = [
        "python", str(phenotype_script),
        "--output", f"results/{output_vcf}",
        "--max_distance", str(max_distance),
        "--small_distance", str(small_distance),
        "--similarity_threshold", str(similarity_threshold),
        "--sample_name", sample_name
    ]
    
    cmd_extended = False
    if "snv" in variants or "all" in variants:
        cmd.extend(["--hap1_snv", f"{hap1_dir}/sortsnv.vcf.gz"])
        cmd.extend(["--hap2_snv", f"{hap2_dir}/sortsnv.vcf.gz"])
        cmd_extended = True
    if "ins" in variants or "del" in variants or "all" in variants:
        cmd.extend(["--hap1_indel", f"{hap1_dir}/sortindel.vcf.gz"])
        cmd.extend(["--hap2_indel", f"{hap2_dir}/sortindel.vcf.gz"])
        cmd_extended = True
    if "inv" in variants or "all" in variants:
        cmd.extend(["--hap1_sv", f"{hap1_dir}/sortSV.vcf.gz"])
        cmd.extend(["--hap2_sv", f"{hap2_dir}/sortSV.vcf.gz"])
        cmd_extended = True

    if cmd_extended:
        subprocess.run(cmd, check=True)

    vcf2bedgt_script = Path(os.path.dirname(os.path.abspath(__file__))) / "src/vcf/vcf2bedGT.sh"
    
    cmd = [
        "bash", str(vcf2bedgt_script),
        f"results/{output_vcf}",
        f"results/LGvarhap1.bed",
        f"results/LGvarhap2.bed",
        f"results/LGvarall.bed"
    ]
    
    subprocess.run(cmd, check=True)
    logging.info("Complete!")

def plot(paf, output, minimum, distance, file):
    plot_script = Path(os.path.dirname(os.path.abspath(__file__))) / "src/pyFiles/plotpaf.py"
    cmd = [
        "python", str(plot_script),
        "-p", paf,
        "-o", output,
        "-m", str(minimum),
        "-d", str(distance),
        "-f", file
    ]
    subprocess.run(cmd, check=True)
    logging.info("Done!")

def pair(paf, length, output_prefix):
    pair_script = Path(os.path.dirname(os.path.abspath(__file__))) / "src/scripts/chrom_pair.sh"
    cmd = [
        "bash", str(pair_script),
        paf,
        str(length),
        output_prefix
    ]
    subprocess.run(cmd, check=True)

def reverse(paf, reverse, genome, new_genome):
    reverse_script = Path(os.path.dirname(os.path.abspath(__file__))) / "src/scripts/rc_chrom.sh"
    cmd = [
        "bash", str(reverse_script),
        paf,
        reverse,
        genome,
        new_genome
    ]
    subprocess.run(cmd, check=True)

    