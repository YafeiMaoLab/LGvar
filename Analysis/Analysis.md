# LGvar: Large-scale Genomic VARiation

## Quick Start
* [Assembly-based SV calling](#SVcalling)
    * [LGvar](#LGvar)
    * [PAV](#PAV)
    * [SyRI](#SyRI)
    * [SVIM-asm](#SVIM-asm)
* [Simulation test](#Simulation)
* [Benchmarking in human population datasets](#Human)
    * [Benchmarking in T2T-HG002 vs. T2T-CHM13](#T2T-HG002)
    * [Benchmarking in NA12878 pedigree vs. hg38](#NA12878)
    * [Inversion benchmark](#Inversion)
* [Benchmarking in cross-species datasets](#Cross-species)
    * [Macaque](#Macaque)
    * [Arabidopsis thaliana](#Tair)
* [Getting Help](#getting-help)
* [Citing LGvar](#citation)

## Assembly-based SV calling <a id="SVcalling"></a>
* **[LGvar](https://github.com/YafeiMaoLab/LGvar) <a id="LGvar"></a>**
    ```Bash
    /usr/bin/time -v LGVAR run \
        -r ${REF}.fa \
        -q1 ${QUERY}.hap1.fa \
        -q2 ${QUERY}.hap2.fa \
        -cp1 ${QUERY}.hap1.cp.tsv \
        -cp2 ${QUERY}.hap2.cp.tsv \
        -m sensitive \
        -dv 0.02 \
        -d 500 \
        -s ${QUERY}
    ```
* **[PAV](https://github.com/BeckLaboratory/pav) (v2.3.4) <a id="PAV"></a>**
    ```Bash
    #!/bin/bash
    /usr/bin/time -v singularity run --bind "$(pwd):$(pwd)" --bind $(pwd)/:/mnt/h1,$(pwd):/mnt/h1 pav_latest.sif -c 16
    ```
    `assemblies.tsv`: 
    | NAME  | HAP1     | HAP2  |
    | ------- | -------- | ---------- |
    | ${QUERY} | /mnt/h1/${QUERY}.hap1.fa    | /mnt/h1/${QUERY}.hap2.fa | 

    `config.json`:
    ```Bash
    {
	    "reference": ${REF}.fa
    }
    ```
* **[SyRI](https://github.com/schneebergerlab/syri) <a id="SyRI"></a>**
    
    `syri.sh`:
    ```Bash
    #!/bin/bash
    ## Minimap2 (same parameters with LGvar, PAV)
    minimap2 -t 12 -cx asm20 -m 10000 -z 10000,50 -r 50000 --end-bonus=100 -O 5,56 -E 4,1 -B 5 --secondary=no --eqx ${REF}.fa ${QUERY}.hap1.fa > align.hap1.paf
    minimap2 -t 12 -cx asm20 -m 10000 -z 10000,50 -r 50000 --end-bonus=100 -O 5,56 -E 4,1 -B 5 --secondary=no --eqx ${REF}.fa ${QUERY}.hap2.fa > align.hap2.paf

    ## Running SyRI
    for hap in hap1 hap2; do
        mkdir ${hap} && cd ${hap}
        tail -n +2 ${QUERY}.${hap}.cp.tsv | while read line; do
            ref=$(echo "$line" | cut -f1)
            queries=$(echo "$line" | cut -f2 | tr ',' ' ')
            for query in $queries; do
                mkdir -p "$query"
                cd "$query" || exit
                less -S ../../align.${hap}.paf | awk -v q="$query" -v c="$ref" '$1 == q && $6 == c' > "${ref}.paf"
                samtools faidx ${REF}.fa "$ref" > "$ref.fa"
                samtools faidx ${QUERY}.${hap}.fa "$query" > "$query.fa"
                nohup syri -c "$ref.paf" -r "$ref.fa" -q "$query.fa" -k -F P > "${ref}_syri.log" 2>/dev/null &
                cd ../
            done
        done
        cd ../
    done

    ## Merging VCF
    for hap in hap1 hap2; do
        cd ${hap}
        bgzip ./*/syri.vcf
        for dir in ./*;do
            tabix ${dir}/syri.vcf.gz
        done
        bcftools concat -o ${hap}.vcf.gz -Oz ./*/syri.vcf.gz -a 
        tabix -p vcf ${hap}.vcf.gz
        cd ..
    done

    bcftools merge -m none -Oz -o diploid.vcf.gz --force-samples ../hap1/hap1.vcf.gz ../hap2/hap2.vcf.gz
    zcat diploid.vcf.gz | awk 'BEGIN{FS=OFS="\t"} 
        /^##/ {print; next} 
        /^#CHROM/ { for(i=1;i<=10;i++) printf "%s%s", $i, (i==10?ORS:OFS); next } 
        { $10 = $10 "|" $11; NF=10; print }' | awk '$3 !~/SYN/ && $3 !~/AL/ && $3 !~/INV/' > diploid.cg.vcf
    bgzip diploid.cg.vcf && tabix diploid.cg.vcf.gz
    ```
    Using `/usr/bin/time -v bash syri.sh` to calculate resource consumption.
* **[SVIM-asm](https://github.com/eldariont/svim-asm) (v1.0.3) <a id="SVIM-asm"></a>**

    `svimasm.sh`:
    ```Bash
    #!/bin/bash
    ## Minimap2 (same parameters with LGvar, PAV, SyRI)
    minimap2 -x asm20 -m 10000 -z 10000,50 -r 50000 --end-bonus=100 -O 5,56 -E 4,1 -B 5 --secondary=no -a -t 12 --eqx -Y ${REF}.fa ${QUERY}.hap1.fa > hap1.align.sam
    samtools sort -m4G -@4 -o hap1.align.sorted.bam hap1.align.sam
    rm hap1.align.sam
    minimap2 -x asm20 -m 10000 -z 10000,50 -r 50000 --end-bonus=100 -O 5,56 -E 4,1 -B 5 --secondary=no -a -t 12 --eqx -Y ${REF}.fa ${QUERY}.hap2.fa > hap2.align.sam
    samtools sort -m4G -@4 -o hap2.align.sorted.bam hap2.align.sam
    rm hap2.align.sam
    samtools index hap1.align.sorted.bam
    samtools index hap2.align.sorted.bam

    ## Running SVIM-asm
    svim-asm diploid ./ hap1.align.sorted.bam hap2.align.sorted.bam ${REF}.fa
    ```
    Using `/usr/bin/time -v bash svimasm.sh` to calculate resource consumption.

## Simulation test <a id="Simulation"></a>
* We used hg38 chromosome 1 as the target genome, and employed **[VISOR](https://github.com/davidebolo1993/VISOR) (v1.1.2)** to construct the query genome by introducing SVs, INDELs and SNVs. 
    ```Bash
    ## Download SNP data
    curl -LO ftp://ftp.1000genomes.ebi.ac.uk:21/vol1/ftp/data_collections/1000_genomes_project/release/20181203_biallelic_SNV/ALL.chr1.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz

    ## Randomly select one sample, and add its SNV into hg38 chr1
    bcftools view -O b -o HG00134.snp.bcf -s HG00134 -m2 -M2 -c 1 -C 1 ALL.chr1.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz
    bcftools index HG00134.snp.bcf
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%SAMPLE=%GT]\n' HG00134.snp.bcf | grep "1|0" | awk 'OFS=FS="\t"''{print $1, ($2 -1), $2, "SNP", $4, "0"}' > HG00134.snps.h1.bed

    ## Then use VISOR to integrate this SNVs
    VISOR HACk -g hg38_chr1.fa -b HG00134.snps.h1.bed -o chr1withsnp
    cut -f1,2 hg38_chr1.fa.fai > chr1.dim.tsv

    ## Generate some Large SVs
    Rscript "/home/ffzhou/sv/VISOR/scripts/randomregion.r" -d chr1.dim.tsv -n 500 -l 20000 -s 10000 -v 'deletion,insertion,inversion,tandem duplication,translocation cut-paste' -r '30:30:30:5:5' | sortBed > chr1.sv.bed
    ## Generate some small SVs
    Rscript "/home/ffzhou/sv/VISOR/scripts/randomregion.r" -d chr1.dim.tsv -n 100 -l 30 -s 10 -v 'deletion,insertion' -r '50:50' | sortBed > chr1.ssv.bed
    cat chr1.ssv.bed chr1.sv.bed > chr1.variation.bed

    ## Then integrate these SVs
    VISOR HACk -g chr1withsnp/h1.fa -b chr1.variation.bed -o chr1
    ```
    The simulated query genome is stored in `chr1/h1.fa`. 
    We generated three query genomes by adjusting `-n`, `-l` and `-s`. 
    * Simulation 1: `-l` 10000, `-s` 5000, no SNVs and INDELs
    * Simulation 2: `-l` 20000, `-s` 10000, 81,145 SNVs and 100 INDELs (`-l` 30, `-s` 10)
    * Simulation 2: `-l` 50000, `-s` 25000, 153,210 SNVs and 200 INDELs (`-l` 50, `-s` 10)
    Details regarding sequence divergence are presented in `Supplementary Figures 3 and 4`.
* **SV calling:** We used the same command as in the [Assembly-based SV calling](#SVcalling) section, with `-dv` set to `0.05` in the LGvar command.
* **Metrix calculating:** We employed **[bedtools](https://github.com/arq5x/bedtools2) (v2.31.1)** intersect to obtain overlapping true-positive variants.
    ```Bash
    for tool in LGvar PAV SyRI SVIM-asm;do
        for sv in ins del inv;do
            bedtools intersect -a chr1_${sv}.bed -b del.bed -wa -wb -f 0.5 -F 0.5 | sort -k2,2n | cut -f1-3 | uniq >> tp.bed
        done
    done
    ```
## Benchmarking in human population datasets <a id="Human"></a>

### 1. Benchmaking in T2T-HG002 vs. T2T-CHM13 <a id="T2T-HG002"></a>
* **dataset:** https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/analysis/NIST_HG002_DraftBenchmark_defrabbV0.019-20241113/
* **SV truthset:** [CHM13v2.0_HG2-T2TQ100-V1.1_stvar.vcf.gz](https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/analysis/NIST_HG002_DraftBenchmark_defrabbV0.019-20241113/CHM13v2.0_HG2-T2TQ100-V1.1_stvar.vcf.gz)
* **SV high-condifent region:** [CHM13v2.0_HG2-T2TQ100-V1.1_stvar.benchmark.bed](https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/analysis/NIST_HG002_DraftBenchmark_defrabbV0.019-20241113/CHM13v2.0_HG2-T2TQ100-V1.1_stvar.benchmark.bed)
* **Small SV truthset:** [CHM13v2.0_HG2-T2TQ100-V1.1_smvar.vcf.gz](https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/analysis/NIST_HG002_DraftBenchmark_defrabbV0.019-20241113/CHM13v2.0_HG2-T2TQ100-V1.1_smvar.vcf.gz)
* **Small SV high-confident region:** [CHM13v2.0_HG2-T2TQ100-V1.1_smvar.benchmark.bed](https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/analysis/NIST_HG002_DraftBenchmark_defrabbV0.019-20241113/CHM13v2.0_HG2-T2TQ100-V1.1_smvar.benchmark.bed)
* **SV calling:** We used the same command as in the [Assembly-based SV calling](#SVcalling) section, with `-dv` set to `0.05` in the LGvar command.
* Then use [truvari](https://github.com/acenglish/truvari) and [hap.py](https://github.com/Illumina/hap.py) to do benchmark test.
    ```Bash
    ##1. SV benchmark
    for tool in LGvar PAV SyRI SVIM-asm;do
        truvari bench --pctseq 0.8 --pctsize 0.7 -r 1000 --pick multi --sizemax 200000 -b CHM13v2.0_HG2-T2TQ100-V1.1_stvar.vcf.gz -c ${sp}.${tool}.vcf.gz --reference chm13v2.0.fa --includebed CHM13v2.0_HG2-T2TQ100-V1.1_stvar.benchmark.bed -o bench-0.80.7-1000r-multi/
    done

    ##2. Small SV benchmark
    for tool in LGvar PAV SyRI;do
        singularity exec hap.py_latest.sif /opt/hap.py/bin/hap.py CHM13v2.0_HG2-T2TQ100-V1.1_smvar.vcf.gz ${sp}.${tool}.ssv.vcf.gz -f CHM13v2.0_HG2-T2TQ100-V1.1_smvar.benchmark.bed -r chm13v2.0.fa --threads=64 -o sSV
    done
    ```

    The result are `${sp}-bench-0.80.7-1000r-multi/summary.json` and `sSV.summary.csv`.

    <details><summary><b>Resource consumption</b></summary>
    <center>

    | Tool  | CPU time (h)    | Peak memory (GB)  | 
    | ------- | -------- | ---------- | 
    | **LGvar** | 2.5590    | 40.6697 | 
    | PAV | 43.8794      | 40.6562 | 
    | SyRI | 2.0370     | 40.8745 |
    | SVIM-asm | 3.2235 | 40.8219 |

    </details>
    </center>

### 2. Benchmaking in NA12878 pedigree vs. hg38 <a id="NA12878"></a>
* **dataset:** https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Datasets
* **SV truthset:** merged_hg38.svs.sort.oa.vcf.gz
* **Small SV truthset:** CEPH1463.GRCh38.family-truthset.ov.vcf.gz
    ```Bash
    for sp in NA12877 NA12878 NA12879 NA12881 NA12882 NA12885 NA12886;do
        ##SV
        bcftools view -c1 -s ${sp} -Oz -o ${sp}.SV.vcf.gz merged_hg38.svs.sort.oa.vcf.gz
        ##Small SV
        bcftools view -c1 -s ${sp} -Oz -o ${sp}.ssv.vcf.gz CEPH1463.GRCh38.family-truthset.ov.vcf.gz
    done
    ```
* Using [dipcall](https://github.com/lh3/dipcall) to generate high confident region for each sample.
    ```Bash
    for sp in NA12877 NA12878 NA12879 NA12881 NA12882 NA12885 NA12886;do
        ~/software/dipcall.kit/run-dipcall ${sp} GRCh38.no_alt_analysis_set.fa ${sp}.hap2.fa ${sp}.hap1.fa > ${sp}.mak
        make -j2 -f ${sp}.mak
    done
    ```
    Using ${sp}.dip.bed as high confident region to benchmark.
* **SV calling:** We used the same command as in the [Assembly-based SV calling](#SVcalling) section.

* Then use [truvari](https://github.com/acenglish/truvari) and [hap.py](https://github.com/Illumina/hap.py) to do benchmark test.
    ```Bash
    ##1. SV benchmark
    for sp in NA12877 NA12878 NA12879 NA12881 NA12882 NA12885 NA12886;do
        for tool in LGvar PAV SyRI SVIM-asm;do
            truvari bench --pctseq 0.8 --pctsize 0.8 -r 1000 --pick multi --sizemax 200000 -b ${sp}.SV.vcf.gz -c ${sp}.${too}.vcf.gz --reference GRCh38.no_alt_analysis_set.fa --includebed ${sp}.dip.bed -o ${sp}-bench-0.80.7-1000r-multi/
        done
    done

    ##2. Small SV benchmark
    for sp in NA12877 NA12878 NA12879 NA12881 NA12882 NA12885 NA12886;do
        for tool in LGvar PAV SyRI;do
            singularity exec hap.py_latest.sif /opt/hap.py/bin/hap.py ${sp}.ssv.vcf.gz ${sp}.${too}.ssv.vcf.gz -f ${sp}.dip.bed -r GRCh38.no_alt_analysis_set.fa --threads=64 -o sSV
        done
    done
    ```
    The result are `${sp}-bench-0.80.7-1000r-multi/summary.json` and `sSV.summary.csv`.

    <details><summary><b>Resource consumption (Average results of 7 samples)</b></summary>
    <center>

    | Tool  | CPU time (h)    | Peak memory (GB)  | 
    | ------- | -------- | ---------- | 
    | **LGvar** | 2.5181    | 34.168 | 
    | PAV | 75.9507      | 40.2392 | 
    | SyRI | 1.9212     | 34.1936 |
    | SVIM-asm | 3.4726 | 34.3886 |

    </details>
    </center>

### 3. Inversion benchmark <a id="Inversion"></a>
* **Inversion calling:**  We used the same command as in the [Assembly-based SV calling](#SVcalling) section. 
* **benchmark:** we conducted a targeted benchmark across five human genomes (HG002, HG00733, HG02818, HG03486, and NA19240) provided in the GitHub (https://github.com/jamesc99/INV-Benchmark/tree/main/data_zenodo/bed)
## Benchmarking in cross-species datasets <a id="Cross-species"></a>
### 1. Macaque <a id="Macaque"></a>
* **Inversion calling:**  We used the same command as in the [Assembly-based SV calling](#SVcalling) section. with minimap2 parameters set to `-a -x asm20 --eqx --cs -K 500M -k 15 -m 10 -A 1 -B 2 -O 2,12 -n 2 -g 100 -r 200,100000 --secondary=no -s 1000 -o aln.sam` in the LGvar, SyRI and SVIM-asm command.
* **benchmark:** We used a curated set of large-scale genomic rearrangements reported by [Zhang et al.](https://www.nature.com/articles/s41586-025-08596-w) as the benchmark reference. Overlap between tool-derived variant calls and the curated events was assessed using ` bedtools intersect -a truthset.bed -b query.bed -f 0.5 -F 0.5`, requiring at least 50% reciprocal overlap to define a true positive.

### 2. Arabidopsis thaliana <a id="Tair"></a>
* We used TAIR10 [GCA_000001735.3](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001735.3/) as the reference and the Ler [GCA_900660825.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_900660825.1/) as the query.

## Getting help <a id="getting-help"></a>

* If you have questions about benmchmark test, please contact zhoufeifei@sjtu.edu.cn.

## Citing LGvar <a id="citation"></a>

* If you use LGvar in your research, please cite our repository (and forthcoming paper): (https://github.com/YafeiMaoLab/LGvar)
