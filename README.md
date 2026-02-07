# LGvar: Large-scale Genomic VARiation
[![Conda version](https://img.shields.io/badge/conda-v1.1.0-green)](https://anaconda.org/channels/zhoufeifei/packages/lgvar/overview)  [![Appatiner Image Version](https://img.shields.io/badge/singularity-v1.1.0-blue)](https://cloud.sylabs.io/library/feifeizhou/tool/lgvar) [![License](https://img.shields.io/badge/license-MIT-yellow)](https://github.com/YafeiMaoLab/LGvar/blob/main/LICENSE)  [![release](https://img.shields.io/badge/releases-Feb%202026-purple)](https://github.com/YafeiMaoLab/LGvar/releases/tag/v1.1.0)

**LGvar** is a comprehensive toolset designed for large-scale structural variant (SV) detection based on genome assemblies. It excels in **cross-species variant identification**, particularly demonstrating superior performance in **inversion detection** compared to existing tools.

##

**Quick Start**
* [Installation](#installation)
* [Tool Modules](#tool-modules)
* [Usage & Examples](#Usage-Examples)
* [Output Files](#output-files)
* [Detailed Parameters](#detailed-parameters)
* [Getting Help](#getting-help)
* [Citing LGvar](#citation)

##

**Installation <a id="installation"></a>**

**Option 1: Conda (Recommended)**

```Bash
# You can install LGvar from anaconda.org
conda create -n lgvar python=3.9 -y
conda activate lgvar
conda install -c conda-forge -c bioconda zhoufeifei::lgvar
```
```Bash
# Or you can directly create a new environment
# Clone the repository  
git clone https://github.com/YafeiMaoLab/LGvar.git 
cd LGvar
# Create and activate environment  
conda env create -f env.yml  
conda activate LGvar
```
**Option 2: Singularity**
```Bash
# If conda is not available, use our pre-built images:
singularity pull --arch amd64 library://feifeizhou/tool/lgvar:v1.1.0 
# Execution example:  
singularity exec lgvar_v1.1.0.sif /LGVAR/LGVAR [subcommands]
```
##

**Tool Modules <a id="tool-modules"></a>**

LGvar provides several subcommands to streamline the genomic analysis pipeline:

| Command | Function |
| :---- | :---- |
| **run** | **Run SV identification.** Detects SDRs, DELs, INSs, INVs, SDRs, DUPs, INDELs, SNVs etc. |
| **plot** | **Syntenic Plot.** Visualizes the alignment between genomes. Blue regions represent syntenic blocks, while orange regions represent inversions. |
| **pair** | **Homology Mapping.** Generates homologous chromosome pair files. |
| **reverse** | **RC-Check.** Identifies and reverse-complements misoriented chromosomes. |

##

**Usage & Examples <a id="Usage-Examples"></a>**

**1\. Prepare Query Genome (reverse)**

Before calling SVs, check if query chromosomes need reverse-complementing to ensure accuracy:
```Bash
# Recommended: align first then check orientation  
LGVAR reverse -p align.paf -r rev.txt -g original.genome.fa -n new.genome.fa
# Or use the script
bash /src/scripts/rc_chrom.sh align.paf rev.txt original.genome.fa new.genome.fa
```
**2\. Run Test Example**

The test data is located in examples/. We'll use Human (CHM13) as reference and Chimpanzee (PTR) as query.

```Bash

# 1. Setup working directory  
mkdir LGvar_work && cd LGvar_work
# 2. Decompress test genomes  
cd ../examples/genome  
gunzip *.gz  
cd ../LGvar_work
# 3. Execution (Example with existing PAF alignments)  
LGVAR run \  
  -r /examples/genome/chr21.chm13.fa \ 
  -q1 /examples/genome/chr21.ptr.hap1.fa \  
  -q2 /examples/genome/chr21.ptr.hap2.fa \  
  -p1 /examples/align/align_hap1.paf \  
  -p2 /examples/align/align_hap2.paf \
  -cp1 /examples/data/PTR_hap1_pairs.tsv \
  -cp2 /examples/data/PTR_hap2_pairs.tsv \
  -cen /examples/data/chm13_cen.tsv \
  -telo /examples/data/chm13_telo.tsv \
  -m cts \
  -s PTR
```
##

**Output Files <a id="output-files"></a>**

Results are saved in the **${work_dir}/results** folder.
* **sortLGvar\_all.vcf**: Combined results of SNVs, INDELs and SVs for both haplotypes.
* **LGvarhap1(2).vcf**: Haplotype-specific SV details.
* **LGvar.bed**: Merged variant regions (Hap1 \+ Hap2).  
* **LGvarhap1(2).bed**: Haplotype-specific SV details:

| Label | Description |
| :---- | :---- |
| **SNVs** | Single Nucleotide Variation |
| **INDELs** | Small insertions and deletions (Length < 50bp) |
| **DELs/INSs** | Deletions and Insertions (filtered at 50bp threshold) |
| **INVs/INV-INV** | Simple and Nested Inversions |
| **DUP** | Duplications |
| **TRANS** | Translocations |
| **SDRs** | Structure Divergent Regions |

##

**Detailed Parameters <a id="detailed-parameters"></a>**

### **run subcommand**
```Bash
╔════════════════════════════════════════════════╗
║                                                ║
║   ██╗      ██████╗ ██╗   ██╗ █████╗ ██████╗    ║
║   ██║     ██╔════╝ ██║   ██║██╔══██╗██╔══██╗   ║
║   ██║     ██║  ███╗██║   ██║███████║██████╔╝   ║
║   ██║     ██║   ██║║██  ██╔╝██╔══██║██╔══██╗   ║
║   ███████╗╚██████╔╝╚╗ ██ ╔╝ ██║  ██║██║  ██║   ║
║   ╚══════╝ ╚═════╝  ╚════╝  ╚═╝  ╚═╝╚═╝  ╚═╝   ║
║                                                ║
║                    L G V A R                   ║
╚════════════════════════════════════════════════╝

Large-Scale Genetic VARiation caller

usage: LGVAR run [-h] -r REF -q1 HAP1 [-q2 HAP2] [-p1 PAF1] [-p2 PAF2] -cp1 PAIRS1 [-cp2 PAIRS2]
                 [-c CLUSTER] [-dl DELLENGTH] [-inv INVCLUSTER] [-cen CENTROMERE] [-telo TELOMERE]
                 -m {ctn,cts} [-s SAMPLE_NAME] [-v VARIANT_TYPE] [-d DISTANCE] [-t THREADS]
                 [-k CHUNK_SIZE] [-f FRACTION] [-mdist MAX_DISTANCE] [-sdist SMALL_DISTANCE]
                 [-sim SIMILARITY_THRESHOLD]

optional arguments:
  -h, --help                    show this help message and exit

Input Files:
  -r, --ref                     Reference genome for variants calling
  -q1, --hap1                   One query genome (Which is one haplotype of one species genome and
                                needs to be scaffolded to chromosome level using RagTag to ensure
                                better genome quality for more reliable variant calling results.)
  -q2, --hap2                   Another query genome (Which is another haplotype of the species
                                genome and needs to be scaffolded to chromosome level using RagTag
                                to ensure better genome quality for more reliable variant calling
                                results.)
  -p1, --paf1                   Alignment of haplotype1 (Which contains the CIGAR infomation)
                                [Recommend mapping tool: minimap2].
  -p2, --paf2                   Alignment for haplotype2 (Which contains the CIGAR infomation)
                                [Recommend mapping tool: minimap2].
  -cp1, --pairs1                Homologous chromsome pairs of query genome (hap1) and reference.
  -cp2, --pairs2                Homologous chromsome pairs of query genome (hap2) and reference.
  -c, --cluster                 Clustering parameter for filtering chaos alignments, where a smaller
                                value results in a stricter filter [200000].
  -dl, --dellength              A desired deletion length for alignments, where a larger value
                                enforces a stricter filter [300000].
  -inv, --invcluster            Clustering parameter for inversion calling [700000].
  -cen, --centromere            A centromere file which is used to filter out the alignment of
                                complex regions that may not be well aligned [False].
  -telo, --telomere             A telomere file which is used to filter out the alignment of complex
                                regions that may not be well aligned [False].
  -m, --mode                    Analysis mode: ctn (do not remove centromere and telomere
                                alignments) or cts (remove) [ctn].
  -s, --sample_name             Sample name used to generate vcf.
  -v, --variant_type            Process line of each thread in inversion re-identification. [all]

Personalized arguments:
  -d, --distance                Parameters used to identify INS and DEL: Variations where the
                                distance between the start and end positions of the REF or QUERY is
                                less than (or equal with) d will be identified as SVs. A lower value
                                indicates stricter criteria for SV identification [0bp].
  -t, --threads                 Multi threads for inversion re-identification from SDR. [4]
  -k, --chunk_size              Process line of each thread in inversion re-identification. [50]
  -f, --fraction                Specifies the minimum alignment threshold for the realignment step.
                                The input integer is scaled by a factor of 0.1 to determine the
                                actual percentage of bases required to align. For example, setting
                                --fraction 5 corresponds to 50%. [5]

Additional arguments:
  -mdist, --max_distance        Max reference distance for two allele to merge [500bp].
  -sdist, --small_distance      Max reference distance for SSV (small variants) merge [10bp].
  -sim, --similarity_threshold  The similarity of variants, used to merge the variants of two
                                haplotypes [0.8].
```
                                
### **plot subcommand**

You can visualize whole-genome alignments using PAF and PAIR files.
```Bash
LGVAR plot -p align.paf -f pair.tsv -o align.pdf
```
This is an example using T2T-CHM13 as reference and chimpanzee hap1 as query genome.
![https://github.com/YafeiMaoLab/LGvar/images/alignment.png](https://github.com/YafeiMaoLab/LGvar/blob/main/images/alignment.png)
##

**Getting Help <a id="getting-help"></a>**

If you encounter any issues or have questions about specific parameters, please [Open an Issue](https://github.com/YafeiMaoLab/LGvar).

##

**Citing LGvar <a id="citation"></a>**

If you use LGvar in your research, please cite our repository (and forthcoming paper):

(https://github.com/YafeiMaoLab/LGvar)

