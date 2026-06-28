# LGvar: Large-scale Genomic VARiation
[![Conda version](https://img.shields.io/badge/conda-v1.3.0-green)](https://anaconda.org/channels/zhoufeifei/packages/lgvar/overview)  [![Appatiner Image Version](https://img.shields.io/badge/singularity-v1.3.0-blue)](https://cloud.sylabs.io/library/feifeizhou/tool/lgvar) [![License](https://img.shields.io/badge/license-MIT-yellow)](https://github.com/YafeiMaoLab/LGvar/blob/main/LICENSE)  [![release](https://img.shields.io/badge/releases-June%202026-purple)](https://github.com/YafeiMaoLab/LGvar/releases/tag/v1.3.0)

## Quick Start
* [Why LGvar?](#WhyLGvar)
* [Installation](#installation)
* [Usage & Examples](#Usage-Examples)
    * [Reverse-complement the genome](#reverse)
    * [Generate orthologous chromosome pairs](#pair)
    * [Whole genome SV identification](#run)
    * [Output files](#output)
    * [Plot whole genome alignment synteny](#plot)
* [Benchmarking](#benchmark)
* [Getting Help](#getting-help)
* [Citing LGvar](#citation)

## Why LGvar <a id="WhyLGvar"></a>
**LGvar** is a comprehensive toolset designed for large-scale structural variant (SV) detection based on genome assemblies. It excels in **cross-species variant identification**, particularly demonstrating superior performance in **inversion detection** compared to existing tools.

- **Easily use:** Only FASTA and chromosome pairs (.tsv) are needed for input (users can alternatively provide pre-aligned PAF files).
- **Comprehensive and fast:** LGvar detects not only insertions, deletions, inversions, translocations, duplications, INDELs and SNVs, but also detects structurally divergent regions (SDRs) and nested-inversions (INV-INV).
- **More efficient in detection of inversions:** LGvar performs realignment within SDRs to detect inversions missed in initial alignment.

## Installation <a id="installation"></a>

**Option 1: Conda (Recommended)**

```Bash
# You can install LGvar from anaconda.org
conda install -c conda-forge -c bioconda lgvar
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
singularity pull library://feifeizhou/tool/lgvar:v1.3.0 
# Execution example:  
singularity exec lgvar_v1.3.0.sif /LGVAR/LGVAR [subcommands]
```

## Usage & Examples <a id="Usage-Examples"></a>

* LGvar provides several subcommands to streamline the genomic analysis pipeline:
    <details>
    <summary><b>Subcommands List (click to expand)</b></summary>

    | Command | Function |
    | :------ | :------- |
    | **run** | **Run SV identification.** Detects SDRs, DELs, INSs, INVs, DUPs, INDELs, SNVs etc. |
    | **plot** | **Syntenic Plot.** Visualizes the alignment between genomes. Blue regions represent syntenic blocks, while orange regions represent inversions. |
    | **pair** | **Homology Mapping.** Generates homologous chromosome pair files. |
    | **reverse** | **RC-Check.** Identifies and reverse-complements misoriented chromosomes. |

    </details>


**1\. Reverse-complement the genome <a id="reverse"></a>**

* Before calling SVs, using `LGvar reverse` to check if query chromosomes need reverse-complementing to ensure inversion calling accuracy, this step identifies a chromosome as needing reverse complementation by checking if the total length of alignments mapped to the reverse strand exceeds `half` the chromosome’s total length in the PAF file.
    ```Bash
    # Recommended: align first then check orientation  
    LGVAR reverse -p align.paf -r reverse.txt -g original.genome.fa -n new.genome.fa
    # Or use the script
    bash /src/scripts/rc_chrom.sh align.paf reverse.txt original.genome.fa new.genome.fa
    #Then use new.genome.fa as query genome to do SV calling
    ```
    `reverse.txt` contains the chromosomes that need to be reverse-complemented.

**2\. Generate orthologous chromosome pairs <a id="pair"></a>**

* LGvar identifies variants between each pair of orthologous chromosomes. If you do not have a tab-delimited TSV file containing `ref   query` columns, you can generate the orthologous chromosome TSV file using the `LGvar pair` command.
    ```Bash
    LGVAR pair -p align.paf -l 1000000 -o pair
    ```
    where the `-l` parameter specifies the minimum alignment length (in base pairs) required to define a reliable orthologous chromosome pair; only alignments meeting or exceeding this `1Mbp` threshold will be used to determine the corresponding relationships between chromosomes. It is recommended to increase this value when analyzing two species with a large evolutionary divergence (e.g. human vs. lemur, `-l` can be set to `5Mbp`).

**3\. Whole genome SV identification  <a id="run"></a>**

* Once you have prepared the reference genome, query genome and orthologous chromosome pair files, you can use `LGVAR run` to perform genome-wide variant detection. You may optionally provide files containing `centromere and telomere` coordinates of the reference genome to mask highly repetitive regions, thus improving running efficiency.
    ```Bash
    LGVAR run \
        -r ref.fa \
        -q1 query.hap1.fa \
        -q2 query.hap2.fa \
        -cp1 query.hap1.pair.tsv \
        -cp2 query.hap2.pair.tsv \
        -m sensitive \
        -s query
    ```
    Centromere and telomere coordinate files for `T2T-CHM13` are available under the `examples/data` directory, namely `chm13_cen.tsv` and `chm13_telo.tsv`.

* Here is an running example. The test data is located in `examples/`. We'll use Human chromosome 21(T2T-CHM13) as reference and Chimpanzee (PTR) as query.
    ```Bash
    # 1. Setup working directory  
    mkdir LGvar_work
    # 2. Decompress test genomes  
    cd examples/genome && gunzip *.gz  
    cd LGvar_work
    # 3. Execution (Example with existing PAF alignments)  
    LGVAR run \  
        -r /examples/genome/chr21.chm13.fa \ 
        -q1 /examples/genome/chr21.ptr.hap1.fa \  
        -q2 /examples/genome/chr21.ptr.hap2.fa \  
        -p1 /examples/align/align_hap1.paf \  
        -p2 /examples/align/align_hap2.paf \
        -cp1 /examples/data/PTR_hap1_pairs.tsv \
        -cp2 /examples/data/PTR_hap2_pairs.tsv \
        -m sensitive \
        -s PTR
    ```

**4\. Output Files <a id="output"></a>**

* Results are saved in the `${work_dir}/results` folder:
    <details>
    <summary><b>results</b></summary>

    | File | Function |
    | :------ | :------- |
    | **sortLGvar_all.vcf** | Combined results of SNVs, INDELs and SVs for both haplotypes. |
    | **LGvarhap1(2).vcf** | Haplotype-specific SV details. |
    | **LGvarhap1(2).bed** | Haplotype-specific SV details (with SDRs, DUPs, TRANs). |
    | **LGvar.bed** | Merged variant regions (Hap1 + Hap2). |

    </details>

    | Label | Description |
    | :---- | :---- |
    | **SNVs** | Single Nucleotide Variation |
    | **INDELs** | Small insertions and deletions (Length < 50bp) |
    | **DELs/INSs** | Deletions and Insertions (filtered at 50bp threshold) |
    | **INVs/INV-INV** | Simple and Nested Inversions |
    | **DUP** | Duplications |
    | **TRANS** | Translocations |
    | **SDRs** | Structurally Divergent Regions |

**5\. Plot whole genome alignment synteny<a id="plot"></a>**

* You can visualize whole-genome alignments using `PAF` and `PAIR` files.
    ```Bash
    LGVAR plot -p align.paf -f pair.tsv -o align.pdf
    ```
    This is an example using T2T-CHM13 as reference and chimpanzee hap1 as query genome.
    ![https://github.com/YafeiMaoLab/LGvar/images/alignment.png](https://github.com/YafeiMaoLab/LGvar/blob/main/images/alignment.png)

## Benchmarking<a id="benchmark"></a>
* We benchmarked LGvar in many genomes include simulated, population and cross-species data (see forthcoming paper). See [Analysis.md](https://github.com/YafeiMaoLab/LGvar/blob/main/README.md):


## Getting help <a id="getting-help"></a>

* If you encounter any issues or have questions about specific parameters, please [Open an Issue](https://github.com/YafeiMaoLab/LGvar) or contact zhoufeifei@sjtu.edu.cn.


## Citing LGvar <a id="citation"></a>

* If you use LGvar in your research, please cite our repository (and forthcoming paper): (https://github.com/YafeiMaoLab/LGvar)
