# LGvar: Large-scale Genomic VARiation
[![PyPI version](https://badge.fury.io/py/lgvar.svg)](https://pypi.org/project/lgvar/)
[![Conda version](https://img.shields.io/badge/conda-v1.4.0-green)](https://anaconda.org/channels/zhoufeifei/packages/lgvar/overview)  [![Appatiner Image Version](https://img.shields.io/badge/singularity-v1.4.0-blue)](https://cloud.sylabs.io/library/feifeizhou/tool/lgvar) [![License](https://img.shields.io/badge/license-MIT-yellow)](https://github.com/YafeiMaoLab/LGvar/blob/main/LICENSE)  [![release](https://img.shields.io/badge/releases-July%202026-purple)](https://github.com/YafeiMaoLab/LGvar/releases/tag/v1.4.0)

## Quick Start
* [Why LGvar?](#WhyLGvar)
* [Installation & Examples](#installation)
* [Usage](#Usage-Examples)
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

**If you are using an `osx-arm64` platform, you may install LGVAR via `Option 3` with singularity.**

**If you are using an `linux-64` platform, we recommend you install LGVAR via `Option 1`.**

**Option 1: Pip (Recommended)**
```Bash
# You can install LGvar from pypi
pip install lgvar==1.4.0
# Then install the corresponding dependencies
LGVAR setup --spawn
conda activate LGVAR
```

Use provided example datasets to test, the `example/` folder is stored under the LGVAR installation directory.
You can run `pip show lgvar` to retrieve this path, which typically follows the pattern: `*/anaconda3/lib/pythonX.X/site-packages/lgvar/examples/`.
```Bash
# Test
LGVAR run \  
    -r */anaconda3/lib/pythonX.X/site-packages/lgvar/examples/genome/test.ref.fa \ 
    -q1 */anaconda3/lib/pythonX.X/site-packages/lgvar/examples/genome/test.hap1.fa \  
    -q2 */anaconda3/lib/pythonX.X/site-packages/lgvar/examples/genome/test.hap2.fa \  
    -p1 */anaconda3/lib/pythonX.X/site-packages/lgvar/examples/align/align_hap1.paf \  
    -p2 */anaconda3/lib/pythonX.X/site-packages/lgvar/examples/align/align_hap2.paf \
    -cp1 */anaconda3/lib/pythonX.X/site-packages/lgvar/examples/data/hap1.tsv \
    -cp2 */anaconda3/lib/pythonX.X/site-packages/lgvar/examples/data/hap2.tsv \
    -m sensitive \
    -s test
```
**Option 2: Conda**

```Bash
# Or you can directly create a new environment
# Clone the repository  
cd ${workdir}
git clone https://github.com/YafeiMaoLab/LGvar.git 
# Create and activate environment  
conda env create -f LGvar/environment.yml  
conda activate LGVAR

# Use
${workdir}/LGvar/LGVAR run \  
    -r ${workdir}/examples/genome/test.ref.fa \ 
    -q1 ${workdir}/examples/genome/test.hap1.fa \  
    -q2 ${workdir}/examples/genome/test.hap2.fa \  
    -p1 ${workdir}/examples/align/align_hap1.paf \  
    -p2 ${workdir}/examples/align/align_hap2.paf \
    -cp1 ${workdir}/examples/data/hap1.tsv \
    -cp2 ${workdir}/examples/data/hap2.tsv \
    -m sensitive \
    -s test
```
**Option 3: Singularity**
```Bash
# Here is our pre-built image:
singularity pull library://feifeizhou/tool/lgvar:v1.4.0 
# Execution example:  
singularity exec lgvar_v1.4.0.sif /LGVAR/LGVAR run \
    -r /LGVAR/examples/genome/test.ref.fa \ 
    -q1 /LGVAR/examples/genome/test.hap1.fa \  
    -q2 /LGVAR/examples/genome/test.hap2.fa \  
    -p1 /LGVAR/examples/align/align_hap1.paf \  
    -p2 /LGVAR/examples/align/align_hap2.paf \
    -cp1 /LGVAR/examples/data/hap1.tsv \
    -cp2 /LGVAR/examples/data/hap2.tsv \
    -m sensitive \
    -s test
```

## Usage<a id="Usage-Examples"></a>

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
* We benchmarked LGvar in many genomes include simulated, population and cross-species data (see forthcoming paper). See [Analysis/Analysis.md](https://github.com/YafeiMaoLab/LGvar/blob/main/Analysis/Analysis.md):


## Getting help <a id="getting-help"></a>

* If you encounter any issues or have questions about specific parameters, please [Open an Issue](https://github.com/YafeiMaoLab/LGvar) or contact zhoufeifei@sjtu.edu.cn.


## Citing LGvar <a id="citation"></a>

* If you use LGvar in your research, please cite our repository (and forthcoming paper): (https://github.com/YafeiMaoLab/LGvar)
