# mtDNApipe
**A robust pipeline for ultra-sensitive detection of low-frequency mitochondrial DNA variants from paired-end NGS data.**

---

## Overview

`mtDNApipe` is a comprehensive bioinformatics pipeline designed for the high-fidelity detection of mitochondrial DNA (mtDNA) heteroplasmy from Next-Generation Sequencing (NGS) data. The workflow starts with raw paired-end FASTQ files and implements a series of rigorous steps to identify low-frequency variants with high confidence.

Key features of the pipeline include:
- **Stringent Quality Control**: A custom sliding-window filter followed by `fastp` for robust read trimming.
- **Advanced Error Correction**: A custom C++ tool to correct sequencing errors in overlapping read pairs.
- **Best-Practice Alignment**: Standard BWA-MEM alignment followed by GATK 3.8 Indel Realignment to minimize alignment artifacts.
- **Custom Deduplication**: A sophisticated consensus-based deduplication strategy (`con_rmdup_v1.5.py`) that significantly outperforms standard methods for deep sequencing data.
- **Rigorous Variant Calling**: A final calling script that filters for terminal base-pair artifacts, oxidative damage signatures, and known problematic repeat regions in the mitochondrial genome.

---

## Installation

We strongly recommend using `Conda` to manage all software and package dependencies.

### Step 1: Clone the Repository
```bash
git clone https://github.com/FMMU-Xing-Lab/mtDNApipe.git
cd mtDNApipe
```

### Step 2: Create and Activate the Conda Environment

This single command will install all required tools (BWA, fastp, GATK 3.8, etc.) and Python packages (pysam, pandas, tqdm).

```bash
conda env create -f environment.yml
conda activate mtDNApipe
 ```
    
    

### Step 3: Compile C++ Helper Tool
The pipeline uses a custom C++ tool for correcting read overlaps. Compile it using the provided Makefile.
```bash
cd scripts/overlap_corrector/
make
```


### Step 4:  Prepare Reference Genome
You will need a mitochondrial reference genome in FASTA format. You can download one from resources like UCSC or Ensembl. ***You must index it for BWA and Samtools before the first run.***
```bash
# This is a one-time setup for your reference genome
bwa index /path/to/your/hg19_mt.fa
samtools faidx /path/to/your/hg19_mt.fa
# You will also need to create a dictionary file for GATK
java -jar $(which picard.jar) CreateSequenceDictionary R=/path/to/your/hg19_mt.fa O=/path/to/your/hg19_mt.dict
```
   

## 3. Usage

The entire workflow is orchestrated by the `run_pipeline.sh` script.

### Command
```bash
bash run_pipeline.sh \
  -i <sample_id> \
  -1 <path/to/read1.fq.gz> \
  -2 <path/to/read2.fq.gz> \
  -r <path/to/hg19_mt.fa> \
  -o <path/to/output_dir> \
  -t <threads>
```

### Arguments
| Flag | Description | Required | Default |
| ----- | ----- | ----- | ----- |
| `-i` | Sample identifier (e.g., `Sample01`) | Yes | - |
| `-1` | Path to Read 1 FASTQ file (`.fq.gz`)| Yes | - |
| `-2` | Path to Read 2 FASTQ file (`.fq.gz`) | Yes | - |
| `-r` | Path to the mitochondrial reference genome FASTA file.	 | Yes | - |
| `-o` | Path to the main output directory. | Yes | - |
| `-t` | Number of threads to use for BWA alignment. | No | 8 |
| `-h` | Display the help message. | No | - |

## 4. Input and Output

### Input Files

Paired-end FASTQ files: Gzipped (`_R1.fq.gz`, `_R2.fq.gz`).
Reference Genome: A FASTA file containing the mitochondrial sequence (`hg19_mt.fa`), pre-indexed.

### Output Directory Structure
All results for a given sample will be organized in a dedicated sub-directory: `<output_dir>/<sample_id>/.`

### Key Output Files
| Flag | Description |
| ----- | ----- |
| `SAMPLE_ID.conrmdup.v5.end10.mutations.txt` | **Final high-confidence mutation list.** A tab-separated file with position, ref, alt, VAF, and depth. |
| `SAMPLE_ID.mt.no.softclip.bam` | The final processed, consensus-deduplicated BAM file used for variant calling. Ready for IGV inspection. |
| `SAMPLE_ID_analysis.log` | A comprehensive log file detailing every step of the pipeline execution. Essential for troubleshooting. |
| `SAMPLE_ID_fastp.html` | An interactive quality report from `fastp`. |

## Example
To ensure your installation is working correctly, run the provided test case. It uses a small subset of reads to execute the full pipeline.
```bash
# Make sure you have activated the conda environment
# conda activate mtDNApipe

# Run the test script. You must provide the path to your reference genome.
bash example/run_test.sh /path/to/your/hg19_mt.fa
```
The test should complete in a few minutes and generate an output directory named `test_output/`. You can inspect the final mutation file at `test_output/test_sample/test_sample.conrmdup.v5.end10.mutations.txt`.

## Contributing & Issues
We welcome bug reports, feature requests, and pull requests. Please feel free to:
Open an [Issue](https://github.com/FMMU-Xing-Lab/mtDNApipe/issues) to report a bug or suggest a feature.
    

## License
This project is licensed under the [MIT License](https://github.com/FMMU-Xing-Lab/mtDNApipe/blob/main/LICENSE).
