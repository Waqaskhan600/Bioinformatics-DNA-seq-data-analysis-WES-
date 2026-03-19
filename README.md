# Bioinformatics-DNA-seq-data-analysis-WES-
# Basic Whole Exome Sequencing (WES) Pipeline

This repository contains a simple, linear Next-Generation Sequencing (NGS) pipeline using bash for paired-end FastQ data processing, optimized for Whole Exome Sequencing (WES).

## Pipeline Steps

The pipeline incorporates the following standard industry practices:
0. **Quality Control & Trimming:** `fastqc` assesses read quality and `trim_galore` removes sequencing adapters and low-quality sequences to minimize noise.
1. **Alignment:** Trimmed reads are aligned to the reference genome (`hg38.fa`) using `bwa mem`, then sorted via `samtools sort`.
2. **Duplicate Marking:** Eliminates PCR or optical duplicates using `gatk MarkDuplicates` to prevent biased variant calling.
3. **Variant Calling:** Scans the alignment to identify SNPs and Indels using `gatk HaplotypeCaller`.
4. **Variant Filtration:** Applies rigorous, fixed thresholds adapted specifically for high-coverage WES data (e.g., Depth, Quality by Depth, Mapping Quality, etc.) using `gatk VariantFiltration`.
5. **Quality Extraction:** Excludes heavily filtered calls, leaving only high-confidence (`PASS`) variants using `bcftools view`.
6. **CSV Conversion:** Converts the final `PASS` variants VCF into an Excel-compatible Comma-Separated Values (`.csv`) file using `gatk VariantsToTable` and basic `awk` formatting.

## Installation & Setup

We highly recommend using [Conda](https://docs.conda.io/en/latest/miniconda.html) to manage your bioinformatics tools.

### 1. Create a Conda Environment
Included in this repository is an `environment.yaml` file containing all the necessary dependencies (`fastqc`, `trim-galore`, `bwa`, `samtools`, `gatk4`, and `bcftools`).

Create and activate the environment quickly via:
```bash
conda env create -f environment.yaml
conda activate wes_pipeline
```

### 2. Obtain a Reference Genome (hg38)
You need to download and index a reference genome before running the scripts. The Broad Institute provides the standard `hg38` resource bundle.

Run the following to download and prepare your reference:
```bash
# Create a folder for your references
mkdir -p resources/reference
cd resources/reference

# Download the Fasta genome sequence
wget https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta -O hg38.fa

# Create a BWA Index (Warning: this step can take 2-3 hours)
bwa index hg38.fa

# Create a Fasta sequence dictionary for GATK
gatk CreateSequenceDictionary -R hg38.fa

# Create a samtools expected FASTA index (.fai)
samtools faidx hg38.fa
```
*Note: Make sure that `REFERENCE` inside `ngs_pipeline.sh` points to the absolute path of the `hg38.fa` file you just downloaded.*

## Usage

*Note: The pipeline includes built-in **Pre-Run Validation Checks** that will ensure you have correctly downloaded the reference genome and placed your input data before executing any bioinformatics tools. If anything is missing, the script will halt immediately and provide the exact commands you need to fix the issue!*

1. **Configure the pipeline:** Open `ngs_pipeline.sh` and define the `PROJECT_DIR` variable to point to your project's absolute path. If you are running the script directly from within your project directory, you can simply set it to `PROJECT_DIR="$PWD"`. All other directories (data, outputs, etc.) will branch off from this folder.
2. Place your raw input reads in the `data/` directory inside your project folder (e.g., `sample_1_R1.fastq.gz`, `sample_1_R2.fastq.gz`).
3. Update the `REFERENCE` and `SAMPLE` variables inside `ngs_pipeline.sh` to match your reference genome and sample name.
4. Make the script executable:
   ```bash
   chmod +x ngs_pipeline.sh
   ```
5. Run the script:
   ```bash
   ./ngs_pipeline.sh
   ```

A companion `ngs_config.yaml` has also been created. Future versions of this pipeline can be modified to parse variables straight out of the YAML using tools like `yq`, making configuration much easier across multiple samples.

## Output

The pipeline automatically generates a structured output directory system branching off your `PROJECT_DIR`:
- **`qc/`**: FastQC HTML quality reports for the raw input reads.
- **`trimmed/`**: Cleaned paired-end fastq files produced by trim_galore.
- **`results/aligned/`**: Sorted BAMs, GATK duplicate metrics, and duplicate-marked BAMs.
- **`results/variants/`**: Unfiltered VCFs, hard-filtered VCFs, final `PASS` VCFs, and the final converted Excel-compatible `.csv` variant sheets.
