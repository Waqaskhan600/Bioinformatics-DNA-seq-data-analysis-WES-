# Basic Whole Exome Sequencing (WES) Pipeline

This repository contains a simple, linear Next-Generation Sequencing (NGS) pipeline using bash for paired-end FastQ data processing, optimized for Whole Exome Sequencing (WES).

## Pipeline Steps

The pipeline incorporates the following standard industry practices:
0. **Quality Control & Trimming:** `fastqc` assesses read quality and `trim_galore` removes sequencing adapters and low-quality sequences to minimize noise.
1. **Alignment:** Trimmed reads are aligned to the reference genome (`hg38.fa`) using `bwa mem`, then sorted via `samtools sort`.
2. **Duplicate Marking:** Eliminates PCR or optical duplicates using `gatk MarkDuplicates` to prevent biased variant calling.
3. **Base Quality Score Recalibration (BQSR):** Detects and corrects systematic errors made by the sequencing machine in estimating base quality scores using `gatk BaseRecalibrator` and `ApplyBQSR`.
4. **WES Quality Metrics:** Calculates on-target capture efficiency and multiplexed depth coverage statistics using `gatk CollectHsMetrics` and `gatk DepthOfCoverage`.
5. **Variant Calling:** Scans the recalibrated alignment to identify SNPs and Indels using `gatk HaplotypeCaller`.
6. **Variant Filtration:** Applies rigorous, fixed thresholds adapted specifically for high-coverage WES data using `gatk VariantFiltration`.
7. **Quality Extraction:** Excludes heavily filtered calls, leaving only high-confidence (`PASS`) variants using `bcftools view`.
8. **CSV Conversion:** Converts the final `PASS` variants VCF into an Excel-compatible Comma-Separated Values (`.csv`) file.
9. **Variant Annotation:** Structurally and functionally annotates the resulting variants to understand their biological impact.

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

## WES Target Capture Kits

Whole Exome Sequencing explicitly targets specific regions of the genome using probes or "baits". **Critical Point:** Always use the exact target region file (`.bed` or `.interval_list`) that matches your capture kit. Using the wrong target file will lead to incorrect quality metrics, missed variants, and suboptimal variant calling. 

You must place your target regions BED file at `resources/reference/target_regions.bed` (or configure `TARGET_REGIONS` in the shell script). Common panels include:
* **Agilent SureSelect** – Uses RNA baits, excellent uniformity
* **Illumina Nextera** – DNA baits, good for degraded samples
* **Roche SeqCap** – DNA baits, customizable designs
* **Twist Bioscience** – Synthetic DNA baits, latest technology

**GATK Interval List Format:**
Critical quality control tools like `CollectHsMetrics` require the BED file to be converted into a GATK-specific `.interval_list` format with a standard SAM/BAM header. If you only have the `.bed` file, you can easily convert it using:
```bash
# Create interval list for GATK
gatk BedToIntervalList \
    -I resources/reference/target_regions.bed \
    -O resources/reference/target_regions.interval_list \
    -SD resources/reference/hg38.dict
```

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

## Variant Annotation

Once you have your final VCF, you can heavily annotate your variants to predict their clinical or biological consequences through various external tools such as:
* [wANNOVAR](https://wannovar.wglab.org/) 
* [Ensembl VEP](https://www.ensembl.org/Tools/VEP)
* SnpEff
* GATK Funcotator

### Why Use Multiple Annotation Tools?
No single tool provides complete information. Each has unique strengths:

* **GATK Funcotator** – GATK’s native annotator with a clinical database focus
* **Ensembl VEP** – Comprehensive consequence prediction with extensive options
* **SnpEff** – Fast annotation with powerful filtering capabilities
* **ANNOVAR** – Excellent for population frequency analysis

By using multiple tools, you gain different perspectives on each variant and build confidence in your interpretations.

## Understanding WES Quality Metrics

When evaluating the `results/qc/*_hs_metrics.txt` and `coverage` files, keep an eye on these critical WES-specific metrics:

1. **ON_TARGET_BASES:** Percentage of sequenced bases that fall within your target capture regions.
   * Excellent: >80%
   * Good: 70-80%
   * Poor: <70%
2. **MEAN_TARGET_COVERAGE:** Average coverage depth across target regions.
   * Clinical applications: >100x
   * Research applications: >50x
3. **PCT_TARGET_BASES_10X/20X/30X:** Percentage of targets successfully covered at various depths.
   * Clinical quality: >95% at 20x coverage
4. **FOLD_ENRICHMENT:** How much more coverage the targeted regions received versus the genome average.
   * Good enrichment: >40-fold
   * Poor enrichment: <20-fold
5. **AT_DROPOUT and GC_DROPOUT:** Coverage uniformity across different GC content.
   * Good uniformity: <10% dropout
   * Problematic: >20% dropout

## Best Practices for WES Analysis

### Quality Control Standards
**Coverage Requirements:**
* **Mean target coverage:** >100x for clinical applications, >50x for research
* **Coverage uniformity:** >95% of targets at 20x depth
* **On-target rate:** >80% for optimal results

**Performance Benchmarks:**
* **Variant calling sensitivity:** >95% for variants >20x coverage
* **Specificity:** >99% after proper filtering
* **Reproducibility:** >98% concordance between technical replicates

## Common Pitfalls and Solutions

### Technical Issues

1. **Low On-Target Rate (<70%)**
   * **Causes:** Wrong target file, poor capture, sample degradation
   * **Solutions:** Verify target file, check capture kit lot, assess DNA quality
   * **Prevention:** Use exact manufacturer target files, validate protocols

2. **Uneven Coverage Across Targets**
   * **Causes:** GC bias, capture probe efficiency, PCR artifacts
   * **Solutions:** GC bias correction, optimize PCR conditions
   * **Monitoring:** Track coverage uniformity metrics

3. **High False Positive Rate**
   * **Causes:** Insufficient filtering, capture artifacts, systematic biases
   * **Solutions:** Apply WES-specific filters, use population databases
   * **Validation:** Confirm suspicious variants with Sanger sequencing

### Analytical Issues

1. **Missing Expected Variants**
   * **Causes:** Low coverage in specific exons, variant in non-captured regions
   * **Solutions:** Manual review of coverage, consider WGS for negative cases
   * **Prevention:** Review capture kit content for genes of interest

2. **CNV False Positives**
   * **Causes:** Capture bias, reference sample selection, systematic artifacts
   * **Solutions:** Use WES-specific CNV tools, validate with orthogonal methods
   * **Prevention:** Include adequate normal samples in reference

3. **Batch Effects**
   * **Causes:** Different capture lots, library prep variations, sequencing platforms
   * **Solutions:** Batch correction methods, standardized protocols
   * **Prevention:** Process samples consistently, track batch information

4. **Somatic Mutation Artifacts**
   * **Causes:** DNA degradation, FFPE artifacts, low tumor purity
   * **Solutions:** Optimize filtering parameters, use orientation bias models
   * **Prevention:** Use high-quality samples, validate tumor content
