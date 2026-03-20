# Basic Whole Exome Sequencing (WES) Pipeline

This repository contains a simple, linear Next-Generation Sequencing (NGS) pipeline using bash for paired-end FastQ data processing, optimized for Whole Exome Sequencing (WES).

## Pipeline Steps

The pipeline incorporates the following standard industry practices:
0. **Quality Control & Trimming:** `fastqc` assesses read quality and `trim_galore` removes sequencing adapters and low-quality sequences (defaulting to `--quality 20` and `--length 50`) to minimize noise. For advanced trimming options, consult the [Trim Galore User Guide](https://github.com/FelixKrueger/TrimGalore/blob/master/Docs/Trim_Galore_User_Guide.md).
1. **Alignment:** Trimmed reads are aligned to the reference genome (`hg38.fa`) using `bwa mem`, then sorted via `samtools sort`.
2. **Duplicate Marking:** Eliminates PCR or optical duplicates using `gatk MarkDuplicates` to prevent biased variant calling.
3. **Base Quality Score Recalibration (BQSR):** Detects and corrects systematic errors made by the sequencing machine in estimating base quality scores using `gatk BaseRecalibrator` and `ApplyBQSR`.
4. **WES Quality Metrics:** Calculates general alignment statistics, insert sizes, and on-target capture efficiency using `gatk CollectAlignmentSummaryMetrics`, `CollectInsertSizeMetrics`, `CollectHsMetrics` and `DepthOfCoverage`.
5. **Variant Calling:** Scans the recalibrated alignment to identify SNPs and Indels using `gatk HaplotypeCaller`.
6. **Variant Filtration:** Applies rigorous, fixed thresholds adapted specifically for high-coverage WES data using `gatk VariantFiltration`.
7. **Quality Extraction:** Excludes heavily filtered calls, leaving only high-confidence (`PASS`) variants using `bcftools view`.
8. **CSV Conversion:** Converts the final `PASS` variants VCF into an Excel-compatible Comma-Separated Values (`.csv`) file.
9. **Variant Annotation:** Structurally and functionally annotates the resulting variants to understand their biological impact.

## Installation & Setup

### 0. Clone the Repository
First, clone this pipeline repository to your local machine and navigate into it:
```bash
git clone https://github.com/Waqaskhan600/Bioinformatics-DNA-seq-data-analysis-WES-.git
cd Bioinformatics-DNA-seq-data-analysis-WES-
```

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

### What if I don't have a Target BED file?

While target files are **highly recommended** for WES to calculate capture efficiency and restrict variant calling to valid exons, they are not strictly mandatory to generate a VCF. If you absolutely do not have the BED file from the sequencing center:
1. You **must** comment out or remove `gatk CollectHsMetrics` and `gatk DepthOfCoverage` in the script.
2. You **must** remove the `-L ${TARGET_REGIONS}` flag from `gatk HaplotypeCaller` (Step 5). 

*Note: Without the `-L` flag, HaplotypeCaller will blindly search the entire genome for variants instead of just the exome, which will take significantly longer. Furthermore, you will lose all metrics on whether the physical capture kit actually worked and whether you achieved enough depth to trust your variants!*

## Usage

*Note: The pipeline includes built-in **Pre-Run Validation Checks** that will ensure you have correctly downloaded the reference genome, installed all software dependencies, and placed your input data before executing any bioinformatics tools. If anything is missing (for example, if you forgot to run `conda activate wes_pipeline`), the script will catch the missing dependencies, halt immediately, and provide the exact commands you need to fix the issue!*

1. **Configure the pipeline:** Open `ngs_pipeline.sh` and define the `PROJECT_DIR` variable to point to your project's absolute path. If you are running the script directly from within your project directory, you can simply set it to `PROJECT_DIR="$PWD"`. All other directories (data, outputs, etc.) will branch off from this folder.
2. Place your raw input reads in the `data/` directory inside your project folder (e.g., `sample_1_R1.fastq.gz`, `sample_1_R2.fastq.gz`).
   * **Tip on Storage:** While this pipeline dynamically accepts uncompressed `.fastq` files automatically, they are incredibly large. It is considered an industry best practice to compress them yourself to save immense disk space (70-80% reduction) by running `gzip input/*.fastq` before starting.
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
- **`aligned/`**: Sorted BAMs, GATK duplicate metrics, and duplicate-marked BAMs.
- **`variants/`**: Unfiltered VCFs, hard-filtered VCFs, final `PASS` VCFs, and Excel-compatible `.csv` variant sheets.
- **`logs/`**: Detailed terminal console outputs for every single step of the pipeline. If a step fails or produces unexpected results, check these log files for exact tool error messages.

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

When evaluating the `results/qc/` reports (`_hs_metrics.txt`, `_alignment_summary.txt`, and `_insert_size_metrics.txt`), keep an eye on these critical WES-specific metrics:

1. **PCT_PF_READS_ALIGNED (Alignment Summary):** Percentage of reads aligned to the reference.
   * Standard: >95% for reliable human read data.
2. **PF_MISMATCH_RATE (Alignment Summary):** Rate of mismatches to the reference.
   * Standard: <5% (high mismatch rates imply DNA degradation or species contamination).
3. **INSERT_SIZE (Insert Size Metrics):** Distance between paired reads.
   * Expected: WES target kits typically fragment DNA smaller than WGS to tightly fit exons. Expect a clean bell curve around **150-250bp**.
4. **ON_TARGET_BASES (HsMetrics):** Percentage of sequenced bases that fall within your target capture regions.
   * Excellent: >80%
   * Good: 70-80%
   * Poor: <70%
5. **MEAN_TARGET_COVERAGE:** Average coverage depth across target regions.
   * Clinical applications: >100x
   * Research applications: >50x
6. **PCT_TARGET_BASES_10X/20X/30X:** Percentage of targets successfully covered at various depths.
   * Clinical quality: >95% at 20x coverage
7. **FOLD_ENRICHMENT:** How much more coverage the targeted regions received versus the genome average.
   * Good enrichment: >40-fold
   * Poor enrichment: <20-fold
8. **AT_DROPOUT and GC_DROPOUT:** Coverage uniformity across different GC content.
   * Good uniformity: <10% dropout
   * Problematic: >20% dropout

## Best Practices for WES Analysis

### Quality Control Standards
**Coverage Requirements (Review your `qc/*_hs_metrics.txt`):**
* **Mean target coverage:** >100x for clinical applications, >50x for research
* **Coverage uniformity:** >95% of targets at 20x depth
* **On-target rate:** >80% for optimal results

**Performance Benchmarks:**
* **Variant Calling Sensitivity:** >95% for variants >20x coverage (via `gatk HaplotypeCaller`)
* **Specificity:** >99% after applying the hard filters generated by `gatk VariantFiltration`

## Common Pitfalls and Solutions

### Technical Issues

1. **Low On-Target Rate (<70%)**
   * **Causes:** Wrong interval list, poor physical capture, sample degradation.
   * **Solutions:** Verify your `TARGET_INTERVALS` and `TARGET_REGIONS` variables point precisely to the correct manufacturer capture kit.
   * **Prevention:** Always validate protocols using exact manufacturer BED files.

2. **Uneven Coverage Across Targets**
   * **Causes:** GC bias, capture probe efficiency, PCR artifacts.
   * **Solutions:** Review the AT/GC dropout metrics in the QC folder `hs_metrics`.
   * **Monitoring:** Track coverage uniformity metrics output by `gatk DepthOfCoverage`.

3. **High False Positive Rate**
   * **Causes:** Insufficient filtering, capture artifacts, sequence biases.
   * **Solutions:** Adjust the `gatk VariantFiltration` parameters specifically for your sample's coverage depth. The defaults we configured (`QD < 2.0`, `DP < 20`) are a baseline for high-coverage WES.
   * **Validation:** Confirm highly suspicious output "PASS" variants with Sanger sequencing.

4. **Poor Raw Read Quality or Technical Bias at Read Extremes**
   * **Indicator:** The raw FastQC report shows erratic "Per Base Sequence Content" at the 5' end, or adapters are failing to be removed.
   * **Solutions:** The pipeline defaults to standard, safe trimming (`--quality 20` and `--length 50`). If you spot persistent issues in FastQC, you can add advanced flags to the `trim_galore` command in `Step 0.2` of the script:
     * **`--illumina` or `--nextera`**: By default, `trim_galore` is incredibly smart and will auto-detect which adapter your sequencing center used. However, if you know for an absolute fact you used a Nextera Exome kit, you can force it with `--nextera`.
     * **`--clip_R1 5` and `--clip_R2 5`**: If the "Per Base Sequence Content" graph is chaotic for the very first 5 to 10 base pairs (very common in Illumina WES), this flag acts like a guillotine and mindlessly chops off the first 5 bases of every single read, regardless of quality, to remove that technical bias.

### Analytical Issues

1. **Missing Expected Germline Variants**
   * **Causes:** Low coverage in specific exons, or the variant is located in a non-captured region outside the baits.
   * **Solutions:** Manually review the `DepthOfCoverage` file for the specific gene coordinate to ensure it was properly sequenced.
   * **Prevention:** Always review your capture kit's exact probe content for genes of interest beforehand.

## Contact & Collaboration

Thank you for using this pipeline!

If you found this useful and informative, or if you are looking for collaboration and research assistance, please feel free to reach out to me on LinkedIn:
[Waqas Khan - LinkedIn](https://www.linkedin.com/in/waqas-khan-3b937b184/)