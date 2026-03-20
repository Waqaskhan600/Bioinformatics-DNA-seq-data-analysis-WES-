#!/bin/bash
set -e

# ==========================================
# Basic NGS Pipeline: FastQ to Filtered VCF
# ==========================================

# Configuration variables (Ideally passed via command line or a config file)
PROJECT_DIR="/media/user/New_Volume/ngs" # Path to the project directory
REFERENCE="${PROJECT_DIR}/resources/reference/hg38.fa" # Path to the reference genome
KNOWN_SITES_DBSNP="${PROJECT_DIR}/resources/reference/Homo_sapiens_assembly38.dbsnp138.vcf.gz" # Database for BQSR
KNOWN_SITES_INDELS="${PROJECT_DIR}/resources/reference/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz" # Database for BQSR
TARGET_REGIONS="${PROJECT_DIR}/resources/reference/target_regions.bed" # BED file for WES capture kit
TARGET_INTERVALS="${PROJECT_DIR}/resources/reference/target_regions.interval_list" # Interval list for WES QC metrics
RESULTS_DIR="${PROJECT_DIR}/results" # Path to the results directory
DATA_DIR="${PROJECT_DIR}/input" # update Path to the data directory
SAMPLE="demo_1" # update Sample name
THREADS=6 # update Number of threads to use

# ==========================================
# Pre-Run Validation Checks
# ==========================================

# 0. Check Software Dependencies
MISSING_DEPS=0
REQUIRED_TOOLS=("fastqc" "trim_galore" "bwa" "samtools" "gatk" "bcftools" "awk")

echo -e "\n[INFO] Validating software dependencies..."
for tool in "${REQUIRED_TOOLS[@]}"; do
    if command -v $tool &> /dev/null; then
        echo -e "  [OK] $tool is installed."
    else
        echo -e "  [MISSING] $tool is NOT found in your PATH."
        MISSING_DEPS=1
    fi
done

if [ $MISSING_DEPS -eq 1 ]; then
    echo -e "\n[ERROR] Required bioinformatics tools are missing."
    echo "Please ensure you have activated your Conda environment before running:"
    echo "    conda activate wes_pipeline"
    echo -e "If you haven't created it yet, run: conda env create -f environment.yaml\n"
    exit 1
fi

# 1. Check if Reference Genome exists
if [ ! -f "${REFERENCE}" ]; then
    echo -e "\n[ERROR] Reference genome not found at: ${REFERENCE}"
    echo "You must download and prepare a reference genome before running the pipeline."
    echo -e "\nSuggested commands to run from your project directory (${PROJECT_DIR}):"
    echo "--------------------------------------------------------"
    echo "mkdir -p reference && cd reference"
    echo "wget https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta -O hg38.fa"
    echo "bwa index hg38.fa"
    echo "gatk CreateSequenceDictionary -R hg38.fa"
    echo "samtools faidx hg38.fa"
    echo "cd .."
    echo "--------------------------------------------------------"
    echo -e "\nOnce the reference is ready, come back and run this script again!"
    exit 1
fi

# 1.5. Check if Known Sites exist
if [ ! -f "${KNOWN_SITES_DBSNP}" ] || [ ! -f "${KNOWN_SITES_INDELS}" ]; then
    echo -e "\n[ERROR] Known sites databases for BQSR not found."
    echo "Expected:"
    echo "  - ${KNOWN_SITES_DBSNP}"
    echo "  - ${KNOWN_SITES_INDELS}"
    echo "Please place the vcf.gz files into your reference directory."
    exit 1
fi

# 1.6. Check if Target Regions BED file exists
if [ ! -f "${TARGET_REGIONS}" ]; then
    echo -e "\n[ERROR] Target regions BED file not found."
    echo "Expected: ${TARGET_REGIONS}"
    echo "WES variant calling must be restricted to the capture kit targets."
    echo "Please download your specific kit's BED file (e.g. Agilent, Twist, Illumina) and place it in resources/reference/."
    exit 1
fi

# 1.7. Check if Target Interval List exists
if [ ! -f "${TARGET_INTERVALS}" ]; then
    echo -e "\n[ERROR] Target interval list not found."
    echo "Expected: ${TARGET_INTERVALS}"
    echo "CollectHsMetrics requires a Picard-style interval_list file, not just a BED file."
    echo -e "\nPlease convert your BED file using the following command from your project directory:"
    echo "--------------------------------------------------------"
    echo "gatk BedToIntervalList \\"
    echo "    -I ${TARGET_REGIONS} \\"
    echo "    -O ${TARGET_INTERVALS} \\"
    echo "    -SD resources/reference/hg38.dict"
    echo "--------------------------------------------------------"
    echo -e "\nOnce created, come back and run this script again!"
    exit 1
fi

# 2. Check if Input FastQ files exist
if [ ! -f "${DATA_DIR}/${SAMPLE}_R1.fastq.gz" ] || [ ! -f "${DATA_DIR}/${SAMPLE}_R2.fastq.gz" ]; then
    echo -e "\n[ERROR] Input data not found for sample: ${SAMPLE}"
    echo "Please ensure your raw paired-end reads are placed in the data directory: ${DATA_DIR}"
    echo "Expected files:"
    echo "  - ${SAMPLE}_R1.fastq.gz"
    echo "  - ${SAMPLE}_R2.fastq.gz"
    echo -e "\nIf the data folder doesn't exist yet, you can create it with:"
    echo "  mkdir -p ${DATA_DIR}"
    echo -e "\nOnce your files are in the correct location, come back and run this script again!"
    exit 1
fi

# Ensure output directories exist
mkdir -p ${RESULTS_DIR}/aligned
mkdir -p ${RESULTS_DIR}/variants
mkdir -p ${RESULTS_DIR}/qc
mkdir -p ${RESULTS_DIR}/trimmed

echo -e "\n[SUCCESS] All checks passed! Starting basic NGS pipeline for ${SAMPLE}...\n"

# 0. Quality Control & Trimming
echo "Step 0.1: Running FastQC on raw reads..."
fastqc -t ${THREADS} -o ${RESULTS_DIR}/qc \
    ${DATA_DIR}/${SAMPLE}_R1.fastq.gz ${DATA_DIR}/${SAMPLE}_R2.fastq.gz

echo "Step 0.2: Trimming adapters and low-quality bases..."
trim_galore --paired --fastqc --cores ${THREADS} \
    -o ${RESULTS_DIR}/trimmed \
    ${DATA_DIR}/${SAMPLE}_R1.fastq.gz ${DATA_DIR}/${SAMPLE}_R2.fastq.gz

# trim_galore default output naming for paired reads
TRIMMED_R1="${RESULTS_DIR}/trimmed/${SAMPLE}_R1_val_1.fq.gz"
TRIMMED_R2="${RESULTS_DIR}/trimmed/${SAMPLE}_R2_val_2.fq.gz"

# 1. Alignment with BWA-MEM
echo "Step 1: Aligning reads..."
bwa mem -t ${THREADS} -M -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:ILLUMINA\tLB:lib1" \
    ${REFERENCE} \
    ${TRIMMED_R1} ${TRIMMED_R2} | \
    samtools sort -@ ${THREADS} -o ${RESULTS_DIR}/aligned/${SAMPLE}_sorted.bam

samtools index ${RESULTS_DIR}/aligned/${SAMPLE}_sorted.bam

# 2. Mark Duplicates with GATK
echo "Step 2: Marking PCR and optical duplicates..."
gatk MarkDuplicates \
    -I ${RESULTS_DIR}/aligned/${SAMPLE}_sorted.bam \
    -O ${RESULTS_DIR}/aligned/${SAMPLE}_marked_dup.bam \
    -M ${RESULTS_DIR}/aligned/${SAMPLE}_dup_metrics.txt

samtools index ${RESULTS_DIR}/aligned/${SAMPLE}_marked_dup.bam

# 3. Base Quality Score Recalibration (BQSR)
echo "Step 3: Generating BQSR recalibration table..."
gatk BaseRecalibrator \
    -R ${REFERENCE} \
    -I ${RESULTS_DIR}/aligned/${SAMPLE}_marked_dup.bam \
    --known-sites ${KNOWN_SITES_DBSNP} \
    --known-sites ${KNOWN_SITES_INDELS} \
    -O ${RESULTS_DIR}/aligned/${SAMPLE}_recal_data.table

echo "Step 3.1: Applying BQSR..."
gatk ApplyBQSR \
    -R ${REFERENCE} \
    -I ${RESULTS_DIR}/aligned/${SAMPLE}_marked_dup.bam \
    --bqsr-recal-file ${RESULTS_DIR}/aligned/${SAMPLE}_recal_data.table \
    -O ${RESULTS_DIR}/aligned/${SAMPLE}_bqsr.bam

samtools index ${RESULTS_DIR}/aligned/${SAMPLE}_bqsr.bam

# 4. WES Specific Quality Metrics
echo "Step 4: Calculating On-Target capture metrics..."
gatk CollectHsMetrics \
    -I ${RESULTS_DIR}/aligned/${SAMPLE}_bqsr.bam \
    -O ${RESULTS_DIR}/qc/${SAMPLE}_hs_metrics.txt \
    -R ${REFERENCE} \
    -BAIT_INTERVALS ${TARGET_INTERVALS} \
    -TARGET_INTERVALS ${TARGET_INTERVALS}
 
echo "Step 4.1: Calculating detailed coverage statistics..."
gatk DepthOfCoverage \
    -R ${REFERENCE} \
    -I ${RESULTS_DIR}/aligned/${SAMPLE}_bqsr.bam \
    -O ${RESULTS_DIR}/qc/${SAMPLE}_coverage \
    -L ${TARGET_REGIONS} \
    --summaryCoverageThreshold 10 \
    --summaryCoverageThreshold 20 \
    --summaryCoverageThreshold 30 \
    --summaryCoverageThreshold 50 \
    --summaryCoverageThreshold 100

echo "Step 4.2: Collecting general alignment summary metrics..."
gatk CollectAlignmentSummaryMetrics \
    -R ${REFERENCE} \
    -I ${RESULTS_DIR}/aligned/${SAMPLE}_bqsr.bam \
    -O ${RESULTS_DIR}/qc/${SAMPLE}_alignment_summary.txt
 
echo "Step 4.3: Collecting paired-end insert size metrics..."
gatk CollectInsertSizeMetrics \
    -I ${RESULTS_DIR}/aligned/${SAMPLE}_bqsr.bam \
    -O ${RESULTS_DIR}/qc/${SAMPLE}_insert_size_metrics.txt \
    -H ${RESULTS_DIR}/qc/${SAMPLE}_insert_size_histogram.pdf

# 5. Variant Calling with GATK HaplotypeCaller
echo "Step 5: Variant Calling..."
gatk HaplotypeCaller \
    -R ${REFERENCE} \
    -I ${RESULTS_DIR}/aligned/${SAMPLE}_bqsr.bam \
    -L ${TARGET_REGIONS} \
    -O ${RESULTS_DIR}/variants/${SAMPLE}_raw_variants.vcf.gz

# 6. Apply variant filters optimized for high-coverage WES data
echo "Step 6: Filtering variants..."
gatk VariantFiltration \
    -R ${REFERENCE} \
    -V ${RESULTS_DIR}/variants/${SAMPLE}_raw_variants.vcf.gz \
    --filter-expression "QD < 2.0" --filter-name "LowQD" \
    --filter-expression "FS > 60.0" --filter-name "HighFS" \
    --filter-expression "MQ < 40.0" --filter-name "LowMQ" \
    --filter-expression "SOR > 3.0" --filter-name "HighSOR" \
    --filter-expression "MQRankSum < -12.5" --filter-name "LowMQRankSum" \
    --filter-expression "ReadPosRankSum < -8.0" --filter-name "LowReadPosRankSum" \
    --filter-expression "DP < 20" --filter-name "LowDepth" \
    --filter-expression "DP > 500" --filter-name "HighDepth" \
    -O ${RESULTS_DIR}/variants/${SAMPLE}_filtered.vcf.gz
 
# 7. Extract high-quality variants
echo "Step 7: Extracting PASS variants..."
bcftools view -f PASS -Oz \
    -o ${RESULTS_DIR}/variants/${SAMPLE}_pass.vcf.gz \
    ${RESULTS_DIR}/variants/${SAMPLE}_filtered.vcf.gz
 
bcftools index -t ${RESULTS_DIR}/variants/${SAMPLE}_pass.vcf.gz
# 8. Convert VCF to CSV format for Excel
echo "Step 8: Converting PASS variants to CSV format..."
gatk VariantsToTable \
    -V ${RESULTS_DIR}/variants/${SAMPLE}_pass.vcf.gz \
    -F CHROM -F POS -F ID -F REF -F ALT -F QUAL -F FILTER \
    -F DP -F QD -F FS -F MQ -F SOR -F MQRankSum -F ReadPosRankSum \
    -O ${RESULTS_DIR}/variants/${SAMPLE}_pass.tsv

# Convert the TSV output to true comma-separated CSV
awk 'BEGIN { FS="\t"; OFS="," } {$1=$1; print}' ${RESULTS_DIR}/variants/${SAMPLE}_pass.tsv > ${RESULTS_DIR}/variants/${SAMPLE}_pass.csv

echo -e "\n[SUCCESS] Pipeline complete! High-quality variants for ${SAMPLE} are saved at:"
echo " - VCF (for tools): ${RESULTS_DIR}/variants/${SAMPLE}_pass.vcf.gz"
echo " - CSV (for Excel): ${RESULTS_DIR}/variants/${SAMPLE}_pass.csv"
