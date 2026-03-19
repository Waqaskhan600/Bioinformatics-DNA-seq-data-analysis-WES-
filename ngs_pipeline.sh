#!/bin/bash
set -e

# ==========================================
# Basic NGS Pipeline: FastQ to Filtered VCF
# ==========================================

# Configuration variables (Ideally passed via command line or a config file)
PROJECT_DIR="/media/user/New_Volume/ngs" # update Path to the project directory
REFERENCE="${PROJECT_DIR}/resources/reference/hg38.fa" # update Path to the reference genome
RESULTS_DIR="${PROJECT_DIR}/results" # update Path to the results directory
DATA_DIR="${PROJECT_DIR}/input" # update Path to the data directory
SAMPLE="demo_1" # update Sample name
THREADS=6 # update Number of threads to use according your system

# ==========================================
# Pre-Run Validation Checks
# ==========================================

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

# 3. Variant Calling with GATK HaplotypeCaller
echo "Step 3: Variant Calling..."
gatk HaplotypeCaller \
    -R ${REFERENCE} \
    -I ${RESULTS_DIR}/aligned/${SAMPLE}_marked_dup.bam \
    -O ${RESULTS_DIR}/variants/${SAMPLE}_raw_variants.vcf.gz

# 4. Apply variant filters optimized for high-coverage WES data
echo "Step 4: Filtering variants..."
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
 
# 5. Extract high-quality variants
echo "Step 5: Extracting PASS variants..."
bcftools view -f PASS -Oz \
    -o ${RESULTS_DIR}/variants/${SAMPLE}_pass.vcf.gz \
    ${RESULTS_DIR}/variants/${SAMPLE}_filtered.vcf.gz
 
bcftools index -t ${RESULTS_DIR}/variants/${SAMPLE}_pass.vcf.gz
# 6. Convert VCF to CSV format for Excel
echo "Step 6: Converting PASS variants to CSV format..."
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
