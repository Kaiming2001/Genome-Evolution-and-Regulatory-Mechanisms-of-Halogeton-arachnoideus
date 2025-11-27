#!/bin/bash

# =========================================================
# RNA-seq analysis pipeline:
#   1. Quality filtering and adapter trimming (Fastp)
#   2. Read alignment to reference genome (HISAT2)
#   3. Transcript quantification (StringTie)
# =========================================================

# --- 0. Define directories  ---
WORKDIR="./raw_data"             # raw fastq files directory: The raw fastq data has been uploaded to NGDC GSA (CRA031549) under BioProject PRJCA041342.
CLEAN_DIR="./clean_data"         # output from fastp
ALIGN_DIR="./hisat2_alignments"  # output BAMs
ASSEMBLY_DIR="./stringtie"       # transcript quantification results
REF_DIR="./reference"            # reference genome and annotation: The relevant assembly results have been uploaded to FigShare (DOI: 10.6084/m9.figshare.30353017).
INDEX_PREFIX="${REF_DIR}/genome" # HISAT2 index prefix
REF_GFF="${REF_DIR}/genome.gff3"

# --- 1. Load software (or ensure PATH is configured) ---
FASTP=$(which fastp)
HISAT2=$(which hisat2)
SAMTOOLS=$(which samtools)
STRINGTIE=$(which stringtie)

# --- 2. Create output directories ---
mkdir -p ${CLEAN_DIR} ${ALIGN_DIR} ${ASSEMBLY_DIR}

# =========================================================
# STEP 1: Quality control and adapter trimming (Fastp)
# =========================================================
echo "Step 1: Running Fastp..."

for SAMPLE in $(ls ${WORKDIR}/*_R1.fastq.gz | sed 's/_R1.fastq.gz//' | xargs -n1 basename); do
    echo "  Processing ${SAMPLE}..."
    mkdir -p ${CLEAN_DIR}/${SAMPLE}

    ${FASTP} \
        -i ${WORKDIR}/${SAMPLE}_R1.fastq.gz \
        -I ${WORKDIR}/${SAMPLE}_R2.fastq.gz \
        -o ${CLEAN_DIR}/${SAMPLE}/${SAMPLE}_clean_R1.fq.gz \
        -O ${CLEAN_DIR}/${SAMPLE}/${SAMPLE}_clean_R2.fq.gz \
done

# =========================================================
# STEP 2: Alignment with HISAT2
# =========================================================
echo "Step 2: Running HISAT2 alignment..."

# Build index if not exist
if [ ! -e "${INDEX_PREFIX}.1.ht2" ]; then
    echo "  Building HISAT2 index..."
    ${HISAT2}-build -p 8 ${REF_DIR}/genome.fa ${INDEX_PREFIX}
fi

for SAMPLE in $(ls ${CLEAN_DIR}); do
    echo "  Aligning ${SAMPLE}..."
    mkdir -p ${ALIGN_DIR}/${SAMPLE}
    
    ${HISAT2} -p 8 \
        -x ${INDEX_PREFIX} \
        -1 ${CLEAN_DIR}/${SAMPLE}/${SAMPLE}_clean_R1.fq.gz \
        -2 ${CLEAN_DIR}/${SAMPLE}/${SAMPLE}_clean_R2.fq.gz \
        -S ${ALIGN_DIR}/${SAMPLE}/${SAMPLE}.sam \
        --summary-file ${ALIGN_DIR}/${SAMPLE}/${SAMPLE}.hisat2.log
    
    ${SAMTOOLS} sort -@ 8 -o ${ALIGN_DIR}/${SAMPLE}/${SAMPLE}.sorted.bam ${ALIGN_DIR}/${SAMPLE}/${SAMPLE}.sam
    ${SAMTOOLS} index ${ALIGN_DIR}/${SAMPLE}/${SAMPLE}.sorted.bam
    rm ${ALIGN_DIR}/${SAMPLE}/${SAMPLE}.sam
done

# =========================================================
# STEP 3: Quantify transcript abundance with StringTie (TPM)
# =========================================================
echo "Step 3: Running StringTie quantification..."

for SAMPLE in $(ls ${ALIGN_DIR}); do
    echo "  Quantifying ${SAMPLE}..."
    mkdir -p ${ASSEMBLY_DIR}/${SAMPLE}

    ${STRINGTIE} \
        -p 8 \
        -e -G ${REF_GFF} \
        -o ${ASSEMBLY_DIR}/${SAMPLE}/${SAMPLE}.gtf \
        -A ${ASSEMBLY_DIR}/${SAMPLE}/${SAMPLE}_gene_abundance.tsv \
        ${ALIGN_DIR}/${SAMPLE}/${SAMPLE}.sorted.bam
done

echo "RNA-seq pipeline completed successfully!"
