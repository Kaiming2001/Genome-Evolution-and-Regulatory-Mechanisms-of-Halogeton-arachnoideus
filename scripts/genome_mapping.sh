#!/bin/bash
# Usage: bash run_haphic_pipeline.sh contigs.fa read1.fq.gz read2.fq.gz nchrs prefix
contigs=$1
read1=$2
read2=$3
nchrs=$4
prefix=$5
# Step 1: Index contigs
bwa index $contigs
# Step 2: Align Hi-C reads
bwa mem -5SP $contigs $read1 $read2 > ${prefix}.sam
# Step 3: Deduplicate with samblaster
samblaster < ${prefix}.sam > ${prefix}.dedup.sam
# Step 4: Convert to BAM and filter
samtools view -@ 14 -S -h -b -F 3340 -o ${prefix}-HiC.bam ${prefix}.dedup.sam
# Step 5: Filter BAM by MAPQ
/path/to/HapHiC/utils/filter_bam ${prefix}-HiC.bam 1 --nm 3 --threads 14 | samtools view - -b -@ 14 -o ${prefix}-HiC.filtered.bam
# Step 6: Clustering
/path/to/HapHiC/haphic cluster $contigs ${prefix}-HiC.filtered.bam $nchrs
# Step 7: Reassignment
/path/to/HapHiC/haphic reassign $contigs full_links.pkl mcl_inflation_1.6.clusters.txt paired_links.clm --nclusters $nchrs
# Step 8: Ordering and orientation
/path/to/HapHiC/haphic sort $contigs HT_links.pkl split_clms reassignment/final_groups/group*.txt --processes 8
# Step 9: Build pseudomolecules
/path/to/HapHiC/haphic build $contigs $contigs ${prefix}-HiC.filtered.bam ordering-orientation/group*.tour
# Step 10: juicebox
bash juicebox.sh
