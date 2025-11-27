#!/bin/bash
# Usage: bash genome_assembly.sh input.bam
inbam=$1
prefix=$(basename "$inbam" .bam)
# Step 1: BAM -> FASTA by samtools v1.10
samtools view "$inbam" | awk '{OFS="\t"; print ">"$1"\n"$10}' - > "${prefix}.fa"
# Step 2: Run hifiasm assembly with Hifiasm v0.20.0-r639
hifiasm -o "${prefix}.asm" -t 20 "${prefix}.fa"
# Step 3: Extract primary contigs
awk '/^S/{print ">"$2;print $3}' "${prefix}.asm.bp.p_ctg.gfa" > "${prefix}.p_ctg.fa"
