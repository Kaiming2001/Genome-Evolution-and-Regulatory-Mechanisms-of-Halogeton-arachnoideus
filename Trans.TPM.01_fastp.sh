#!/bin/bash
workdir=/data/users/xukm2024/hara/TE/new
outdir=/data/users/xukm2024/hara/TE/new_fastp
fastp=/home/xukm2024/software/anaconda3/envs/fastp/bin/fastp

for dir in $workdir/*
do
	if [ -d $dir ]; then
	name=${dir##*/}
	r1="${name}_combined_R1.fastq.gz"
	r2="${name}_combined_R2.fastq.gz"
	fi
echo -e "mkdir -p $outdir/${name};cd $outdir/${name};$fastp -i $workdir/${name}/$r1 -o $outdir/${name}/${name}.clean_1.fq.gz -I $workdir/${name}/$r2 -O $outdir/${name}/${name}.clean_2.fq.gz " >> fastp.job.sh  
done
