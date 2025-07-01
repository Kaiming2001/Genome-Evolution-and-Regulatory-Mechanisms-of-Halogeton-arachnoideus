hisat2=/home/xukm2024/software/anaconda3/envs/hisat2/bin/hisat2
Indir=/home/xukm2024/workspace/TE/fastq/24h
Outdir=/home/xukm2024/workspace/TE/hisat2
ref_dir=/home/xukm2024/workspace/TE/index
samtools=/home/xukm2024/software/anaconda3/envs/samtools/bin/samtools

#1. build ref index

#cd $ref;$hisat2-build -p 8 $ref/genome.fasta genome

#2.get TPM/FPKM 

for dir in $Indir/*
do
if [ -d $dir ]; then
name=${dir##*/}
echo -e "mkdir -p $Outdir/$name; cd $Outdir/$name;$hisat2 -p 5 -x $ref_dir/hara_chr -1 $Indir/$name/$name.clean_1.fq.gz -2 $Indir/$name/$name.clean_2.fq.gz -S $Outdir/$name/$name.sam --summary-file $Outdir/$name/$name.hisat.out;$samtools sort $Outdir/$name/$name.sam -O BAM -T $name -l 3 -o $Outdir/$name/$name.sort.bam  && rm $Outdir/$name/$name.sam">> hisat24h.job.sh
fi
done
