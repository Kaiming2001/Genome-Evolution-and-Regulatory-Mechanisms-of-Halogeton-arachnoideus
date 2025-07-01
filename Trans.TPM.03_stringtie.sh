samtools=/home/xukm2024/software/anaconda3/envs/samtools/bin/samtools
stringtie=/home/xukm2024/software/anaconda3/envs/stringtie/bin/stringtie
Indir=/home/xukm2024/workspace/TE
gff=/home/xukm2024/workspace/Hara.GeneModels.gff3


    for dir in $Indir/hisat2/*
    do
	if [ -d $dir ]; then
	name=${dir##*/}
	fi
    echo -e "mkdir -p $Indir/Stringtie/$name;cd $Indir/Stringtie/$name;$stringtie -p 5 -e -G $gff -o $Indir/Stringtie/$name/$name.gtf $Indir/hisat2/$name/$name.sort.bam" >> stringtie.job.sh
    done
