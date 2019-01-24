

for sample in run/*
do
    sample_name=`basename $sample`

    for bam in $sample/*_1.bam
    do  
        p1=`basename $bam`
        p2=`echo $p1 | sed 's#_1.bam#_2.bam#g'`
        out=`echo $p1 | sed 's#_1.bam##g'`
        collect_valid_pairs.py $sample/$p1 $sample/$p2 $out.bedpe #>$out.out
    done
done



