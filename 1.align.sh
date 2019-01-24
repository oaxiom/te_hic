
if [ ! -d run ]
then
    mkdir run
fi

# supports SE p33, p64, PE p33
for sample in samples/*
do
    sample_name=`basename $sample`

    for fq in $sample/*.fastq.gz
    do  
        name=`basename $fq | sed 's#.fastq.gz##g'`
        if [ ! -d run/$sample_name ]
        then
            mkdir run/$sample_name
        else
            rm -r run/$sample_name
            mkdir run/$sample_name
        fi
 
        qsub -N tehic.$name -v fq=$fq,out=run/$sample_name/$name align.pbs
        
    done
done

