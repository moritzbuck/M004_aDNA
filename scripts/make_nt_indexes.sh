fold=/home/moritz/scratch/dbs/nt

for f in `ls $fold | grep "[013-9]" `
do
    id=${f##nt.fna.}
#    id=`printf "%03d" $id`
    if [ ! -f  $fold/bowtie2/nt.${id}.1.bt2l ]
       then
	   bowtie2-build  --threads 20  nt.fna.$id  bowtie2/nt.$id > /dev/null
    else
	echo $id done
    fi
done


for f in `ls $fold | grep "[013-9]" `
do
    id=${f##nt.fna.}
    if [ ! -f  $fold/bowtie2/nt.${id}.1.bt2l ]
    then
bowtie2 -x $fold/bowtie2/nt.$id -U $samp_path/${samp}.fastq.gz -D 20 -R 3 -N 1 -L 20 -i S,1,0.50 --end-to-end --no-unal  --threads 20  -S ${samp}_1019.$id.sam -t 
    else
	echo $id done
    fi
done
 
	 

