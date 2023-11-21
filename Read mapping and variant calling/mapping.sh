#!/bin/bash

ref=DM_v6.1.fa
threads=12

cd /media/scratchpad_01/li322/work/02_137_atlas_potato_reseq
###### mapping

bwa index $ref
gatk CreateSequenceDictionary -R $ref -O DM_v6.1.dict
samtools faidx $ref

rm -rf mapping; mkdir mapping
#!/bin/bash

###### bwa mem -R "@RG\tID:VTN62-33-3\tSM:VTN62-33-3\tPL:Illumina"
ls reseq_data_clean | while read i
do
	name=`ls reseq_data_clean/${i}/*_R1.clean.fastq.gz | awk -F '/' '{print $NF}'|awk -F '_' '{print $1}'`
	echo """#!/bin/bash
	cd /media/scratchpad_01/li322/work/02_137_atlas_potato_reseq
	
	fastp -i li/${sample}-${rep}_*R1_001_AH5N5CDSX5.filt.fastq.gz -I li/${sample}-${rep}_*R2_001_AH5N5CDSX5.filt.fastq.gz -o li_rnaseq_fastp/${sample}-${rep}_1.clean.fq.gz -O li_rnaseq_fastp/${sample}-${rep}_2.clean.fq.gz -w 16
    bwa mem -R \"@RG\\\tID:${name}\\\tSM:${name}\\\tPL:Illumina\" -t $threads $ref <(zcat reseq_data_clean/${i}/*_R1.clean*gz) <(zcat reseq_data_clean/${i}/*_R2.clean*gz) | samtools sort -O bam -@ $threads -m 2G > mapping/${name}.sort.bam 
    gatk --java-options \"-Xmx10G -XX:ParallelGCThreads=$threads\" MarkDuplicates -I mapping/${name}.sort.bam -O mapping/${name}.sort.markdup.bam -M mapping/${name}.sort.markdup.metrics.txt
    samtools index mapping/${name}.sort.markdup.bam && rm mapping/${name}.sort.bam
	""" > mapping_${name}.sh && chmod 755 mapping_${name}.sh
done


rm mapping_rate_stat.xls
ls mapping/*bam | while read i
do
    samtools flagstat $i > ${i}.stat
    name=`echo $i | awk -F '/' '{print $2}' | awk -F '.' '{print $1}'`
    rate=`awk '$4=="mapped"' ${i}.stat | awk '{print $5}' | sed 's/(//g'`
    echo -e "${name}\t${rate}" >> mapping_rate_stat.xls
done

