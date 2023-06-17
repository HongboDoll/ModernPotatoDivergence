#!/bin/bash

cd /media/scratchpad_02/li322/work/04_137_atlas_population_analyses/structure

vcf=166_potato_Allchr.gatk.filter.SNP.filterMissMAF.4d.vcf
end=`grep -v '#' ${vcf} | awk '{print NF}'|head -1`
n_sample=`expr $end - 9`
n_loci=`grep -v '#' ${vcf}|wc -l `
ran_loci=10000

rm -rf genotype $(basename ${vcf} .vcf).structure.summary.xls; mkdir genotype
for i in `seq 14 14` # repeat 20 times
do
cat <(grep '#CHROM' ${vcf}) <(grep -v '#' ${vcf}) | ./randomly_select_loci_from_vcf.py $n_loci $ran_loci > $(basename ${vcf} .vcf).genotype${i} ## randomly select 10000 loci from 4d vcf
	for n in `seq 10 $end`
	do
	cut -f 4,5,$n $(basename ${vcf} .vcf).genotype${i} | ./convert_vcf_2_structure.py > genotype/$(basename ${vcf} .vcf)_${n}.structure${i}
	done
cd genotype

cat *structure${i} | sort -k1,1 > $(basename ${vcf} .vcf).structure${i}

cd -
for k in `seq 1 20` # K = 1 ~ 20
do
structure -m mainparams -i genotype/$(basename ${vcf} .vcf).structure${i} -o genotype/$(basename ${vcf} .vcf).structure${i}_K${k}.out -K $k -L $ran_loci -N $n_sample
lnP=`grep 'Estimated Ln Prob of Data' genotype/*.structure${i}_K${k}.out_f | awk '{print $NF}'` 
echo -e "${lnP}\t\c" >> $(basename ${vcf} .vcf).structure.summary.xls${i}
done
echo -e "" >> $(basename ${vcf} .vcf).structure.summary.xls${i}

done

## re-run STRUCTURE using the whole set of SNPs on optimal K value (2-6?)
cat <(grep '#CHROM' ${vcf}) <(grep -v '#' ${vcf}) > $(basename ${vcf} .vcf).genotype

rm -rf struc_tmp; mkdir struc_tmp
for n in `seq 10 $end`
    do
    cut -f 4,5,$n $(basename ${vcf} .vcf).genotype | ./convert_vcf_2_structure.py > struc_tmp/$(basename ${vcf} .vcf)_${n}.structure
done

cat struc_tmp/*structure | sort -k1,1 > $(basename ${vcf} .vcf).structure

for i in `seq 2 14` ### K values
do
    echo """#!/bin/bash
    cd /media/scratchpad_02/li322/work/04_137_atlas_population_analyses/structure
    structure -m mainparams -i $(basename ${vcf} .vcf).structure -o $(basename ${vcf} .vcf).structure.outK${i} -K $i -L $n_loci -N $n_sample
    """ > run.sh${i} && chmod 755 run.sh${i}
    qsub -cwd -N struc${i} -q seq.long.q -j y -V -b y -S /bin/bash "/media/scratchpad_02/li322/work/04_137_atlas_population_analyses/structure/run.sh${i}"
    
done

for k in `seq 2 14`
do
    num=`expr $k + 5`
    grep -A 167 'Inferred ancestry' 166_potato_Allchr.gatk.filter.SNP.filterMissMAF.4d.structure.outK${k}_f|awk '$1~/[1-9]/'|awk -v n=$num '{printf $2" ";for(i=6;i<=n;i++)printf $i""FS;print""}' | sed 's/ /\t/g' > 166_potato_structure_K${k}.xls
done

awk '{print $1}' 166_potato_Allchr.gatk.filter.SNP.filterMissMAF.4d.structure > 166_potato_name.xls
for k in `seq 2 5`
do
    paste 166_potato_name.xls <(cut -f '2-99' 166_potato_structure_K${k}.xls) > t && mv t 166_potato_structure_K${k}.xls
    rm 166_potato_structure_K${k}_order.xls
    cat order | while read i
    do
        awk '$1=="'"$i"'"' 166_potato_structure_K${k}.xls >> 166_potato_structure_K${k}_order.xls
    done
done

