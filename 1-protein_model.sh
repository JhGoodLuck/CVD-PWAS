######Obtain the genotype data of the protein
Pro=${1}
gene=`awk -v ID=${Pro} '$1==Pro{print $2}' coding143.tsv |cut -f2 |cut -f1 -d ';'`
mkdir -p /WJH/data/UKB_PWAS_weight/data/${gene}
folder=/WJH/data/UKB_PWAS_weight/data/${gene}
awk -F ',' -v n=${n} 'NR!=1 && $n!=""{print $1,$1,$n}' ${folder}/olink_protein.csv >> ${folder}/${gene}_exp.txt
awk '{print $1,$2}' ${folder}/${gene}_exp.txt > ${folder}/${gene}_IID.list

chr=`awk -v g=${gene} '$4==g{print}' /home/WJH/project/WHA/test/gencode.v26.GRCh37.gene.bed | cut -f1 |sed s#chr##g`
start=`awk -v g=${gene} '$4==g{print}' /home/WJH/project/WHA/test/gencode.v26.GRCh37.gene.bed | cut -f2`
end=`awk -v g=${gene} '$4==g{print}' /home/WJH/project/WHA/test/gencode.v26.GRCh37.gene.bed | cut -f3`
((start=start-500000))
((end=end+500000))

if [ -z ${chr} ]
then
	echo There is no information for protein ${gene}
	exit
fi

if [ ${start} -lt 0 ]
then
	start=0
fi

if [ ! -f ${folder}/${gene}.bed ]
then
	plink2 --pfile ${folder}/ukb_imp_chr${chr} \
		--chr ${chr} --from-bp ${start} --to-bp ${end} \
		--keep ${folder}/${gene}_IID.list \
		--extract ${folder}/w_hm3.snplist \
		--mind 0.05 --maf 0.01 --hwe 1e-8 --make-bed \
		--out ${folder}/${gene} --force-intersect
fi

#####Build the model using FUSION pipeline
fusion_folder=/WJH/data/UKB_PWAS_weight/fusion/
Rscript ${fusion_folder}/FUSION.compute_weights.R \
	--PATH_gcta ${fusion_folder}/gcta/gcta_nr_robust \
	--PATH_plink ${fusion_folder}/plink/plink \
	--PATH_gemma ${fusion_folder}/gemma/gemma-0.98.5-linux-static-AMD64 \
	--bfile ${folder}/${gene}  \
	--pheno  ${folder}/${gene}_exp.txt \
	--tmp  ${folder}/tmp/tmp_${gene} \
	--out /WJH/data/UKB_PWAS_weight/output/${gene} \
	--models top1,lasso,enet --save_hsq --verbose 1 \
	--covar /WJH/data/UKB_PWAS_weight/data/cov.txt


