######Obtain the genotype data of the protein
ID=${1}
Pro=`awk -v ID=${ID} '$1==Pro{print $2}' coding143.tsv |cut -f2 |cut -f1 -d ';'`
mkdir -p /WJH/data/UKB_PWAS_weight/data/${Pro}
folder=/WJH/data/UKB_PWAS_weight/data/${Pro}
awk -F ',' -v n=${n} 'NR!=1 && $n!=""{print $1,$1,$n}' ${folder}/olink_protein.csv >> ${folder}/${Pro}_exp.txt
awk '{print $1,$2}' ${folder}/${Pro}_exp.txt > ${folder}/${Pro}_IID.list

chr=`awk -v g=${Pro} '$4==g{print}' /WJH/data/gencode.bed | cut -f1 |sed s#chr##g`
start=`awk -v g=${Pro} '$4==g{print}' /WJH/data/gencode.bed | cut -f2`
end=`awk -v g=${Pro} '$4==g{print}' /WJH/data/gencodebed | cut -f3`
((start=start-500000))
((end=end+500000))

if [ -z ${chr} ]
then
	echo There is no information for protein ${Pro}
	exit
fi

if [ ${start} -lt 0 ]
then
	start=0
fi

if [ ! -f ${folder}/${Pro}.bed ]
then
	plink2 --pfile ${folder}/ukb_imp_chr${chr} \
		--chr ${chr} --from-bp ${start} --to-bp ${end} \
		--keep ${folder}/${Pro}_IID.list \
		--extract ${folder}/w_hm3.snplist \
		--mind 0.05 --maf 0.01 --hwe 1e-8 --make-bed \
		--out ${folder}/${Pro} --force-intersect
fi
#####Estimate the heritability of the protein
/WJH/data/UKB_PWAS_weight/fusion/plink/plink \
	--allow-no-sex \
	--bfile ${folder}/${Pro} \
	--make-grm-bin \
	--out ${folder}/${Pro}_grm

/WJH/data/UKB_PWAS_weight/fusion/gcta_1.93.2beta/gcta64 \
	--grm ${folder}/${Pro}_grm \
	--pheno ${folder}/${Pro}_exp.txt \
	--out ${folder}/${Pro}_reml \
	--reml --reml-no-constrain --reml-lrt 1 --thread-num 15

#####Build the model using FUSION pipeline
fusion_folder=/WJH/data/UKB_PWAS_weight/fusion/
Rscript ${fusion_folder}/FUSION.compute_weights.R \
	--PATH_gcta ${fusion_folder}/gcta/gcta_nr_robust \
	--PATH_plink ${fusion_folder}/plink/plink \
	--PATH_gemma ${fusion_folder}/gemma/gemma-0.98.5-linux-static-AMD64 \
	--bfile ${folder}/${Pro}  \
	--pheno  ${folder}/${Pro}_exp.txt \
	--tmp  ${folder}/tmp/tmp_${Pro} \
	--out /WJH/data/UKB_PWAS_weight/output/${Pro} \
	--models top1,lasso,enet --save_hsq --verbose 1 \
	--covar /WJH/data/UKB_PWAS_weight/data/cov.txt
