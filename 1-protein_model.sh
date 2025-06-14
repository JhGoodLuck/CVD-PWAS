### This script uses directory paths defined in configure.env.
### Please ensure the configuration file is correctly set up before running.
source configure.env

######Obtain the genotype data of the protein
ID=PROTEIN
Pro=`awk -v ID=${ID} '$1==Pro{print $2}'  ${PRO_FOLDER}/coding143.tsv |cut -f2 |cut -f1 -d ';'`
mkdir -p ${PRO_FOLDER}/${Pro}

echo -e "FID\tIID\tEXP" >  ${PRO_FOLDER}/${Pro}_exp.txt
awk -F ',' -v n=${n} 'NR!=1 && $n!=""{print $1,$1,$n}'  ${PRO_FOLDER}/olink_protein.csv >>  ${PRO_FOLDER}/${Pro}_exp.txt
awk '{print $1,$2}'  ${PRO_FOLDER}/${Pro}_exp.txt >  ${PRO_FOLDER}/${Pro}_IID.list

chr=`awk -v g=${Pro} '$4==g{print}'  ${PRO_FOLDER}/gencode.bed | cut -f1 |sed s#chr##g`
start=`awk -v g=${Pro} '$4==g{print}'  ${PRO_FOLDER}/gencode.bed | cut -f2`
end=`awk -v g=${Pro} '$4==g{print}'  ${PRO_FOLDER}/gencode.bed | cut -f3`
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
	${SOFTWARE_FOLDER}/plink2/bin/plink2 \
        --pfile ${GENO_FOLDER}/ukb_imp_chr${chr} \
		--chr ${chr} --from-bp ${start} --to-bp ${end} \
		--keep ${PRO_FOLDER}/${Pro}_IID.list \
		--extract ${GENO_FOLDER}/w_hm3.snplist \
		--mind 0.05 --maf 0.01 --hwe 1e-8 --make-bed \
        --force-intersect \
		--out ${GENO_FOLDER}/${Pro} 
fi

#####Estimate the heritability of the protein using gcta
${SOFTWARE_FOLDER}/plink/bin/plink \
	--allow-no-sex \
	--bfile ${GENO_FOLDER}/${Pro} \
	--make-grm-bin \
	--out ${GENO_FOLDER}/${Pro}_grm

${SOFTWARE_FOLDER}/gcta/bin/gcta64 \
	--grm ${GENO_FOLDER}/${Pro}_grm \
	--pheno ${PRO_FOLDER}/${Pro}_exp.txt \
	--reml --reml-no-constrain \
    --reml-lrt 1 --thread-num 20 \
    --out ${GENO_FOLDER}/${Pro}_reml 

#####Build the model using FUSION pipeline
Rscript ${SOFTWARE_FOLDER}/FUSION/FUSION.compute_weights.R \
	--PATH_gcta ${SOFTWARE_FOLDER}/gcta/bin/gcta64 \
	--PATH_plink ${SOFTWARE_FOLDER}/plink/bin/plink \
	--PATH_gemma ${SOFTWARE_FOLDER}/gemma/bin/gemma-0.98.5-linux-static-AMD64 \
	--bfile ${GENO_FOLDER}/${Pro}  \
	--pheno ${PRO_FOLDER}/${Pro}_exp.txt \
	--tmp ${OUTPUT_FOLDER}/tmp_${Pro} \
	--models top1,lasso,enet --save_hsq --verbose 1 \
	--covar ${GENO_FOLDER}/cov.txt \
    --out ${OUTPUT_FOLDER}/UKB_Olink/${Pro} 
