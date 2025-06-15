### This script uses directory paths defined in configure.env.
### Please ensure the configuration file is correctly set up before running.
source ./configure.env

##### Using LDSC to Convert GWAS Summary Statistics Format
PHE = CVD               ###Phenotype
N = 10000               ###Total sample size of GWAS
chr = 1                 ###Chromosome
pT = 1.92E-5            ###The threshold of P-value

SOFTWARE_FOLDER/python/python SOFTWARE_FOLDER/ldsc/munge_sumstats.py \
    --sumstats GWAS_FOLDER/${PHE}.ma \
    --N-col N \
    --out GWAS_FOLDER/${GWAS}

#####PWAS analysis#####
Rscript SOFTWARE_FOLDER/FUSION/FUSION.assoc_test.R \
    --sumstats GWAS_FOLDER/${GWAS}.sumstats.gz \
    --weights PWAS_WEIGHT_FOLDER/UKB_Olink.pos \
    --weights_dir PWAS_WEIGHT_FOLDER/UKB_Olink/ \
    --GWASN ${N} \
    --chr ${chr} \
    --ref_ld_chr REF_FOLDER/1000G.EUR. \
    --out OUTPUT_FOLDER/PWAS/${PHE}.${chr}.dat

#####Conditional analysis#####
cat OUTPUT_FOLDER/PWAS/${PHE}.${chr}.dat | awk "$20 < ${pT} || NR!=1{print}" | cut -f1-20 > OUTPUT_FOLDER/Conditional/${PHE}.${chr}.top
Rscript SOFTWARE_FOLDER/FUSION/FUSION.post_process.R \
    --sumstats GWAS_FOLDER/${GWAS}.sumstats.gz \
    --input OUTPUT_FOLDER/Conditional/${PHE}.${chr}.top \
    --ref_ld_chr REF_FOLDER/1000G.EUR.\
    --chr ${chr} --plot \
    --locus_win 1000000 \
    --out OUTPUT_FOLDER/Conditional/${PHE}.${chr}.top.analysis
    
#####SMR analysis#####
cut -f 2 OUTPUT_FOLDER/Conditional/${PHE}.${chr}.top.analysis.joint_included.dat > prob.list
SOFTWARE_FOLDER/SMR/bin/smr-1.3.1 \
    -–bfile REF_FOLDER/1000G.EUR.${chr} \
    --gwas-summary GWAS_FOLDER/${PHE}.ma \
    --beqtl-summary SOFTWARE_FOLDER/SMR/besd/UKB_Olink \
    --extract-probe OUTPUT_FOLDER/SMR/ \
    --maf 0.01 --smr-multi --ld-multi-snp 0.001 \
    --out OUTPUT_FOLDER/SMR/${PHE}.${chr}
