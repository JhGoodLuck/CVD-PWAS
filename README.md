# CVD-PWAS
This repository contains the primary codes used in the research article:
“Blood plasma proteome-wide association study implicates novel proteins in the pathogenesis of multiple cardiovascular diseases.”
The analyses were conducted using data from the UK Biobank (UKB) cohort and publicly available genome-wide association studies (GWAS) datasets. In addition, we provide a plasma protein genetic prediction model trained on the UKB Pharma Proteomics Project (UKB-PPP) data. This model enables further exploration and hypothesis-driven research into the molecular underpinnings of cardiovascular diseases.

# Required Tools & Dependencies
To reproduce the analyses in this project, the following tools and software are required:
| Tool / Package        | Description                                                                                 | Installation                                                |
| --------------------- | --------------------------------------------------------------------------------------------| ----------------------------------------------------------- |
| FUSION                | A framework for transcriptome/proteome-wide association studies (TWAS/PWAS)                 | http://gusevlab.org/projects/fusion/                        |
| PLINK2                | PLINK is a free, open-source whole genome association analysis toolset                      | https://www.cog-genomics.org/plink/2.0/                     |
| LDSC                  | Using LDSC to Convert GWAS Summary Statistics Format                                        | https://github.com/bulik/ldsc                               |
| SMR                   | Implements SMR & HEIDI to assess pleiotropic associations using GWAS and pQTL summary data. | https://yanglab.westlake.edu.cn/software/smr/#Overview      |
| GCTA                  | GCTA was used to estimate the heritability of plasma protein abundance based on individual-level genotype data. | https://yanglab.westlake.edu.cn/software/gcta/#Overview     |
| PALMO                 | PALMO is a platform for anayzing longitudinal data from bulk as well as single cell. | https://github.com/aifimmunology/PALMO |
| CARET                 | An R package that simplifies building, tuning, and evaluating machine learning models with a consistent interface across various algorithms. | https://cran.r-project.org/web/packages/caret/index.html  |

# Repository Contents
This repository includes the scripts and data files used in our PWAS analysis pipeline. Below is a brief description of the key components:
| File                      | Description                                   |
|---------------------------|-----------------------------------------------|
| 1-protein_model.sh        | Estimates the SNP-based heritability of plasma proteins using GCTA and builds protein prediction models via the FUSION pipeline, based on data from the UKB-PPP. |
| 2-find_causal_associon.sh | Performs PWAS using the FUSION framework and identifies putatively causal proteins using SMR and HEIDI analysis. |
| 3-CVD_sub_pre.R           | Builds disease diagnostic models using proteins identified as putatively causal for CVDs.   |
| 4-My_PALMO.R              | Evaluates the stability of CVD causal proteins across bulk plasma and scRNA using the PALMO package. |
| 5-cellenrich.R            | Identifies cell type–specific expression patterns of candidate proteins or genes based on single-cell RNA-seq data. |

# Citations
1. Gusev A, Ko A, Shi H, Bhatia G, Chung W, Penninx BW, et al. Integrative approaches for large-scale transcriptome-wide association studies. Nat Genet. 2016;48(3):245-52
2. Chang CC, Chow CC, Tellier LC, Vattikuti S, Purcell SM, Lee JJ. Second-generation PLINK: rising to the challenge of larger and richer datasets. Gigascience. 2015;4:7.
3. Bulik-Sullivan BK, Loh PR, Finucane HK, Ripke S, Yang J, Schizophrenia Working Group of the Psychiatric Genomics C, et al. LD Score regression distinguishes confounding from polygenicity in genome-wide association studies. Nat Genet. 2015;47(3):291-5.
4. Zhu Z, Zhang F, Hu H, Bakshi A, Robinson MR, Powell JE, et al. Integration of summary data from GWAS and eQTL studies predicts complex trait gene targets. Nat Genet. 2016;48(5):481-7.
5. Yang J, Lee SH, Goddard ME, Visscher PM. GCTA: a tool for genome-wide complex trait analysis. Am J Hum Genet. 2011;88(1):76-82.
6. Vasaikar SV, Savage AK, Gong Q, Swanson E, Talla A, Lord C, et al. A comprehensive platform for analyzing longitudinal multi-omics data. Nat Commun. 2023;14(1):1684.
