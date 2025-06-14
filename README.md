# CVD-PWAS
This repository contains the primary codes used in the research article:
“Blood plasma proteome-wide association study implicates novel proteins in the pathogenesis of multiple cardiovascular diseases.”
The analyses were conducted using data from the UK Biobank (UKB) cohort and publicly available genome-wide association studies (GWAS) datasets. In addition, we provide a plasma protein genetic prediction model trained on the UKB Pharma Proteomics Project (UKB-PPP) data. This model enables further exploration and hypothesis-driven research into the molecular underpinnings of cardiovascular diseases.

To reproduce the analyses in this project, the following tools and software are required:
| Tool / Package        | Description                                                                                 | Installation                                                |
| --------------------- | --------------------------------------------------------------------------------------------| ----------------------------------------------------------- |
| [FUSION]              | A framework for transcriptome/proteome-wide association studies (TWAS/PWAS)                 | http://gusevlab.org/projects/fusion/                        |
| [PLINK]               | PLINK is a free, open-source whole genome association analysis toolset                      | https://github.com/chrchang/plink-ng                        |
| [LDSC]                | Formatting GWAS summary statistics                                                          | https://github.com/bulik/ldsc                               |
| [SMR]                 | Implements SMR & HEIDI to assess pleiotropic associations using GWAS and pQTL summary data. | https://yanglab.westlake.edu.cn/software/smr/#Overview      |
| [gcta]                |                                                                                             | https://yanglab.westlake.edu.cn/software/gcta/#Overview     |
