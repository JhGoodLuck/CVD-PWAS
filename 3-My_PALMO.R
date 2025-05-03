library("PALMO")
library("Hmisc")
library("ggpubr")
library("cowplot")
library("data.table")
library("circlize")
library("ComplexHeatmap")
library("ggrepel")
library("pbapply")

CVD_list <- fread('CVD_gene.list',data.table=F)
#####bulk protein#####
load("AIFI-Olink-NPX_log2_Protein.Rda")
load("AIFI-Metadata.Rda")
data <- data[which(row.names(data)%in%CVD_list$gene),]
palmo_obj <- createPALMOobject(anndata=ann, data=data)
palmo_obj<- annotateMetadata(data_object=palmo_obj,sample_column= "Sample", donor_column= "PTID",time_column= "Time")
palmo_obj <- mergePALMOdata(data_object=palmo_obj, datatype="bulk")
palmo_obj <- checkReplicates(data_object=palmo_obj, mergeReplicates = T)
palmo_obj <- naFilter(data_object=palmo_obj, na_cutoff=0.4)
featureSet <- c("PTID","Time")
Genelmo_obj <- lmeVariance(data_object=palmo_obj, featureSet=featureSet,meanThreshold=1, fileName="CVDbulk_10")
var_decomp <- palmo_obj@result$variance_decomposition
#plots <- variancefeaturePlot(vardata=var_decomp, featureSet=featureSet,Residual=T)
palmo_obj <- cvCalcBulk(data_object=palmo_obj, meanThreshold=1, cvThreshold=10,fileName="CVDbulk_10",gene_list=CVD_list)


#####sc-RNA#####
pbmc <- readRDS("data/AIFI-scRNA-PBMC-FinalData.RDS")
#CVD_list <- fread('CVD-stable.list',data.table=F)
metaData <- pbmc@meta.data
pbmc@meta.data$Sample <- pbmc@meta.data$orig.ident
pbmc@meta.data$celltype <- gsub(" ", "_", pbmc@meta.data$celltype)
load("data/AIFI-Metadata.Rda")
library("Seurat")
vgGroup <- "celltype"
cell_type <- sort(unique(pbmc@meta.data$celltype))
celltype_oi <- c("CD4-Naive","CD4-TEM","CD4-TCM","CD4-CTL","CD8-Naive","CD8-TEM","CD8-TCM","Treg","MAIT","gdT","NK", "NK-CD56bright","B-naive", "B-memory", "B-intermediate","CD14-Mono","CD16-Mono","cDC2","pDC")
palmo_obj <- createPALMOobject(anndata=ann, data=pbmc)
palmo_obj <- annotateMetadata(data_object=palmo_obj,sample_column= "Sample",donor_column= "PTID",time_column= "Time")
palmo_obj <- mergePALMOdata(data_object=palmo_obj, datatype="singlecell")
palmo_obj <- avgExpCalc(data_object=palmo_obj,assay="RNA", group_column="celltype")
#palmo_obj <- cvCalcSCProfile(data_object=palmo_obj,housekeeping_genes=c("GAPDH", "ACTB"),fileName="CVDsc-10")
palmo_obj <- cvCalcSCProfile(data_object=palmo_obj, meanThreshold = 0.1,housekeeping_genes=c("GAPDH", "ACTB"),fileName="CVDsc-10_all")
#cvSCsampleprofile(data_object=palmo_obj, meanThreshold = 0.1,cvThreshold = 10)
featureSet <- c("PTID", "Time","celltype")
palmo_obj <- lmeVariance(data_object=palmo_obj,featureSet=featureSet,meanThreshold=0.1, cl=4,fileName="CVDsc-10_all")
var_decomp <- palmo_obj@result$variance_decomposition
palmo_obj <- cvCalcSC(data_object=palmo_obj,meanThreshold=0.1, cvThreshold=10,fileName="CVDsc-10_all")
palmo_obj <- StableFeatures(data_object=palmo_obj, group_oi=celltype_oi,cvThreshold=10,donorThreshold=4, groupThreshold=40,topFeatures=25,fileName="CVDsc-10_all",gene_list=CVD_list)

