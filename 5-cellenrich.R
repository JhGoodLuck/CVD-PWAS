### This script uses directory paths defined in configure.env.
### Please ensure the configuration file is correctly set up before running.

library(data.table)

readRenviron("configure.env")
Stable_gene_FOLDER <- Sys.getenv("Stable_gene_FOLDER")
OUTPUT_FOLDER <- Sys.getenv("OUTPUT_FOLDER")

genelist <- fread(file.path(Stable_gene_FOLDER,'CVD_gene.list'),data.table=F)
pbmc <- readRDS(file.path(Stable_gene_FOLDER,"AIFI-scRNA-PBMC-FinalData.RDS"))
celltype_oi <- c("CD4-Naive","CD4-TEM","CD4-TCM","CD4-CTL","CD8-Naive","CD8-TEM","CD8-TCM","Treg","MAIT","gdT","NK", "NK-CD56bright","B-naive", "B-memory", "B-intermediate","CD14-Mono","CD16-Mono","cDC2","pDC")

meta <- pbmc@meta.data
meta <- meta[which(meta$celltype%in%celltype_oi),]
matrix <- pbmc@assays$RNA@counts
matrix <- matrix[,which(colnames(matrix)%in%row.names(meta))]
matrix <- matrix[which(row.names(matrix)%in%gene_list$Gene),]

Run_select <- function(matrix,meta,gene,cell)
{
	counts1 <- matrix[gene,which(meta$celltype == cell), drop = F]
	counts2 <- matrix[gene,which(meta$celltype != cell), drop = F]
	return(list(counts1,counts2))
}
cells <- unique(meta$celltype)
output <- data.frame()
for (i in 1:length(genelist))
{
	gene <- genelist[i]
	print(gene)
	tempout <- data.frame(matrix(nrow = length(cells),ncol = 4))
	colnames(tempout) <- c("gene","cell","logfc","pval")
	tempout$gene <- gene
	tempout$cell <- cells
	for (j in 1:length(cells))
	{
        cell <- cells[j]
        l1 <- Run_select(matrix,meta,gene,cell)[[1]]
		l2 <- Run_select(matrix,meta,gene,cell)[[2]]
		tempout[j,3] <- log2(rowMeans(l1))-log2(rowMeans(l2))
		tempout[j,4] <- wilcox.test(as.numeric(l1),as.numeric(l2),exact = F, correct=T,alternative = "two.sided")$p.value
	}
	output <- rbind(output,tempout)
}
write.table(output,file.path(OUTPUT_FOLDER,"Type_Difference.genes.logfc.txt"),row.names = F,col.names = T,sep = "\t",quote = F)
