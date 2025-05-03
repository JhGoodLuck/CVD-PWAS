library(data.table)
genelist <- fread('CVD_gene.list',data.table=F)
pbmc <- readRDS("data/scRNA.RDS")
meta <- pbmc@meta.data
matrix <- matrix[,which(colnames(matrix)%in%row.names(meta))]
matrix <- matrix[which(row.names(matrix)%in%gene_list$gene),]
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
write.table(output,"Type_Difference.genes.logfc.txt",row.names = F,col.names = T,sep = "\t",quote = F)



