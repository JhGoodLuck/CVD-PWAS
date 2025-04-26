library(e1071)
library(caret)
library(pROC)
library(data.table)
library(doParallel)

olink_pro <- fread("olink_protein.csv" , data.table = F)
pro_anno <- fread("coding143.tsv" , data.table = F)
gene_list <- fread("gene.list" , data.table = F)
control_list <- fread("Control.csv" , data.table = F)
pro_anno$Gene <- sapply(strsplit(pro_anno$meaning, split = ";"), "[", 1)
m <- match(names(olink_pro) , pro_anno$coding)
names(olink_pro) <- c("IID",pro_anno$Gene[m][-1])
olink_pro <- olink_pro[,which(colnames(olink_pro) %in% c("IID",gene_list$Gene))]
olink_pro$num_NA <- rowSums(is.na(olink_pro))
olink_pro <- olink_pro[which(olink_pro$num_NA<7),1:95]
olink_pro[is.na(olink_pro)]<- 0
Sample <- fread('sample.list',data.table=F)
olink_pro_control <- olink_pro[which(olink_pro$IID%in%control_list$eid),]
set.seed(1111)
cl <- makePSOCKcluster(10)
registerDoParallel(cl)
ctrl <- trainControl(method = "cv",number = 10, classProbs = TRUE,summaryFunction = twoClassSummary)
model_all <- list()
for(i in 1:nrow(Sample))
{
	PHE <- paste0(strsplit(Sample$PHE[i],split="_")[[1]][-length(strsplit(Sample$PHE[i],split="_")[[1]])],collapse="_")
	case_list <- fread(paste0('sample/',Sample$PHE[i],'.list'),data.table=F,header=F)
	dat_case <- olink_pro[which(olink_pro$IID%in%case_list$V1),]
	dat_all <- rbind(dat_case,olink_pro_control)
	dat_all$PHE <- "control"
	dat_all[which(dat_all$IID%in%case_list$V1),"PHE"] <- "case"
	dat_all$PHE <- factor(dat_all$PHE , levels = c("control","case"), labels = c("control","case"))
	set.seed(1111)
	index <- createDataPartition(dat_all$PHE,p = 0.7,list = FALSE)
	dat_train <- dat_all[index,]
	dat_test <- dat_all[-index,]
	svm.model <- train(PHE ~ ., data = dat_train,method = "svmRadial",metric="ROC" , trControl = ctrl)
	pre_prob_svm <- predict(svm.model,newdata = dat_test , type="prob")
	rea <- dat_test[,c("PHE")]
	svm.ROC <- roc(rea , predictor = pre_prob_svm$case , levels = levels(dat_test[,c("PHE")]) , direction = "<")
	
}


