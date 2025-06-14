### This script uses directory paths defined in configure.env.
### Please ensure the configuration file is correctly set up before running.

library(e1071)
library(caret)
library(pROC)
library(data.table)
library(doParallel)

readRenviron("configure.env")
diagnostic_model_FOLDER <- Sys.getenv("diagnostic_model_FOLDER")
OUTPUT_FOLDER <- Sys.getenv("OUTPUT_FOLDER")

olink_pro <- fread(file.path(diagnostic_model_FOLDER, "olink_protein.csv") , data.table = F)
pro_anno <- fread(file.path(diagnostic_model_FOLDER,"coding143.tsv") , data.table = F)
gene_list <- fread(file.path(diagnostic_model_FOLDER,"gene.list") , data.table = F)
control_list <- fread(file.path(diagnostic_model_FOLDER,"Control.csv") , data.table = F)
Sample <- fread(file.path(diagnostic_model_FOLDER,'sample.list'),data.table=F)

pro_anno$Gene <- sapply(strsplit(pro_anno$meaning, split = ";"), "[", 1)
m <- match(names(olink_pro) , pro_anno$coding)
names(olink_pro) <- c("IID",pro_anno$Gene[m][-1])
olink_pro <- olink_pro[,which(colnames(olink_pro) %in% c("IID",gene_list$Gene))]
olink_pro$num_NA <- rowSums(is.na(olink_pro))
olink_pro <- olink_pro[which(olink_pro$num_NA<7),1:95]
olink_pro[is.na(olink_pro)]<- 0
olink_pro_control <- olink_pro[which(olink_pro$IID%in%control_list$eid),]
set.seed(1111)
cl <- makePSOCKcluster(10)
registerDoParallel(cl)
ctrl <- trainControl(method = "cv",number = 10, classProbs = TRUE,summaryFunction = twoClassSummary)
model_all <- list()
for(i in 1:nrow(Sample))
{
	PHE <- paste0(strsplit(Sample$PHE[i],split="_")[[1]][-length(strsplit(Sample$PHE[i],split="_")[[1]])],collapse="_")
    file_name <- paste0(Sample$PHE[i],'.list')
	case_list <- fread(file.path(diagnostic_model_FOLDER,file_name),data.table=F,header=F)
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
	model_all[[i]] <- list(PHE=Sample$PHE[i], preProcess = glm.model$preProcess,bestTune = glm.model$bestTune,control = glm.model$control,method = glm.model$method,coefnames=glm.model$coefnames)
}
saveRDS(model_all, file = file.path(OUTPUT_FOLDER,"CVD_diagnostic_model.RDS"))

