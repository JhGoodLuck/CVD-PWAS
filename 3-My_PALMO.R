library("PALMO")
library("Hmisc")
library("ggpubr")
library("cowplot")
library("data.table")
library("circlize")
library("ComplexHeatmap")
library("ggrepel")
library("pbapply")

cvCalcBulk <- function(data_object, meanThreshold=1,
                       cvThreshold=5,
                       median_cvThreshold=NULL,
                       donorThreshold=NULL,
                       housekeeping_genes=NULL,
                       naThreshold=1,
                       plot_log10=FALSE,
                       selectedFeatures=NULL,
                       median_cv_max=NULL,
                       plotWidth=5, plotHeight=8,
                       fileName=NULL, filePATH=NULL,gene_list=NULL) {
  if(is.null(gene_list)) {
  colnames(gene_list) <- "gene"
  }
  message(date(),": Performing Coefficient of variance analysis")
  if(is.null(fileName)) {
    fileName <- "outputFile"
  }
  if(is.null(filePATH)) {
    filePATH <- data_object@filePATH
  }
  
  ## meanThrehold
  if(is.null(meanThreshold)) {
    meanThreshold <- 0
    message(date(),": Using mean threshold >= 0")
  }
  data_object@meanThreshold <- meanThreshold
  data_object@cvThreshold <- cvThreshold
  
  ## Cumulative cvThrehold (across participants)
  if(is.null(median_cvThreshold)) {
    median_cvThreshold <- cvThreshold
    message(date(),": Using median CV threshold (across donors) same as CV
          threshold at single donor")
    data_object@median_cvThreshold <- median_cvThreshold
  }
  
  ## Assign housekeeping_genes
  if(!is.null(housekeeping_genes)) {
    data_object@housekeeping_genes <- housekeeping_genes
  }
  
  ## get the data
  ann <- data_object@curated$anndata
  mat <- data_object@curated$data
  check_data <- all.equal(row.names(ann), colnames(mat))
  if(check_data == FALSE) {
    stop(date(),": Annotation of samples (rows) and datamatrix columns do
             not match")
  }
  
  ## If features selected before hand
  if(!is.null(selectedFeatures)) {
    message(date(),": Filtering for selected features")
    mat <- mat[selectedFeatures,]
  }
  
  ## CV vs Mean
  unigene <- row.names(mat)
  uniSample <- sort(unique(ann$PTID))
  
  variable_gene <- NA
  stable_gene <- NA
  pdf(paste(filePATH,"/",fileName,"-CV-Sample-Plot.pdf", sep=""),
      width=5, height=5)
  for(i in 1:length(uniSample)) {
    uS <- uniSample[i]
    #print(uS)
    meta_df <- ann[ann$PTID %in% uS,]
    if(nrow(meta_df)>1) {
      df <- mat[unigene, meta_df$Sample]
      df <- data.frame(df, NAs=apply(df,1,function(x){sum(is.na(x))}),
                       Zeros=apply(df,1,function(x){sum(x==0)}),
                       mean=rowMeans(df, na.rm=TRUE),
                       sd=apply(df,1,sd, na.rm=TRUE),
                       var=apply(df,1,var, na.rm=TRUE),
                       stringsAsFactors = FALSE)
      df$CV <- 100*df$sd/df$mean
      dp1 <- df
      dp1 <- dp1[!is.na(dp1$mean),]
      dp1 <- dp1[abs(dp1$mean) >= meanThreshold,]
      
      ## Plot
      plot1 <- ggplot(dp1, aes(x=mean, y=CV)) +
        geom_point(size=0.5, color="grey") +
        labs(title=paste(uniSample[i], " (g=", nrow(dp1),"/",
                         nrow(df),")", sep="")) +
        theme_classic()
      if(plot_log10 ==TRUE) {
        plot1 <- plot1 + scale_x_continuous(trans='log10') +
          scale_y_continuous(trans='log10')
      }
      
      ## Variable genes
      dp2a <- dp1[abs(dp1$CV)> cvThreshold,]
      dp2a <- dp2a[!is.na(dp2a$CV),]
      #dp2a <- dp1 %>% filter(abs(CV)> cvThreshold)
      if(nrow(dp2a)>0) {
        dp2a <- dp2a[order(abs(dp2a$CV), abs(dp2a$mean),
                           decreasing = TRUE),]
        if(nrow(dp2a)>10) {
          dp2a_label <- dp2a[1:10,]
        } else {
          dp2a_label <- dp2a
        }
        plot1 <- plot1 +geom_text_repel(data=dp2a_label,
                                        aes(x=mean, y=CV, label=row.names(dp2a_label)),
                                        col="#D32B2B", size=2, max.overlaps=20)
        
        ## Add variable genes
        variable_gene <- rbind(variable_gene,
                               data.frame(donor=uS,
                                          feature=row.names(dp2a),
                                          dp2a[,c("mean","sd","var", "CV", "NAs",
                                                  "Zeros")],
                                          stringsAsFactors = FALSE))
      }
      
      ## Stable features
      dp2b <- dp1[abs(dp1$CV) <= cvThreshold,]
      #dp2b <- dp1 %>% filter(abs(CV) <= cvThreshold)
      dp2b <- dp2b[!is.na(dp2b$CV),]
      if(nrow(dp2b)>0) {
        dp2b <- dp2b[order(-abs(dp2b$mean), abs(dp2b$CV),
                           decreasing = FALSE),]
        if(nrow(dp2b)>10) {
          dp2b_label <- dp2b[1:10,]
        } else {
          dp2b_label <- dp2b
        }
        plot1 <- plot1 + geom_text_repel(data=dp2b_label,
                                         aes(x=mean, y=CV,
                                             label=row.names(dp2b_label)),
                                         col="#4A92E2", size=2, max.overlaps=20)
        
        #Add stable genes
        stable_gene <- rbind(stable_gene,
                             data.frame(donor=uS,
                                        feature=row.names(dp2b),
                                        dp2b[,c("mean","sd","var","CV",
                                                "NAs", "Zeros")],
                                        stringsAsFactors = FALSE))
        
      }
      #print(plot1)
      
      ## Heatmap for variance
      #df <- df[abs(df$mean) >= meanThreshold,]
      temp <- data.frame(feature=row.names(df), var=df$CV,
                         stringsAsFactors = FALSE)
      colnames(temp)[-1] <- paste(uS, "_", colnames(temp)[-1], sep="")
      if(i==1) {
        res <- temp
      } else {
        res <- merge(res, temp, by="feature", all=TRUE)
      }
    }
    
  }
  dev.off()
  
  ## Remove blank
  if(!is.null(nrow(variable_gene))) {
    variable_gene <- variable_gene[!is.na(variable_gene$donor),]
  }
  if(!is.null(nrow(stable_gene))) {
    stable_gene <- stable_gene[!is.na(stable_gene$donor),]
  }
  
  ## Decompose result into mean and CV (variance)
  if(ncol(res) == 2) {
    res_var <- data.frame(res[,-1], stringsAsFactors = FALSE)
    colnames(res_var) <- colnames(res)[2]
    row.names(res_var) <- res$feature
    colnames(res_var) <- gsub("_var", "", colnames(res_var))
    var_mat <- NA
    stable_mat <- NA
    
  } else if(ncol(res) > 2){
    res_var <- res[,-1]
    row.names(res_var) <- res$feature
    colnames(res_var) <- gsub("_var", "", colnames(res_var))
    
    #May be some samples are bad
    uniSample <- intersect(uniSample, colnames(res_var))
    
    #If donor cutoff NULL then use all donors
    if(is.null(donorThreshold)) {
      #Calculate feature Mean of CV by all donors
      res_var$Mean <- apply(res_var[,uniSample],1,function(x){
        mean(abs(x),na.rm=TRUE)
      })
      #Calculate feature Median of CV by all donors
      res_var$Median <- apply(res_var[,uniSample],1,function(x){
        median(abs(x),na.rm=TRUE)
      })
    } else {
      
      if(donorThreshold > length(uniSample)) {
        message(date(),": donorThreshold greater than number of donors")
      }
      #Calculate feature Mean of CV by donor threshold
      res_var$Mean <- apply(res_var[,uniSample],1,function(x){
        x <- sort(x)[1:donorThreshold]
        mean(abs(x),na.rm=TRUE)
      })
      #Calculate feature Median of CV by donor threshold
      res_var$Median <- apply(res_var[,uniSample],1,function(x){
        x <- sort(x)[1:donorThreshold]
        median(abs(x),na.rm=TRUE)
      })
    }
    
    #Calculate number of NAs in all donors
    res_var$NAs <- apply(res_var[,uniSample],1,function(x){
      sum(is.na(x),na.rm=TRUE)
    })
    
    ## Remove genes with CV cutoff (high noise)
    if(!is.null(median_cv_max)) {
      res_var <- res_var[res_var$Median < median_cv_max,]
    }
    
    ## Calculate CV range and check whether greater than median
    max_cv <- max(abs(res_var$Median), na.rm=TRUE)
    if(median_cvThreshold> max_cv) {
      median_cvThreshold <- max_cv
      message(date(),": Input median_cv_max higher than maximum CV
                    range.")
    }
    
    #Remove NAs
    unisample_naThreshold <- ceiling(length(uniSample)*naThreshold)
    res_var <- res_var[res_var$NAs <= unisample_naThreshold,]
    
    #histogram of CV
    suppressMessages(dx <- melt(res_var[,uniSample]), classes = "message")
    plot1 <- ggplot(dx, aes(x=value)) +
      geom_histogram(aes(y=..density..), colour="black", fill="skyblue",
                     bins=50)+
      labs(x="CV") +
      theme_classic()
    pdf(paste(filePATH,"/",fileName,"-CV-distribution.pdf", sep=""),
        width=5, height=5)
    #print(plot1)
    dev.off()
    #print(plot1)
    
    #color list
    if(min(res_var$Median,na.rm=TRUE)<0) {
      col_list <- colorRamp2(c(-max_cv, -cvThreshold-0.01, -cvThreshold,
                               -cvThreshold/2, 0, cvThreshold/2,
                               cvThreshold, cvThreshold+0.01, max_cv),
                             c("red", "pink", "blue", "skyblue", "black",
                               "skyblue", "blue", "pink", "red"))
    } else {
      col_list <- colorRamp2(c(0, cvThreshold/2, cvThreshold,
                                cvThreshold+(max_cv-cvThreshold)/2, max_cv),
                             c("#4A92E2","#86B6EB", "white","#E17171","#D32B2B"))
    }
    
    ## Plot variable genes CV
    mat <- res_var[order(-(res_var$NAs), abs(res_var$Median),
                         decreasing = TRUE),]
    mat <- mat[abs(mat$Median) > median_cvThreshold,]
    if(nrow(mat)>1) {
      var_mat <- mat
      mat <- mat[,uniSample]
      if(nrow(mat)>20) {
        #mat <- mat[1:50,]
        write.csv(mat, file=paste(filePATH,"/",fileName,
                                  "-CV-Variable-Matrix.csv", sep=""))
      }
      #Avoid boxplot error
      column_mat <- mat[,uniSample]
      column_mat[is.na(column_mat)] <- 0
      check_0 <- unique(as.numeric(unlist(apply(column_mat,2,
                                                function(x) {
                                                  unique(x,na.rm=TRUE)
                                                }) )))
      if(length(check_0) != 1) {
        column_ha <- HeatmapAnnotation(CV = anno_boxplot(column_mat))
        ht1 <- Heatmap(data.matrix(mat), cluster_rows =FALSE,
                       cluster_columns = FALSE,
                       na_col = "grey", col = col_list,
                       row_names_max_width=unit(10, "cm"),
                       column_names_gp = gpar(fontsize = 5),
                       row_names_gp = gpar(fontsize = 6),
                       top_annotation = column_ha,
                       column_title="Variable features",
                       heatmap_legend_param = list(title = "median CV(%)",
                                                   heatmap_legend_side = "right") )
      } else {
        ht1 <- Heatmap(data.matrix(mat), cluster_rows =FALSE,
                       cluster_columns = FALSE,
                       na_col = "grey", col = col_list,
                       row_names_max_width=unit(10, "cm"),
                       column_names_gp = gpar(fontsize = 5),
                       row_names_gp = gpar(fontsize = 6),
                       column_title="Stable features",
                       heatmap_legend_param = list(title = "median CV(%)",
                                                   heatmap_legend_side = "right") )
      }
      #print(ht1)
      pdf(paste(filePATH,"/",fileName,"-CV-Variable-Heatplot.pdf",
                sep=""), width=plotWidth, height=plotHeight)
      draw(ht1)
      dev.off()
    } else {
      var_mat <- NA
      message(date(),": Variable features do not found. Check data
            for missing values or change CV parameters")
    }
    
    ## Plot stable genes CV
    mat <- res_var[order(res_var$NAs, abs(res_var$Median),
                         decreasing = FALSE),]
    mat <- mat[abs(mat$Median) <= median_cvThreshold,]
    m <- match(gene_list$gene,row.names(mat))
    m <- m[!is.na(m)]
    mat <- mat[m,]
    if(nrow(mat)>1) {
      stable_mat <- mat
      mat <- mat[,uniSample]
      if(nrow(mat)>20) {
        #mat <- mat[,]
        write.csv(mat, file=paste(filePATH,"/",fileName,
                                  "-CV-Stable-Matrix.csv", sep=""))
      }
      #Avoid boxplot error
      column_mat <- mat[,uniSample]
      column_mat[is.na(column_mat)] <- 0
      check_0 <- unique(as.numeric(unlist(apply(column_mat,2,
                                                function(x) {
                                                  unique(x,na.rm=TRUE)
                                                }) )))
      if(length(check_0) != 1) {
        column_ha <- HeatmapAnnotation(CV = anno_boxplot(column_mat))
        ht2 <- Heatmap(data.matrix(mat), cluster_rows =FALSE,
                       cluster_columns = FALSE,
                       na_col = "grey", col = col_list,
                       row_names_max_width=unit(10, "cm"),
                       column_names_gp = gpar(fontsize = 5),
                       row_names_gp = gpar(fontsize = 6),
                       top_annotation = column_ha,
                       column_title="Stable features",
                       heatmap_legend_param = list(title = "median CV(%)",
                                                   heatmap_legend_side = "right") )
      } else {
        ht2 <- Heatmap(data.matrix(mat), cluster_rows =FALSE,
                       cluster_columns = FALSE,
                       na_col = "grey", col = col_list,
                       row_names_max_width=unit(10, "cm"),
                       column_names_gp = gpar(fontsize = 5),
                       row_names_gp = gpar(fontsize = 6),
                       column_title="Stable features",
                       heatmap_legend_param = list(title = "median CV(%)",
                                                   heatmap_legend_side = "right") )
      }
      
      #print(ht2)
      pdf(paste(filePATH,"/",fileName,"-CV-Stable-Heatplot.pdf", sep=""),
          width=plotWidth, height=plotHeight)
      draw(ht2)
      dev.off()
    } else {
      stable_mat <- NA
      message(date(),": Stable features do not found. Check data for
            missing values or change CV parameters. Summary of CV range
                    in given data:")
      message(summary(res_var$Median))
    }
    
    ## Housekeeping genes
    if(!is.null(housekeeping_genes)) {
      housekeeping_genes <- intersect(housekeeping_genes,
                                      row.names(res_var))
      mat <- res_var[housekeeping_genes,uniSample]
      #Avoid boxplot error
      column_mat <- mat[,uniSample]
      column_mat[is.na(column_mat)] <- 0
      check_0 <- unique(as.numeric(unlist(apply(column_mat,2,
                                                function(x) {unique(x,na.rm=TRUE) }) )))
      if(length(check_0) != 1) {
        column_ha <- HeatmapAnnotation(CV = anno_boxplot(column_mat))
        ht3 <- Heatmap(data.matrix(mat), cluster_rows =FALSE,
                       cluster_columns = FALSE,
                       na_col = "grey", col = col_list,
                       row_names_max_width=unit(10, "cm"),
                       column_names_gp = gpar(fontsize = 5),
                       row_names_gp = gpar(fontsize = 6),
                       top_annotation = column_ha,
                       column_title="Housekepping features",
                       heatmap_legend_param = list(title = "median CV(%)",
                                                   heatmap_legend_side = "right"))
      } else {
        ht3 <- Heatmap(data.matrix(mat), cluster_rows =FALSE,
                       cluster_columns = FALSE,
                       na_col = "grey", col = col_list,
                       row_names_max_width=unit(10, "cm"),
                       column_names_gp = gpar(fontsize = 5),
                       row_names_gp = gpar(fontsize = 6),
                       column_title="Housekepping features",
                       heatmap_legend_param = list(title = "median CV(%)",
                                                   heatmap_legend_side = "right"))
      }
      
      #print(ht3)
      pdf(paste(filePATH,"/",fileName,
                "-CV-housekeeping_genes-Heatplot.pdf", sep=""),
          width=plotWidth, height=plotHeight)
      draw(ht3)
      dev.off()
    }
    
  }
  
  #All CV
  write.csv(res_var, file=paste(filePATH,"/",fileName,"-CV-result.csv",
                                sep=""))
  rm(res, temp)
  #Variable genes
  write.csv(variable_gene, file=paste(filePATH,"/",
                                      fileName,"-CV-VariableFeatures-result.csv", sep=""))
  #Stable genes
  write.csv(stable_gene, file=paste(filePATH,"/",
                                    fileName,"-CV-StableFeatures-result.csv", sep=""))
  
  ## Add CV result
  data_object@result$cv_all <- res_var
  data_object@result$variable_gene <- variable_gene
  data_object@result$non_variable_gene <- stable_gene
  data_object@result$var_mat <- var_mat
  data_object@result$stable_mat <- stable_mat
  
  message(date(),": Done. Please check output directory for Plots/results.")
  return(data_object)
}

avgExpCalc <- function(data_object, assay = "RNA", group_column) {
  
  message(date(), ": Calculating scRNA Average expression")
  
  ## Get the data
  if (!is.null(data_object@curated$SeuratObj)) {
    anndata <- data_object@curated$anndata
    dataObj <- data_object@curated$SeuratObj
    
    ## Add sample group to metadata Define sample group and Calculate
    ## average expression
    dataObj@meta.data$group <- dataObj@meta.data[, group_column]
    dataObj@meta.data$Sample_group <- paste(dataObj@meta.data$Sample,
                                            dataObj@meta.data$group,
                                            sep = ":")
    dataObj@meta.data$Sample_group <- gsub(" ", "_", dataObj@meta.data$Sample_group)
    metaData <- dataObj@meta.data
    
    ## Check assay
    DefaultAssay(dataObj) <- assay
    ## Average expression on log-scaled data
    dataObj@assays[[assay]]@counts <- dataObj@assays[[assay]]@data
    scrna_avgmat <- AverageExpression(object = dataObj, assays = assay,
                                      slot = "counts",
                                      group.by = "Sample_group",
                                      verbose = TRUE)
    cn <- data.frame(colnames(scrna_avgmat[[assay]]))
    #cn <- gsub("data\\[, 1]", "", as.character(row.names(cn)))
    scrna_avgmat <- data.frame(scrna_avgmat[[assay]], check.names = FALSE,
                               stringsAsFactors = FALSE)
    colnames(scrna_avgmat) <- cn[,1]
    message(date(), ": scRNA Average expression finished")
    
    ## Keep genes with avgExpression > zero
    rowDF <- rowSums(scrna_avgmat)
    rowDF <- rowDF[rowDF > 0]
    mat <- scrna_avgmat[names(rowDF), ]
    message(date(), ": Keeping genes with avg expression >0")
    
    ## Create annotation
    cn <- data.frame(Sample_group = colnames(mat))
    #print(head(cn))
    temp <- data.frame(do.call(rbind, strsplit(cn$Sample_group,
                                               split = "_")),
                       stringsAsFactors = FALSE)
    #print(head(temp))
    cn <- data.frame(cn, Sample = temp$X2, group = temp$X1,
                     stringsAsFactors = FALSE)
    
    row.names(cn) <- cn$Sample_group
    cn <- merge(cn, anndata, by = "Sample", all = TRUE)
    cn <- cn[!is.na(cn$Sample_group), ]
    row.names(cn) <- cn$Sample_group
    ann <- cn
    ann[, ncol(ann) + 1] <- ann$group
    colnames(ann) <- c(colnames(ann)[1:ncol(ann) - 1], group_column)
    ann$Sample_group_i <- paste(ann$group, ann$PTID, sep = ":")
    rm(cn)
    
    ## Add CV result
    data_object@curated$anndata <- ann
    data_object@curated$data <- mat
    data_object@rownames <- row.names(mat)
    data_object@colnames <- colnames(mat)
    return(data_object)
    
  } else {
    stop(date(), ": Seurat object is absent")
  }
  
  
}

cvCalcSC <- function(data_object, meanThreshold = NULL, cvThreshold = NULL,
                     housekeeping_genes = NULL, cl = 2, fileName = NULL,
                     filePATH = NULL) {
  
  message(date(), ": Performing Coefficient of variance analysis")
  
  ## If filename or filepath null
  if (is.null(fileName)) {
    fileName <- "outputFile"
  }
  if (is.null(filePATH)) {
    filePATH <- data_object@filePATH
  }
  
  ## Assign housekeeping_genes
  if (is.null(housekeeping_genes)) {
    housekeeping_genes <- c("ACTB", "GAPDH")
    data_object@housekeeping_genes <- housekeeping_genes
  }
  
  ## meanThrehold
  if (is.null(meanThreshold)) {
    meanThreshold <- 0
    message(date(), ": Using mean threshold >= 0")
  }
  data_object@meanThreshold <- meanThreshold
  
  ## cvThreshold
  if (is.null(cvThreshold)) {
    cvThreshold <- 10
    message(date(), ": Using cv threshold 10")
  }
  data_object@cvThreshold <- cvThreshold
  
  ## get the data
  ann <- data_object@curated$anndata
  mat <- data_object@curated$data
  check_data <- all.equal(row.names(ann), colnames(mat))
  if (check_data == FALSE) {
    stop(date(), ": Annotation of samples (rows) and datamatrix columns do
             not match")
  }
  
  ## Define group
  ann$group_donor <- paste(ann$group, ann$PTID, sep = ":")
  sample_freq <- data.frame(table(ann$group_donor))
  sample_freq <- sample_freq[sample_freq$Freq > 1, ]
  sample_freq <- as.character(sample_freq$Var1)
  if (length(sample_freq) > 0) {
    ann <- ann[ann$group_donor %in% sample_freq, ]
    mat <- mat[, row.names(ann)]
  } else {
    stop(date(), ": Not enough group samples to perform CV analysis\n")
  }
  
  # Calculate CV vs Mean for all genes per celltype
  unigene <- row.names(mat)
  uniDonor <- sort(unique(ann$PTID))
  uniSamplegroup <- as.character(unique(ann$group_donor))
  
  # All genes CV calculations
  message(date(), ": Performing CV calculations")
  op <- pboptions(type = "timer")  # default
  res <- pblapply(uniSamplegroup, cl = cl, function(uS) {
    # print(uS)
    ann_df <- ann[ann$group_donor %in% uS, ]
    df <- mat[unigene, ann_df$Sample_group]
    df <- data.frame(df,
                     na = apply(df, 1, function(x) {sum(x != 0) }),
                     mean = rowMeans(df, na.rm = TRUE),
                     sd = apply(df, 1, sd, na.rm = TRUE),
                     var = apply(df, 1, var, na.rm = TRUE),
                     stringsAsFactors = FALSE)
    df$cv <- 100 * df$sd/df$mean
    return(df$cv)
  })
  pboptions(op)
  cv_res <- do.call(cbind, res)
  row.names(cv_res) <- unigene
  colnames(cv_res) <- uniSamplegroup
  cv_res <- data.frame(cv_res, check.names = FALSE, stringsAsFactors = FALSE)
  rm(res)
  # save result
  save(cv_res, file = paste(filePATH, "/", fileName, "-CV-allgenes-raw.Rda",
                            sep = ""))
  cv_all <- cv_res
  
  # Genes with minimum mean threshold CV calculations
  message(date(), ": Checking Mean Threshold")
  op <- pboptions(type = "timer")  # default
  res <- pblapply(uniSamplegroup, cl = cl, function(uS) {
    # print(uS)
    ann_df <- ann[ann$group_donor %in% uS, ]
    df <- mat[unigene, ann_df$Sample_group]
    df <- data.frame(df,
                     na = apply(df, 1, function(x) {sum(x != 0)}),
                     mean = rowMeans(df,na.rm = TRUE),
                     sd = apply(df, 1, sd, na.rm = TRUE),
                     var = apply(df, 1, var, na.rm = TRUE),
                     stringsAsFactors = FALSE)
    # CV <- 100*df$sd/df$mean return(CV)
    df$cv <- 100 * df$sd/df$mean
    df$cv <- ifelse(df$mean >= meanThreshold, df$cv, NA)
    return(df$cv)
  })
  pboptions(op)
  cv_res <- do.call(cbind, res)
  row.names(cv_res) <- unigene
  colnames(cv_res) <- uniSamplegroup
  cv_res <- data.frame(cv_res, check.names = FALSE, stringsAsFactors = FALSE)
  rm(res)
  # save result
  save(cv_res, file = paste(filePATH, "/", fileName, "-CV-allgenes.Rda",
                            sep = ""))
  
  # Variable genes
  message(date(), ": Performing Variable gene CV analysis")
  op <- pboptions(type = "timer")  # default
  res <- pblapply(uniSamplegroup, cl = cl, function(uS) {
    ann_df <- ann[ann$group_donor %in% uS, ]
    df <- mat[unigene, ann_df$Sample_group]
    df <- data.frame(df,
                     nonZero = apply(df, 1, function(x) {sum(x != 0)}),
                     mean = rowMeans(df, na.rm = TRUE),
                     sd = apply(df, 1, sd, na.rm = TRUE),
                     var = apply(df, 1, var, na.rm = TRUE),
                     stringsAsFactors = FALSE)
    df$CV <- 100 * df$sd/df$mean
    # the CV becomes very high for data with 0
    df <- df[df$mean >= meanThreshold, ]  #minimum expression >2^0.1=1
    dp2a <- df[df$mean >= meanThreshold & df$CV > cvThreshold,
               c("mean", "sd", "var", "CV")]
    dp2a <- dp2a[order(dp2a$CV, dp2a$mean, decreasing = TRUE), ]
    # Find variable genes
    if (nrow(dp2a) >= 1) {
      variable_gene <- data.frame(donor = uS, gene = row.names(dp2a),
                                  dp2a, stringsAsFactors = FALSE)
      return(variable_gene)
    }
  })
  pboptions(op)
  variable_gene <- do.call(rbind, res)
  rm(res)
  # save result
  save(variable_gene, file = paste(filePATH, "/", fileName,
                                   "-CV-Variablegene.Rda", sep = ""))
  
  # Stable genes
  message(date(), ": Performing Stable gene CV analysis")
  op <- pboptions(type = "timer")  # default
  res <- pblapply(uniSamplegroup, cl = cl, function(uS) {
    ann_df <- ann[ann$group_donor %in% uS, ]
    df <- mat[unigene, ann_df$Sample_group]
    df <- data.frame(df,
                     nonZero = apply(df, 1, function(x) { sum(x != 0) }),
                     mean = rowMeans(df, na.rm = TRUE),
                     sd = apply(df, 1, sd, na.rm = TRUE),
                     var = apply(df, 1, var, na.rm = TRUE),
                     stringsAsFactors = FALSE)
    df$CV <- 100 * df$sd/df$mean
    # the CV becomes very high for data with 0
    df <- df[df$mean >= meanThreshold, ]  #minimum expression >2^0.1=1
    dp2b <- df[df$mean >= meanThreshold & df$CV <= cvThreshold,
               c("mean", "sd", "var", "CV")]
    dp2b <- dp2b[order(-dp2b$mean, dp2b$CV, decreasing = FALSE), ]
    # Find stable genes
    if (nrow(dp2b) >= 1) {
      non_variable_gene <- data.frame(donor = uS, gene = row.names(dp2b),
                                      dp2b, stringsAsFactors = FALSE)
      return(non_variable_gene)
    }
  })
  pboptions(op)
  non_variable_gene <- do.call(rbind, res)
  rm(res)
  # save result
  save(non_variable_gene, file = paste(filePATH, "/",fileName,
                                       "-CV-nonVariablegene.Rda", sep = ""))
  
  # Housekeeping genes data
  message(date(), ": Checking Housekeeping-genes CV")
  op <- pboptions(type = "timer")  # default
  res <- pblapply(uniSamplegroup, cl = cl, function(uS) {
    ann_df <- ann[ann$group_donor %in% uS, ]
    df <- mat[unigene, ann_df$Sample_group]
    df <- data.frame(df,
                     nonZero = apply(df, 1, function(x) { sum(x != 0) }),
                     mean = rowMeans(df, na.rm = TRUE),
                     sd = apply(df, 1, sd, na.rm = TRUE),
                     var = apply(df, 1, var, na.rm = TRUE),
                     stringsAsFactors = FALSE)
    df$CV <- 100 * df$sd/df$mean
    # the CV becomes very high for data with 0
    df <- df[df$mean >= meanThreshold, ]  #minimum expression >2^0.1=1
    temp <- df[housekeeping_genes, c("mean", "sd", "var", "CV")]
    thr <- round(1 + max(df$CV, na.rm = TRUE))
    # Housekeeping gene profile
    temp <- data.frame(gene = row.names(temp), donor = uS, temp,
                       stringsAsFactors = FALSE)
    return(temp)
  })
  pboptions(op)
  hg_res <- do.call(rbind, res)
  rm(res)
  # Housekeeping genes
  hg_res <- hg_res[!is.na(hg_res$mean), ]
  
  # Variable genes
  rn <- data.frame(do.call(rbind, strsplit(variable_gene$donor,
                                           split = ":")),stringsAsFactors = FALSE)
  variable_gene$sample <- rn$X2
  variable_gene$group <- rn$X1
  plot2 <- ggplot(variable_gene, aes(x = mean, y = CV)) +
    geom_point() +facet_wrap(~group)
  
  # stable/non-variable genes
  rn <- data.frame(do.call(rbind, strsplit(non_variable_gene$donor,
                                           split = ":")), stringsAsFactors = FALSE)
  non_variable_gene$sample <- rn$X2
  non_variable_gene$group <- rn$X1
  plot3 <- ggplot(non_variable_gene, aes(x = mean, y = CV)) +
    geom_point() + facet_wrap(~group)
  
  if (nrow(hg_res) > 0) {
    # house-keeping genes
    rn <- data.frame(do.call(rbind, strsplit(hg_res$donor,
                                             split = ":")), stringsAsFactors = FALSE)
    hg_res$sample <- rn$X2
    hg_res$group <- rn$X1
    
    plot4 <- ggplot(variable_gene, aes(x = mean, y = CV)) +
      geom_point() +
      geom_point(data=non_variable_gene, aes(x=mean, y=CV), color="red") +
      geom_point(data=hg_res, aes(x=mean, y=CV), color="blue") +
      geom_hline(yintercept=cvThreshold) +
      facet_wrap(~group, scales="free_y")
  } else {
    plot4 <- ggplot(variable_gene, aes(x = mean, y = CV)) +
      geom_point() +
      geom_point(data=non_variable_gene, aes(x=mean, y=CV), color="red") +
      geom_hline(yintercept=cvThreshold) +
      facet_wrap(~group, scales="free_y")
  }
  
  message(date(), ": Saving CV plots in output directory")
  png(paste(filePATH, "/", fileName, "-CV-distribution.png", sep = ""),
      width = 10, height = 10, res = 200, units = "in")
  print(plot4)
  dev.off()
  
  ## Add CV result
  data_object@result$cv_all <- cv_all
  data_object@result$cv_meanthreshold <- cv_res  #minimum mean threshold
  data_object@result$variable_gene <- variable_gene
  data_object@result$non_variable_gene <- non_variable_gene
  
  message(date(), ": Done. Please check output directory for results.")
  return(data_object)
}


StableFeatures <- function(data_object, group_oi=NULL,
                           cvThreshold=NULL, donorThreshold=NULL,
                           housekeeping_genes=NULL,
                           groupThreshold=NULL, topFeatures=25,
                           filePATH=NULL, fileName=NULL,gene_list=NULL) {
  
  message(date(),": Identifying Stable features")
  
  if(is.null(gene_list)) {
    colnames(gene_list) <- "gene"
  }
  ## If filename or filepath null
  if(is.null(fileName)) {
    fileName <- "outputFile"
  }
  if(is.null(filePATH)) {
    filePATH <- data_object@filePATH
  }
  
  #Load data
  if(!is.null(data_object@result$cv_meanthreshold)) {
    cv_res <- data_object@result$cv_meanthreshold
  } else {
    stop(date(),": Please run cvCalcSC function before StableFeatures
             function\n")
  }
  variable_gene <- data_object@result$variable_gene
  non_variable_gene <- data_object@result$non_variable_gene
  non_variable_gene <- non_variable_gene[which(non_variable_gene$gene %in% gene_list$gene),]
  ## get the annotation data
  data_object@curated$anndata$Sample_group_i <- paste(data_object@curated$anndata$group,
                                                      data_object@curated$anndata$PTID, sep=":")
  ann <- data_object@curated$anndata
  
  #group list
  if(is.null(group_oi)) {
    group_oi <- unique(as.character(ann$group))
  }
  
  #Select group of interest
  ann_sub <- ann[ann$group %in% group_oi,]
  if(nrow(ann_sub)<1) {
    stop(date(), ": Group of interest features do not match with annotation
             group eg.",unique(ann$group)[1:3])
  }
  Sample_group <- unique(ann_sub$Sample_group_i)
  Sample_group <- intersect(Sample_group, colnames(cv_res))
  
  if(is.null(donorThreshold)) {
    donorThreshold <- length(unique(ann_sub$PTID))
    message(date(),": Donor threshold defined ",donorThreshold)
  } else if(donorThreshold > length(unique(ann$PTID))) {
    donorThreshold <- length(unique(ann_sub$PTID))
    message(date(),": Donors were larger than unique donors. Donor threshold
                defined ",donorThreshold)
  }
  
  gThr <- round(length(unique(ann_sub$PTID)) *
                  length(unique(ann_sub$group)) * 0.9)
  if(is.null(groupThreshold)) {
    groupThreshold <- round(length(unique(ann_sub$PTID)) *
                              length(unique(ann_sub$group)) * 0.5)
    message(date(),": Groupwise threshold defined ",groupThreshold)
  } else if(groupThreshold > gThr) {
    groupThreshold <- round(length(unique(ann_sub$PTID)) *
                              length(unique(ann_sub$group)) * 0.9)
    message(date(),": Number of groups were larger than unique
                donor x group. Groupwise threshold defined ",groupThreshold)
  }
  
  #Summary of non-variable genes
  temp <- data.frame(do.call(rbind,
                             strsplit(as.character(non_variable_gene$donor), split = ":")),
                     stringsAsFactors = FALSE)
  non_variable_gene$PTID <- temp$X2
  non_variable_gene$group <- temp$X1
  stable_genelist <- non_variable_gene[non_variable_gene$group %in% group_oi,]
  stable_genelist <- data.frame(table(stable_genelist$gene))
  stable_genelist <- stable_genelist[order(stable_genelist$Freq,
                                           decreasing = TRUE),]
  stable_gene <- as.character(stable_genelist$Var1)
  #create Stable matrix
  stable_gene <- intersect(stable_gene, row.names(cv_res))
  stable_mat <- cv_res[stable_gene,]
  #stable_mat[stable_mat > cvThreshold] <- NA
  save(stable_mat, file=paste(filePATH,"/",fileName,"-stableMatrix.Rda",
                              sep=""))
  
  #Define the super-stable genes
  plot1 <- ggplot(stable_genelist, aes(x=Freq)) +
    geom_histogram(binwidth=1) +
    labs(title="Stable genes occurance from each sample")
  ## Atleast in all donor x group
  super_stable1 <- stable_genelist[stable_genelist$Freq >= groupThreshold,]
  super_stable2 <- stable_genelist[stable_genelist$Freq >= donorThreshold, ]
  #plot heatmap (super-stable)
  gn <- as.character(super_stable1$Var1)
  data_mat <- stable_mat[gn,Sample_group]
  rn <- data.frame(do.call(rbind, strsplit(colnames(data_mat), split = ":")),
                   stringsAsFactors = FALSE)
  #set.seed(2020)
  ha_col <- HeatmapAnnotation(df=data.frame(group=rn$X1),
                              annotation_name_gp = gpar(fontsize = 6),
                              simple_anno_size = unit(0.3, "cm"))
  ht1 <- Heatmap(data.matrix(data_mat), cluster_rows =FALSE,
                 cluster_columns = FALSE,
                 column_split = factor(rn$X1, levels = group_oi),
                 na_col = "grey",
                 col = colorRamp2(c(0,cvThreshold,(cvThreshold+0.001),
                                    (cvThreshold+10)), c("blue","white","pink","red")),
                 row_names_max_width=unit(10, "cm"),
                 column_names_gp = gpar(fontsize = 5),
                 row_names_gp = gpar(fontsize = 6),
                 row_title = "Super stable 25",
                 column_title_gp = gpar(fontsize = 4),
                 top_annotation = ha_col,
                 heatmap_legend_param = list(title = "CV",
                                             heatmap_legend_side = "right") )
  pdf(paste(filePATH,"/",fileName,"-Super-stable-Top25.pdf", sep=""),
      width=10, height=3.5)
  print(plot1)
  print(ht1)
  dev.off()
  print(ht1)
  write.csv(data_mat, paste(filePATH,"/",fileName,"-Super-stable-Top25.csv",
                            sep=""))
  
  #Housekeeping genes
  if(is.null(housekeeping_genes)) { housekeeping_genes <- c("ACTB", "GAPDH") }
  housekeeping_genes <- intersect(housekeeping_genes, row.names(cv_res))
  if(length(housekeeping_genes)>0) {
    data_mat_hg <- cv_res[housekeeping_genes,Sample_group]
    rn_hg <- data.frame(do.call(rbind, strsplit(colnames(data_mat_hg),
                                                split = ":")), stringsAsFactors = FALSE)
    #set.seed(2020)
    ha_col <- HeatmapAnnotation(df=data.frame(group=rn_hg$X1),
                                annotation_name_gp = gpar(fontsize = 6),
                                simple_anno_size = unit(0.3, "cm"))
    ht1_hg <- Heatmap(data.matrix(data_mat_hg), cluster_rows =FALSE,
                      cluster_columns = FALSE,
                      column_split = factor(rn_hg$X1, levels = group_oi),
                      na_col = "grey",
                      col = colorRamp2(c(0,cvThreshold,(cvThreshold+0.001),
                                         (cvThreshold+10)), c("blue","white","pink","red")),
                      row_names_max_width=unit(10, "cm"),
                      column_names_gp = gpar(fontsize = 5),
                      row_names_gp = gpar(fontsize = 6),
                      row_title = "Stable (housekeeping) genes",
                      column_title_gp = gpar(fontsize = 4),
                      top_annotation = ha_col,
                      heatmap_legend_param = list(title = "CV",
                                                  heatmap_legend_side = "right") )
    pdf(paste(filePATH,"/",fileName,"-housekeepingGenes.pdf", sep=""),
        width=10, height=3.5)
    print(ht1_hg)
    dev.off()
  }
  
  #Define the stable genes
  stable_list <- as.character(super_stable2$Var1)
  dfx <- c()
  for(i in 1:length(group_oi)) {
    groupName <- group_oi[i]
    df <- non_variable_gene[non_variable_gene$gene %in% stable_list &
                              non_variable_gene$group %in% groupName,]
    df1 <- data.frame(table(df$gene, df$group))
    df1 <- df1[df1$Freq >= donorThreshold &
                 order(df1$Freq, decreasing=TRUE),]
    df <- df[df$gene %in% df1$Var1,]
    dfx <- rbind(dfx, df[,])
  }
  dfx <- dfx[!is.na(dfx$mean),]
  stable_gene <- unique(dfx$gene)
  write.csv(dfx, file=paste(filePATH,"/",fileName,"-stable-genelist.csv",
                            sep=""))
  #plot heatmap
  data_mat <- stable_mat[stable_gene,Sample_group]
  m <- match(gene_list$gene,row.names(data_mat))
  m <- m[!is.na(m)]
  data_mat <- data_mat[m,]
  print(head(data_mat))
  ht2 <- Heatmap(data.matrix(data_mat), cluster_rows =FALSE,
                 cluster_columns = FALSE,
                 column_split = factor(rn$X1, levels = group_oi),
                 na_col = "grey",
                 col = colorRamp2(c(0,cvThreshold,max(data_mat[!is.na(data_mat)])), c("#4A92E2", "white","#D32B2B")),
                 row_names_max_width=unit(10, "cm"),
                 column_names_gp = gpar(fontsize = 5),
                 row_names_gp = gpar(fontsize = 5),
                 row_title = paste("Stable genes:", length(stable_gene)),
                 column_title_gp = gpar(fontsize = 4),
                 top_annotation = ha_col,
                 heatmap_legend_param = list(title = "CV",
                                             heatmap_legend_side = "right") )
  
  pdf(paste(filePATH,"/",fileName,"-stable-Features.pdf", sep=""),
      width=10, height=10)
  print(ht2)
  dev.off()
  print(ht2)
  write.csv(data_mat, paste(filePATH,"/",fileName,"-stable-Features",".csv", sep=""))
  
  ## Add the result
  data_object@result$stable_genes <- dfx
  message(date(),": Check output directory for results")
  
  return(data_object)
  
}

###############

CVD_list <- fread('CVD_gene.list',data.table=F)

#####bulk protein#####
load("data/AIFI-Olink-NPX_log2_Protein.Rda")
load("data/AIFI-Metadata.Rda")
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

