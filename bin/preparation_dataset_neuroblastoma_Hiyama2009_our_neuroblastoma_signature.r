setwd(".")
options(stringsAsFactors = FALSE)
cat("\014")
# set.seed(11)


# For the printed files
num_to_return <- 1
exe_num <- sample(1:as.numeric(10000), num_to_return)

globalYMin <- -1
globalYMax <- 1

# plotTrends()
plotTrends <- function(geneSymbol, thisMeltTable, savePlotToFile) {

    library("ggplot2")
    
    thisTitle <- paste0(geneSymbol, " gene expression: survived versus deceased patients")
    
    numNAs <- thisMeltTable %>% is.na() %>% sum()
    # cat("numNAs =", numNAs, "\n", sep="")
    
    cat("globalYMin = ", globalYMin, "\t")
    cat("globalYMax = ", globalYMax, "\n")
    
    ggplot(thisMeltTable[thisMeltTable$"variable"==geneSymbol,], aes(x = num, y = value)) +
    geom_line(
        data = (thisMeltTable[(thisMeltTable$"variable"==geneSymbol & thisMeltTable$"survival"=="survived"),]),
        mapping = aes(x = num, y = value, colour = survival),
        size = 1,
        alpha = 5 / 6
        ) +
    geom_line(
    data = (thisMeltTable[(thisMeltTable$"variable"==geneSymbol & thisMeltTable$"survival"=="deceased"),]),
        mapping = aes(x = num, y = value, colour = survival),
        size = 1,
        alpha = 5 /6
    ) + labs(title=thisTitle, x="patients number", y="gene expression") + ylim(globalYMin, globalYMax)
    
    
    if(savePlotToFile == TRUE)  {
        plotFileName <- paste0("../results/plot_", geneSymbol, "_rand", exe_num, ".pdf")
        ggsave(plotFileName, width=8, height=4)
        cat("Saved file", plotFileName, "\n", sep="")
    }
  
}

list.of.packages <- c("easypackages", "pacman", "dplyr", "reshape") # other packages
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos="https://utstat.toronto.edu/cran/")

library("pacman")
pacman::p_load("dplyr")

source("utils.r")

verbose <- TRUE



if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
listOfBiocPackages <- c("annotate", "GEOquery", "limma")

bioCpackagesNotInstalled <- which( !listOfBiocPackages %in% rownames(installed.packages()) )
cat("package missing listOfBiocPackages[", bioCpackagesNotInstalled, "]: ", listOfBiocPackages[bioCpackagesNotInstalled], "\n", sep="")

# check there's still something left to install
if( length(bioCpackagesNotInstalled) ) {
    BiocManager::install(listOfBiocPackages[bioCpackagesNotInstalled])
}

library("easypackages")
libraries(list.of.packages)
libraries(listOfBiocPackages)

GSE_code <- "GSE16237" # this line will change for each different dataset
thisGEOplatform <- "GPL570" # this line will change for each different dataset
datasetName <-  "Hiyama2010" # this line will change for each different dataset
cancer_type <- "neuroblastoma" # this line will change for each different dataset

cat("\n\tGSE_code: ", GSE_code, "\n", sep="")
cat("\tthisGEOplatform: ", thisGEOplatform, "\n", sep="")
cat("\tdatasetName: ", datasetName, "\n", sep="")
cat("\tcancer_type: ", cancer_type, "\n", sep="")

patients_data <- NULL

gset <- getGEO(GSE_code,  GSEMatrix =TRUE, getGPL=FALSE)

if (length(gset) > 1) idx <- grep(thisGEOplatform, attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]
gene_expression <- as.data.frame(exprs(gset))

cat("str(gset@phenoData@data)\n")
print(str(gset@phenoData@data))
cat("str(gset@phenoData@data)\n")

LABEL_DETECTED <- TRUE
 
 if(LABEL_DETECTED == TRUE) {
 
    # # # # we add the labels
    library("plyr")
    label_list <- c()
    i <- 1
    for(thisTitle in gset@phenoData@data$"outcome of the patient:ch1") { # this line will change for each different dataset
      
	if(grepl("Died of disease", thisTitle)) { # this line will change for each different dataset
	      label_list[i] <-  0
	 } else if(grepl("Alive", thisTitle)) { # this line will change for each different dataset
	    label_list[i] <- 1
	  }
	    i <- i + 1
  }
      
    cat("label_list:\n")
    print(label_list)
        
    gset_expression_colmeans_mean <- gene_expression %>%  colMeans() %>% mean()
    gset_expression_colmeans_sd <-gene_expression %>%  colMeans() %>% sd()
        
    cat("mean of the average gene expression per patient profile +- standard deviation = ", gset_expression_colmeans_mean, " +- ", gset_expression_colmeans_sd,"\n", sep="")
     
     batchCorrection <- TRUE
     if(batchCorrection == TRUE) {
     
        cat("doing batch correction through limma::removeBatchEffect()\n")
        output_batch_correction <- limma::removeBatchEffect(gene_expression, label_list)
     
        output_batch_correction_mean <- output_batch_correction %>%  colMeans() %>% mean()
        output_batch_correction_sd <-output_batch_correction %>%  colMeans() %>% sd()
        
        cat("mean of the average gene expression per patient profile after batch correction +- standard deviation = ", output_batch_correction_mean, " +- ", output_batch_correction_sd,"\n", sep="")
        
        gene_expression <- output_batch_correction
        
    }
    
    targetName <- "survival"
    
    labels_df_temp <- as.data.frame(label_list)
    
    labels_df <- as.data.frame(t(labels_df_temp))
    colnames(labels_df) <- colnames(gene_expression)
    rownames(labels_df) <- targetName
    gene_expression_with_labels <- rbind(labels_df, gene_expression)
    
    gene_expression_with_labels$ID <- rownames(gene_expression_with_labels)

    AHCY_DPYSL3_NME1_signature <- c("200903_s_at", "201431_s_at", "201577_at" )
    
    cat("Here we use the AHCY-DPYSL3-NME1 signature\n")
    print(AHCY_DPYSL3_NME1_signature)
    cat("\n")
    patients_data_filtered_ourSignature <- (gene_expression_with_labels[gene_expression_with_labels$ID %in% c(AHCY_DPYSL3_NME1_signature, "survival"), ])
    
    patients_data_filtered_ourSignature$ID <- NULL
    tableZerosAndOnes <- as.data.frame(t(patients_data_filtered_ourSignature))
    
    tableZerosAndOnes <- tableZerosAndOnes %>% relocate(survival, .after = last_col())
    
    colnames(tableZerosAndOnes) <- c("AHCY", "DPYSL3", "NME1", "survival")

}

tableZerosAndOnes$num <- seq(1:(tableZerosAndOnes %>% nrow()))
tableZerosAndOnes$patientID <- rownames(tableZerosAndOnes)

tableZerosAndOnes$"survival"  <-  tableZerosAndOnes$"survival" %>% as.integer()

meltTable <- melt(tableZerosAndOnes, id=c("patientID", "survival", "num"))

meltTable$"survivalBoolean" <-  meltTable$"survival" %>% as.logical()
meltTable[meltTable$"survivalBoolean"==TRUE,]$survival <- "survived"
meltTable[meltTable$"survivalBoolean"==FALSE,]$survival <- "deceased"

ratioForRange <- 0.1
minValue <- ((meltTable$value %>% summary())[1] %>% as.numeric())
maxValue <- ((meltTable$value %>% summary())[6] %>% as.numeric())
globalYMin <- minValue - 1100
globalYMax <- maxValue + (maxValue*ratioForRange)

# p_deceased <- ggplot(meltTable[(meltTable$variable=="AHCY" & meltTable$survival==0),], aes(ID, value, colour=variable)) + geom_point() + geom_line()
# p_alive <- ggplot(meltTable[(meltTable$variable=="AHCY" & meltTable$survival==1),], aes(ID, value, colour=variable)) + geom_point() + geom_line()

plotTrends("AHCY", meltTable,  TRUE)
plotTrends("DPYSL3", meltTable,  TRUE)
plotTrends("NME1", meltTable,  TRUE)
