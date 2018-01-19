## Differential gene expression analysis on TCGA
# esfahani@aices.rwth-aachen.de, 19.1.2018

######################## Data Preparation:
allVoomedFiles <- list.files(path = "/Volumes/compbio/compbio_share/data/TCGA-GDAC(OldVersion_ForMoreInfoAskAli)/",
                             full.names = T)
allVoomedFiles <- allVoomedFiles[-grep("clinical",allVoomedFiles)]

TCGABigMat <- matrix(nrow = 20501,ncol = 0)
NumOfSamples <- 0
for(voomFile in allVoomedFiles){
  load(voomFile)
  # colnames(expr) <- paste( colnames(expr), voomFile)
  TCGABigMat <- cbind(TCGABigMat, expr)
  NumOfSamples <- NumOfSamples + ncol(expr)
  print(paste(voomFile,"doNE!!!"))
}

##Removing the sample duplications (There are super-series in TCGA containing exact same samples as sub-series):
TCGABigMat <- TCGABigMat[,!duplicated(colnames(TCGABigMat))]

## didn't convert from Symbol to Entrez cause of loss in conversion while you actually don't need Entrez.

##Making the annotations:
annotBarcodes <- colnames(TCGABigMat)
sampleTypeBarcode <- sapply(X = strsplit(annotBarcodes,split = "-"),
                            function(x) sapply(X = strsplit(x[4],split = ""),
                                               function(n) paste0(n[1],n[2])))

sampleTypeBarcode <- as.integer(sampleTypeBarcode)
sampleTypeConv <- read.csv(file = "./sampleType.csv")
sampleType <- factor(sapply(X = sampleTypeBarcode,function(x) sampleTypeConv$Definition[x==sampleTypeConv$Code]))

tissueSourceSiteBarcode <- sapply(X = strsplit(annotBarcodes,split = "-"),function(x) x[2])
tissueSourceSiteConv <- read.csv(file = "./tissueSourceSite.csv")
# "NA" is interpreted as NA -> messy correction:
tmpTssCode <- as.character(tissueSourceSiteConv$TSS.Code)
tmpTssCode[523] <- "NA"
tissueSourceSiteConv$TSS.Code <- as.factor(tmpTssCode)
#
tissueSourceSite <- factor(sapply(X = tissueSourceSiteBarcode,
                                  function(x) tissueSourceSiteConv$Study.Name[x==tissueSourceSiteConv$TSS.Code]))

sampleTissueAnnotation <- as.factor(paste(sampleType,tissueSourceSite,sep = ", "))

TCGAAnnot <- data.frame(sampleType,tissueSourceSite,sampleTissueAnnotation)
########################
######################## Data Normalization:
TCGABigMatRelative <- TCGABigMat - apply(X = TCGABigMat,MARGIN = 1,FUN = mean) # Centering genes as a basic normalization (Look at the hist of samples before and after and you'll see)
diseaseLevels <- as.factor(TCGAAnnot$sampleTissueAnnotation)
#Some clean-up:
rm(list = setdiff(x = ls(), y = c("TCGABigMatRelative","diseaseLevels")))
gc()
########################
######################## Hard coded annots:(Contains all cancers that have both Primary solid tumor and Solid tissue healthy, other choices are moglich)
grup1 <- c("Solid Tissue Normal, Bladder Urothelial Carcinoma", "Solid Tissue Normal, Breast invasive carcinoma",
           "Solid Tissue Normal, Cervical squamous cell carcinoma and endocervical adenocarcinoma",
           "Solid Tissue Normal, Cholangiocarcinoma", "Solid Tissue Normal, Esophageal carcinoma ",
           "Solid Tissue Normal, Glioblastoma multiforme","Solid Tissue Normal, Head and Neck squamous cell carcinoma",
           "Solid Tissue Normal, Kidney Chromophobe","Solid Tissue Normal, Kidney renal clear cell carcinoma",
           "Solid Tissue Normal, Kidney renal papillary cell carcinoma","Solid Tissue Normal, Liver hepatocellular carcinoma",
           "Solid Tissue Normal, Lung adenocarcinoma","Solid Tissue Normal, Lung squamous cell carcinoma",
           "Solid Tissue Normal, Pancreatic adenocarcinoma","Solid Tissue Normal, Pheochromocytoma and Paraganglioma",
           "Solid Tissue Normal, Prostate adenocarcinoma","Solid Tissue Normal, Sarcoma", "Solid Tissue Normal, Skin Cutaneous Melanoma",
           "Solid Tissue Normal, Stomach adenocarcinoma","Solid Tissue Normal, Thyroid carcinoma","Solid Tissue Normal, Thymoma",
           "Solid Tissue Normal, Uterine Corpus Endometrial Carcinoma")
grup2 <- c("Primary solid Tumor, Bladder Urothelial Carcinoma", "Primary solid Tumor, Breast invasive carcinoma",
           "Primary solid Tumor, Cervical squamous cell carcinoma and endocervical adenocarcinoma",
           "Primary solid Tumor, Cholangiocarcinoma", "Primary solid Tumor, Esophageal carcinoma ",
           "Primary solid Tumor, Glioblastoma multiforme","Primary solid Tumor, Head and Neck squamous cell carcinoma",
           "Primary solid Tumor, Kidney Chromophobe","Primary solid Tumor, Kidney renal clear cell carcinoma",
           "Primary solid Tumor, Kidney renal papillary cell carcinoma","Primary solid Tumor, Liver hepatocellular carcinoma",
           "Primary solid Tumor, Lung adenocarcinoma","Primary solid Tumor, Lung squamous cell carcinoma",
           "Primary solid Tumor, Pancreatic adenocarcinoma","Primary solid Tumor, Pheochromocytoma and Paraganglioma",
           "Primary solid Tumor, Prostate adenocarcinoma","Primary solid Tumor, Sarcoma", "Primary solid Tumor, Skin Cutaneous Melanoma",
           "Primary solid Tumor, Stomach adenocarcinoma","Primary solid Tumor, Thyroid carcinoma","Primary solid Tumor, Thymoma",
           "Primary solid Tumor, Uterine Corpus Endometrial Carcinoma")
########################
######################## Limma fitting:
require(limma)
TCGAAbbrevConv <-
  read.csv("./TCGANamesAbbrevConvTable.csv",
           sep = ";",
           header = F)
for (K in 1:length(grup1)) {
  GLUPH <- diseaseLevels == grup1[K]
  GLUPC <- diseaseLevels == grup2[K]
  dEannots <-
    as.factor(c(rep("Normal", sum(GLUPH)), rep("Tumor", sum(GLUPC))))
  subMat <-
    cbind(TCGABigMatRelative[, GLUPH], TCGABigMatRelative[, GLUPC])
  
  design <- model.matrix( ~ 0 + as.factor(dEannots))
  colnames(design) <- c("Normal", "Tumor")
  
  fit <- lmFit(subMat, design)
  fitC <-
    contrasts.fit(fit, makeContrasts(Tumor - Normal, levels = design))
  
  #Used empirical Bayes: (there are other options available here like treat)
  BayesOuti <- eBayes(fitC) 
  output <- topTable(BayesOuti, number = Inf)
  
  #Making the output file name:
  outputName <- TCGAAbbrevConv$V1[match(substr(x = grup1[K], start = 22, stop = nchar(grup1[K])), TCGAAbbrevConv$V2)]
  #Fixing bug (extra space in name), not my mistake! TCGA came with this extra space:
  if(is.na(outputName)) outputName <- TCGAAbbrevConv$V1[match(substr(x = grup1[K], start = 22, stop = nchar(grup1[K])-1), TCGAAbbrevConv$V2)]
  write.csv(file = paste0("Results/",outputName,".csv"), x = output)
}
########################