####To load in relevant packages:---------------
library(GEOquery)
library(limma)

####To retrieve the GSE as an ExpressionSet:---------------
gset <- getGEO("GSE52763", GSEMatrix =TRUE, getGPL=FALSE)
if (length(gset) > 1) idx <- grep("GSE52763", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]
  #this was obtained from the GEO2R R-script
  #It essentially retrieves the GSEdataset as an ExpressionSet

#To see a summary of how many rows and columns there are:
dim(gset)

#To see the names of the GSM samples: use table(<data name>$<column name>)
table(gset$description)
head(gset)
pData(gset) [,1:3]

####To subset the data to ones we're interested in:---------------
sample <- gset[, 1:7]
  #This will make a variable, that consists of ALL the features (genes), and the first 7 samples (animals)
dim(sample)
  #This will let us see the number of rows/columns in the newly created variable
featureNames(sample)[1:10]
  #This will give the names of the features (in this case, the probe names)
sampleNames(sample)[1:10]
  #This will give us the sample names within our subset of data! 
varLabels(sample)
  #This will yield the different variables/columns that each sample has

####To create the design matrix for limma:---------------
design <- model.matrix(~ 0+factor(c(1,1,1,2,2,2,2)))
  #this essentially tells limma which "group" each sample belongs to
  #ex: 1,1,1,2,2,2,2 tells limma that the first 3 GSMs are from RNA1, and the last 4 GSMs are from RNA2
colnames(design) <- c("sham", "SCI")
  #this renames the groups in the design matrix to whatever you want
  #colnames(design) <- c("group1", "group2", "group3")
fit <- lmFit(sample, design)

####To make pairwise comparisons between your groups:
contrast.matrix <- makeContrasts(sham-SCI, levels=design)
  #contrast.matrix <- makeContrasts(group2-group1, group3-group2, group3-group1, levels=design)
  
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
topTable(fit2, coef=1, adjust="BH")
  #this will generate a table of the top DEGs between your comparison groups!

####To compare TWO GROUPS: ex, 3 sham vs 4 SCI mice
design <- cbind(Sham=c(1,1,1,0,0,0,0),SCI=c(0,0,0,1,1,1,1))
fit <- lmFit(sample, design)
cont.matrix <- makeContrasts(Sham-SCI, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
topTable(fit2, adjust="BH")
