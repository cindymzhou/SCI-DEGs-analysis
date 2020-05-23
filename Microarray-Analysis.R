####To load in relevant packages:---------------
library(GEOquery)
library(limma)
library(annotate)
library(rat2302.db)
  #here, you have to load the appropriate db file for your organism + array of interest! 
  #files can be obtained at: http://www.bioconductor.org/packages/release/data/annotation/

####To retrieve the GSE as an ExpressionSet:---------------
gset <- getGEO("GSE45006", GSEMatrix =TRUE, getGPL=FALSE)
if (length(gset) > 1) idx <- grep("GSE45006", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]
  #this was obtained from the GEO2R R-script
  #It essentially retrieves the GSEdataset as an ExpressionSet

###To see a preview of your GSE data:
  ##Summary of how many rows and columns there are:
dim(gset)

  ##To see the names of the GSM samples: use table(<data name>$<column name>)
table(gset$description)
head(gset)
pData(gset) [,1:3]
length(sampleNames(gset))


##To convert the probe name to GENE names and IDS
  ##Get the transcript cluster IDs from the expressionset
ID <- featureNames(gset)

# Look up the Gene Symbol, name, and Ensembl Gene ID for each of those IDs
Symbol <- getSYMBOL(ID, "rat2302.db")
Name <- as.character(lookUp(ID, "rat2302.db", "GENENAME"))
Ensembl <- as.character(lookUp(ID, "rat2302.db", "ENSEMBL"))

# Make a temporary data frame with all those identifiers
tmp <- data.frame(ID=ID, Symbol=Symbol, Name=Name, Ensembl=Ensembl, stringsAsFactors=F)

# The stringsAsFactors makes "NA" characters. This fixes that problem.
tmp[tmp=="NA"] <- NA

# set the featureData for your expressionset using the data frame you created above.
fData(gset) <- tmp

# Clean up
rm(ID, Symbol, Name, Ensembl, tmp)
  ##Annotating code is courtesy of: https://www.r-bloggers.com/annotating-limma-results-with-gene-names-for-affy-microarrays/

###To subset the data, and to create the design matrix for limma for MULTIPLE TIME POINT COMPARISONS:---------------
# Identifying which samples we want, code from GEO2R:
gsms <- "XXXXCCCCXXXXXXXX11113333"
    #X represents samples we don't want, numbers/letters represent samples we DO want 
sml <- c()
for (i in 1:nchar(gsms)) { sml[i] <- substr(gsms,i,i) }

# eliminate samples marked as "X"
sel <- which(sml != "X")
sml <- sml[sel]
gset <- gset[ ,sel]

###To set names for the individual groups + run the analysis
sml <- paste("G", sml, sep="")  
fl <- as.factor(sml)
design <- model.matrix(~ fl + 0, gset)
  ####design <- model.matrix(~0+f)
    ##0 tells limma that there's no intercept, and +f is the factor, in our case, it's fl
colnames(design) <- levels(fl)
  ##this renames the design column names to our factor names (aka, it's renamed to the diff groups we have)
fit <- lmFit(gset, design)
cont.matrix <- makeContrasts(ControlvsDay1=GC-G1,
                             Day1vsDay3=G1-G3,
                             ControlvsDay3=GC-G3, levels=design)
  #contrast groups must be separated by "-", whereas the name for the contrasts can be given in front of the = 
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
tT <- topTable(fit2, coef=1, adjust="BH", number=20)
  ##the coefficient specifies which comparison you want to look at. Coef=1 refers to "ControlvsDay1", because it is the first contrast


##To get the FULL table of DEGs and a summary of which are up/down reg:
complete_list <- topTable(fit2, coef=1, adjust="BH", number=Inf)
summary(decideTests(fit2))
