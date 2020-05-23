###Load in libraries:
library(edgeR)
library(limma)
library(biomaRt)
library(tidyverse)

#First, import the count file, and make sure to create clear row and column names
counts <- read.delim("GSE133093_counts_IH.txt")
  ##with working directory as "Documents"
  ##GEOQuery cannot download RNA-seq count files :(
counts <- data.frame(counts[,-1], row.names = counts[,1])
  #This helps make the first column into the row name instead

##Now, we use edgeR's DGElist function
dge <- DGEList(counts=counts)

##Making the design matrix:
###To create the design matrix for limma for MULTIPLE TIME POINT COMPARISONS:---------------
# group names for all samples
gsms <- "AAAABBBXXXXXXXXXX"
sml <- c()
for (i in 1:nchar(gsms)) { sml[i] <- substr(gsms,i,i) }

# eliminate samples marked as "X"
sel <- which(sml != "X")
sml <- sml[sel]
dge <- dge[ ,sel]

###To set names for the individual groups 
sml <- paste("G", sml, sep="")  
fl <- as.factor(sml)
design <- model.matrix(~fl+0, dge)
colnames(design) <- levels(fl)

##Filtering out counts that are close to 0
keep <- filterByExpr(y = dge$counts, design = design)
dge <- dge[keep,,keep.lib.sizes=FALSE]

#TMM Normalizing
dge <- calcNormFactors(dge)

#For samples with equal sequencing depth (largest library size:smallest lib size = less than 3), we use limma-trend
logCPM <- cpm(dge, log=TRUE, prior.count=3)
  ##First we must convert to logCPM value
fit <- lmFit(logCPM, design)
  ##Then, we can run the limma pipeline as usual! 
cont.matrix <- makeContrasts(ShamvsSCI = GA-GB, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit, trend=TRUE)
  #The trend = TRUE argument is needed when running limma with logCPM values 
topTable(fit2, coef=1)
  ##coef=ncol(design) makes it so that the coefficients (aka, # of RNA samples) is equivalent to the number of rows/columns in the design matrix
  ##QUESTION:....what does making coef=design column number compare....?
complete_list <- topTable(fit2, number=Inf, adjust="BH", coef=1)
  #topTableF is equivalent to topTable(coef=NULL, or not specifying coef at all)
  ##If you do not specify coef, topTable will test whether the genes vary between the RNA samples in ANY way

##To get the gene IDs into the table:
complete_list <- rownames_to_column(complete_list, var="ID")
  #This will make the ensembl rownames into a column! So that it can be name matched 
mart = useMart('ensembl')
  #This defines which database you want to pull IDs from
#To list all the ensembl database of organisms:
listDatasets(mart) 
  ##This will list ALL of the ones available 
searchDatasets(mart = mart, pattern = "norvegicus")
  ##This will help you search for a particular organism, using close words 
#Next, we need to choose database of  interest:
ensembl = useMart( "ensembl", dataset = "rnorvegicus_gene_ensembl" )  
# choose attributes of your interest
listAttributes(ensembl)
  ##this will give all the attributes (i.e. gene ID, transcript ID)
gene <- getBM( attributes = c("ensembl_gene_id","external_gene_name"),values = complete_list$ID,mart = ensembl)  
#Now we need to match the gene id with ensembl_gene_id
id <- match(complete_list$ID, gene$ensembl_gene_id)
#Add Gene symbol column into the data frame
complete_list$ID <- gene$external_gene_name[id]
head(complete_list)
  ##Converting RNA-Seq names to gene names is courtesy of: https://www.biostars.org/p/337072/ 

##To see a summary of the number of differentially expressed genes:
summary(decideTests(fit2))
