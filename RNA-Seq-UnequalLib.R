###Load in libraries:
library(edgeR)
library(limma)
library(biomaRt)
library(tidyverse)

#First, import the count file, and make sure to create clear row and column names
counts <- read.csv("GSE93249_raw_counts.csv")
  #with working directory set as "Documents"
##To isolate the coding RNA samples:
counts <- counts[,-(2:3)]
counts <- data.frame(counts[,-1], row.names = counts[,1])
  #This helps make the first column into the row name instead
counts <- counts[-(22076:30443),]
  ##GSE93249 has long non-coding RNA data, this line deletes that and only keeps the coding transcripts


##Now, we use edgeR's DGElist function
dge <- DGEList(counts=counts)


##Making the design matrix manually for time-course experiments:
design <- model.matrix(~ 0+factor(c(1,1,1,2,2,2,3,3,3,6,6)))
colnames(design) <- c("control", "m1", "m3", "m6")

##Filtering out counts that are close to 0
keep <- filterByExpr(y = dge$counts, design = design)
dge <- dge[keep,,keep.lib.sizes=FALSE]

#TMM Normalizing
dge <- calcNormFactors(dge)

#For samples with UNEQUAL sequencing depth (largest library size:smallest lib size = > than 3), we use voom:
v <- voom(dge, design, plot=TRUE)
fit <- lmFit(v, design)
  ##Then, we can run the limma pipeline as usual! 
cont.matrix <- makeContrasts(Controlvs1Month = control-m1, 
                             Controlvs3Months= control-m3,
                             Controlvs6Months= control-m6,
                             levels=design)
  ##This contrast matrix asks "Compared to control, which genes are responding at either 1 month, 3 months, OR 6 months?"
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit)
topTable(fit2, coef=1, adjust="BH")
  ##the coef specifies which contrast/comparison you want to look at. Coef=1 here refers to "Controlvs1Month"
  ##coef=ncol(design) makes it so that the coefficients (aka, # of RNA samples) is equivalent to the number of rows/columns in the design matrix 
complete_list <- topTable(fit2, coef=1, number=Inf, adjust="BH")


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
gene <- getBM(attributes = c("ensembl_gene_id","external_gene_name"),values = complete_list$ID,mart = ensembl)  
#Now we need to match the gene id with ensembl_gene_id
id <- match(complete_list$ID, gene$ensembl_gene_id)
#Add Gene symbol column into the data frame
complete_list$ID <- gene$external_gene_name[id]
head(complete_list)
  ##Converting RNA-Seq names to gene names is courtesy of: https://www.biostars.org/p/337072/ 

##To see a summary of the number of differentially expressed genes:
summary(decideTests(fit2))
