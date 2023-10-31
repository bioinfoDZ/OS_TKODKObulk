library(DESeq2)
library(tidyverse)
library(cowplot)
library(MAST)
library(fgsea)
library(extrafont)
library(EnhancedVolcano)
library(readxl)

library(foreach)
library(doParallel)

library("RColorBrewer")
library(pheatmap)
library(ComplexHeatmap)

library(dendextend)
library(dendsort)

library(biomaRt)

set.seed(2020)


#setwd("/Users/ferrenaa/Documents/deyou/target/")


#parse raw files
files <- list.files('data/rawdata/', full.names = T)

#read in one file, get gene names
names <- read.table(files[1], header = T)
names <- names[,1]


#progress bar
total = length(files) - 1
pb <- txtProgressBar(min = 0, max = total, style = 3)

datlist <- list()
for(filedex in 1:length(files)) {
  
  file <- files[filedex]
  
  samp <- basename(tools::file_path_sans_ext(file))
  samp <- strsplit(x=samp, split = '\\.')[[1]][1]
  
  dat <- as.data.frame(read.table(file, header = T)[,4])
  colnames(dat) <- samp
  
  datlist[[filedex]] <- dat
  rm(dat)
  
  setTxtProgressBar(pb, filedex)
  
}

#bind the list to a big df
gem <- dplyr::bind_cols(datlist)
rownames(gem) <- names
rm(datlist, names)

#fiter lowly exp genes
gem <- gem[rowSums(gem) > 10,]


### ensembl ; get gene names, select only coding genes ###


#some really weird gene names...
# the "\\." indicates the version i think...
# underscores?? what does it mean? --> nonstandard... https://www.biostars.org/p/437927/
genes <- rownames(gem)
head(genes, 100)

#remove underscores...
genes <- sapply(str_split(rownames(gem), '_'), FUN = function(x){x[[1]][1]})



#select the protein-coding genes; need to use ensembl / biomart
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# list <- getBM(attributes=c("ensembl_gene_id","ensembl_gene_id_version", "external_gene_name","transcript_biotype"),
#               filters = "ensembl_gene_id_version", values=genes, mart=mart)

#doing this, we get ~ 24-thousand genes out, 2746 protein coding...


#remove the '\\.' suffix and try the stable gene id...
genes <- sapply(str_split(genes, '\\.'), FUN = function(x){x[[1]][1]})

list <- getBM(attributes=c("ensembl_gene_id","external_gene_name","transcript_biotype"),
              filters = "ensembl_gene_id", values=genes, mart=mart)

#removing version number we get ~67 thousand genes, 18641 protein coding

#select non-duplicate coding genes...
coding <- list[list$transcript_biotype == 'protein_coding',]
coding <- coding[!duplicated(coding$external_gene_name),]

rm(list)

#remove non protein coding genes from gem...
rownamesdf <- data.frame(raw = rownames(gem), split = genes)
rownamesdf <- rownamesdf[rownamesdf$split %in% coding$ensembl_gene_id,]

gem <- gem[rownames(gem) %in% rownamesdf$raw,]

#match order of gene names, 
coding <- coding[match(rownamesdf$split, coding$ensembl_gene_id),]

#rename gem rows with the gene symbols
rownames(gem) <- coding$external_gene_name







rm(mart,pb,rownamesdf,genes,files,samp,total,file,filedex)



### try to pair with clinical data...


#sample ID format:
# https://ocg.cancer.gov/sites/default/files/OCG%20Sample%20Codes_finalized%2009%202020.pdf

#TARGET-##-ABC123-##A-.3X-01R

#TARGET - tumorcode (40 = OS) - case ID - tissue code, aliquot - .3 (cell culture?), tissue portion? - Assay, ie DNA/RNA


#TARGET-40-0A4HLD-01A-01R


#get samps
samps <- colnames(gem)
sp <- str_split_fixed(samps,'-',Inf)

#remove samps with 06R code... no definition for nucleic acid type...
table(sp[,5])
samps <- samps[sp[,5]!='06R']
sp <- str_split_fixed(samps,'-',Inf)
gem <- gem[,colnames(gem) %in% samps]


#get caseIDs
cases <- sp[,1:3]
cases <- apply(cases, 1, function(x){ paste(x, collapse = '-') } )


#keep as df
sdf <- data.frame(samps = samps,
                  cases = cases)
rm(samps, cases, sp)

#read in files
clin <- as.data.frame(readxl::read_excel('data/clinical/TARGET_OS_ClinicalData_Discovery_20181009.xlsx'))
sampmat <- as.data.frame(readxl::read_excel('data/clinical/TARGET_OS_SampleMatrix_Discovery_20190808.xlsx'))

#keep the intersect of the three things; 
# the gem colnames (samples); 
# the clinical data excel
# and the samplematrix file...

sdf <- sdf[sdf$cases %in% clin$`TARGET USI`,]
clin <- clin[clin$`TARGET USI` %in% sdf$cases,]

sampmat <- sampmat[sampmat$`Case USI` %in% sdf$cases,]

gem <- gem[,colnames(gem) %in% sdf$samps]

#order properly...
clin <- clin[match(sdf$cases, clin$`TARGET USI`),]
sampmat <- sampmat[match(sdf$cases, sampmat$`Case USI`),]




#for saving; make a list object
l <- list('gem' = gem,
          'samp_case_ids' = sdf,
          'clinical_data' = clin,
          'metadata' = sampmat,
          'ensembl_hugo' = coding)


#for cibersortx, save TPM file
gemx <- as.data.frame(gem)
gemx <- cbind(rownames(gemx), gemx)
colnames(gemx)[1] <- 'GeneSymbol'

write.table(gemx, 'data/parsed/targetgem.tsv',
            row.names = F, quote = F, sep = '\t')


beepr::beep()
# rm(gem,sdf,clin,sampmat,coding)

#save the gem and the coding df (with gene names and ensembl IDs)
saveRDS(l, 'data/parsed/parsed_tpm.rds')

#rm(l)







