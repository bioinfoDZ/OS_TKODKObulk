library(DESeq2)
library(tidyverse)
library(ComplexHeatmap)
library(FerrenaBulkRNAseq)

library(mMCPcounter)

set.seed(2021)


#read in
gem <- readRDS('data/gem.rds')
gencode <- read.csv('data/gencode.vM23.ids_names_types.csv')
md <- readxl::read_excel('data/metadata.xlsx')



#keep only protein coding
gencode <- gencode[gencode$gene_type == 'protein_coding',]
gem <- gem[match(gencode$gene_id, rownames(gem)),]


# for easier time, remove id version info
gencode$id_with_version <- gencode$gene_id
split <- str_split_fixed(string = gencode$gene_id, '\\.', 2)
gencode$gene_id <- split[,1] ; rm(split)

rownames(gem) <- gencode$gene_id





# quantile normalization
# code taken from their github
# https://github.com/wuaipinglab/ImmuCC/blob/master/webserver/MouseHTSeq_counts_stat.R#L73
#data = gene exp matrix, counts, Htseq output
# p = quantile, default 0.75
data <- gem
p <- 0.75

n <- ncol(data)
result <- matrix(nrow=nrow(data), ncol=0)
for (i in seq(n)) {
  expres <- as.numeric(data[, i])
  value <- expres[expres != 0]
  value <- sort(value, decreasing=T)
  value.quantile <- quantile(value, p)
  scale <- as.numeric(value.quantile)/1000
  expres <- ceiling(expres/scale)
  result <- cbind(result, expres)
}
colnames(result) <- colnames(data)
rownames(result) <- rownames(data)


#get quant-normed result
gemnorm <- result




#write out result


write.csv(gemnorm, 'data/bulk_seqimmunecc.csv',
            row.names = T, quote = F)





############ read in results #####################
dc <- read.csv('data/immunecc/result.SVR.csv', row.names = 1)

#try to plot...


#stackedbaplot

#merge with metadata
mdx <- cbind(md, dc)

#plot
celltypes <- colnames(dc)
ctlist <- list()
for(ct in celltypes){
  
  pdf <- data.frame(Sample = mdx$Sample,
                    Condition = mdx$Condition,
                    estimate = mdx[,ct],
                    celltype = ct
  )
  
  
  p <- t.test(pdf[pdf$Condition=='TKO',"estimate"],
              pdf[pdf$Condition=='DKO',"estimate"])$p.value
  
  p <- signif(p, digits = 2)
  
  pdf$celltype_pvalue <- paste0(ct, '\nP = ', p)
  
  ctlist[[ct]] <- pdf
  
  
}
ctres <- dplyr::bind_rows(ctlist)

#plot them

#sort celltypes by global mean?
ctag <- aggregate(estimate ~ celltype, ctres, mean)
ctres$celltype <- factor(ctres$celltype, levels =   ctag[order(ctag$estimate, decreasing = T),"celltype"])

ctag <- aggregate(estimate ~ celltype_pvalue, ctres, mean)
ctres$celltype_pvalue <- factor(ctres$celltype_pvalue, levels =   ctag[order(ctag$estimate, decreasing = T),"celltype_pvalue"])


#rename estimate column based on tool
colnames(ctres)[3] <- 'seqImmuneCC_estimate'

ggplot(ctres, aes(Condition, seqImmuneCC_estimate, fill = Condition))+
  geom_point()+
  geom_boxplot()+
  facet_wrap(~ celltype_pvalue)+
  scale_fill_manual(values = rev(unique(md$Color)))+
  labs(caption = 'Statistics via T test')



#sample barplot


#better color scale...
pal <- c( RColorBrewer::brewer.pal(Inf, "Set1"),
          RColorBrewer::brewer.pal(Inf, "Set2")
)

pal <- sample(pal, size = length(levels(ctres$celltype)), replace = F)


ggplot(ctres, aes(Sample, seqImmuneCC_estimate, fill = celltype))+
  geom_bar(position="stack", stat="identity")+
  theme(axis.text.x = element_text(angle = 45, vjust =1, hjust=1))+
  scale_fill_manual(values= pal)


