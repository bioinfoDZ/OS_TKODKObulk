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





#give it a label...
md$Sample
code <- data.frame(Sample = c("TL724M0T", "TL30F1T",  "TL744F2T", "DL761M2R", "DL579FOR", "DL608F1R"),
                   Code = c('TKO1', 'TKO2', 'TKO3', 'DKO1', 'DKO2', 'DKO3'))

md$Code <- code$Code
colnames(gem) <- code$Code
md$Sample <- code$Code


### prep for mMCP-Counter ###

# input is normalized, log-transformed counts


#normalize
sizefactors <- DESeq2::estimateSizeFactorsForMatrix(gem)
gemnorm <- as.data.frame(t( t(gem) / sizefactors ))

#log transform
gemlog <- log2(gemnorm + 1)
rm(gemnorm)


#estimate deconvolved results...
dc <- mMCPcounter.estimate(gemlog, features = "ENSEMBL.ID")


#rename monocytes/macs to just macs?
rownames(dc)[grepl('Monocytes / macrophages', rownames(dc))] <- 'Macrophages'




#try to plot...
#transpose, samples = rows
dct <- as.data.frame(t(dc))

#merge with metadata
mdx <- cbind(md, dct)

#plot
celltypes <- colnames(dct)
ctlist <- list()
for(ct in celltypes){
  
  pdf <- data.frame(Sample = mdx$Sample,
                    Condition = mdx$Condition,
                    estimate = mdx[,ct],
                    celltype = ct
  )
  
  pdf[pdf$mMCP_estimate==-Inf,"estimate"] <- 0
  
  p <- t.test(pdf[pdf$Condition=='TKO',"estimate"],
              pdf[pdf$Condition=='DKO',"estimate"])$p.value
  
  
  pdf$pval <- p
  
  pdf$meandiff <- mean( pdf[pdf$Condition=='TKO',"estimate"] ) - mean(pdf[pdf$Condition=='DKO',"estimate"])
  
  if(p >= 0.05){
    
    p <- signif(p, digits = 2)
    
    pdf$celltype_pvalue <- paste0(ct, '\nP = ', p)
    
  }
  
  if(p < 0.05 & p >= 0.01) {
    
    p <- signif(p, digits = 2)
    
    pdf$celltype_pvalue <- paste0(ct, '\n* P = ', p, ' *')
    
  }
  
  if(p < 0.01 & p >= 0.001) {
    
    p <- signif(p, digits = 2)
    
    pdf$celltype_pvalue <- paste0(ct, '\n** P = ', p, ' **')
    
  }
  
  if(p < 0.001) {
    
    p <- signif(p, digits = 2)
    
    pdf$celltype_pvalue <- paste0(ct, '\n*** P = ', p, ' ***')
    
  }
  
  
  ctlist[[ct]] <- pdf
  
  
}
ctres <- dplyr::bind_rows(ctlist)

#plot them

#sort celltypes by global mean
ctag <- aggregate(estimate ~ celltype, ctres, mean)
ctres$celltype <- factor(ctres$celltype, levels =   ctag[order(ctag$estimate, decreasing = T),"celltype"])


#sort celltypes pval by lowest pvalue...
ctag <- aggregate(estimate ~ celltype_pvalue, ctres, mean)
ctpval <- data.frame(celltype_pvalue = ctres$celltype_pvalue, 
                     pval = ctres$pval)
ctpval <- ctpval[!duplicated(ctpval$celltype_pvalue),]
ctpval <- ctpval[order(ctpval$pval),]
ctres$celltype_pvalue <- factor(ctres$celltype_pvalue, levels =   ctpval$celltype_pvalue)



#rename estimate column based on tool
colnames(ctres)[3] <- 'mMCP_estimate'


#boxplots comparing conditions

gbox <- ggplot(ctres, aes(Condition, mMCP_estimate, fill = Condition))+
  geom_point()+
  geom_boxplot()+
  facet_wrap(~ celltype_pvalue)+
  scale_fill_manual(values = rev(unique(md$Color)))+
  labs(caption = 'Statistics via T test')+
  theme_linedraw() +
  theme(panel.grid = element_line(color = rgb(235, 235, 235, 195, maxColorValue = 255)))




#sample barplot

#better color scale...
pal <- c( RColorBrewer::brewer.pal(Inf, "Set1"),
          RColorBrewer::brewer.pal(Inf, "Set2")
)

pal <- sample(pal, size = length(levels(ctres$celltype)), replace = F)




gbar <- ggplot(ctres, aes(Sample, mMCP_estimate, fill = celltype))+
  geom_bar(position="stack", stat="identity")+  
  theme_linedraw()+
  theme(axis.text.x = element_text(angle = 45, vjust =1, hjust=1))+
  scale_fill_manual(values = pal)+
  theme(panel.grid = element_line(color = rgb(235, 235, 235, 195, maxColorValue = 255)))





pdf('~/Desktop/AACR_updatedeconvplots_gbox.pdf', height = 5, width = 6)
gbox
dev.off()

pdf('~/Desktop/AACR_updatedeconvplots_gbar.pdf', height = 4, width = 3)
gbar
gbar + Seurat::NoLegend()
dev.off()

pdf('~/Desktop/AACR_updatedeconvplots_gbar_BIGLEGEND.pdf', height = 5, width = 5)
gbar
gbar + Seurat::NoLegend()
dev.off()

# save res
ggsave('results/comparative-de/TKO-vs-DKO/4.downstream/mMCP_counter/barplot.jpg',
       gbar, width = 5, height = 5, dpi = 400)

ggsave('results/comparative-de/TKO-vs-DKO/4.downstream/mMCP_counter/boxplots.jpg',
       gbox, width = 7, height = 7, dpi = 500)


pdf('results/comparative-de/TKO-vs-DKO/4.downstream/mMCP_counter/boxplots.pdf',
    width = 5, height = 5)
gbar

dev.off()



# exclude character formatted celltype + pval...
ctressave <- ctres[,-ncol(ctres)]
write.csv(ctressave,
          'results/comparative-de/TKO-vs-DKO/4.downstream/mMCP_counter/mmcpres.csv', 
          quote = F, row.names = F)






