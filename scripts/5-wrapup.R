library(DESeq2)
library(tidyverse)
library(cowplot)
library(fgsea)
library(readxl)
library("RColorBrewer")
library(ComplexHeatmap)
library(biomaRt)
library(dendsort)
library(msigdbr)
library(FerrenaBulkRNAseq)


set.seed(2021)


### wrap up final analyses


### m1 vs m1, module score
# need TPM


### subtract myeloid from TKO, then check survival...
# need DE res an gsea myeloid
# need to do in target dir


# multivar survival 
# need to do in target dir








### m1 vs m1, module score
# need TPM

deres <- 'results/comparative-de/TKO-vs-DKO/2.deresults/DEresults-normalized_count_matrix.csv'
res <- read.csv(deres)

md <- readxl::read_excel('data/metadata.xlsx')

tpm <- readRDS('data/gem-tpm.rds')


#remove point / version info from ensembl IDs
res$ensembl_gene_id <- str_split_fixed(res$ensembl_gene_id, '\\.', 2)[,1]

rownames(tpm) <- str_split_fixed(rownames(tpm) , '\\.', 2)[,1]



#try mac-spectrum...

# we need TPM matrix, and a single column condition csv
# https://macspectrum.uconn.edu/

head(tpm)

#prep tpm mat
tpmsave <- cbind(rownames(tpm), tpm)

colnames(tpmsave)[1] <- 'Ensembl_ID'


#prep cond df
cond <- data.frame(feature = md$Condition)


#save
macspectrumdir <- 'results/comparative-de/TKO-vs-DKO/4.downstream/macspectrum'

write.csv(tpmsave, paste0(macspectrumdir, '/tpm.csv'),
          quote = F, row.names = F)


write.csv(cond, paste0(macspectrumdir, '/feature.csv'),
          quote = F, row.names = F)




#keep only TKO
# first col is gene ids
tpmsave <- tpmsave[,1:4]
cond <- cond[1:3,]




write.csv(tpmsave, paste0(macspectrumdir, '/tpm_justTKO.csv'),
          quote = F, row.names = F)


write.csv(cond, paste0(macspectrumdir, '/feature_justTKO.csv'),
          quote = F, row.names = F)


rm(tpmsave, cond)

#parse and plot...


mac <- read.csv(paste0(macspectrumdir, '/MPI_AMDI_table.csv'))

macgg <- ggplot(mac, aes(MPI, AMDI, col = Feature, label=Samples))+
  geom_point()+
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0)+
  ggrepel::geom_text_repel()+
  theme_classic()+
  xlab('Macrophage Polarization Index\n(MPI)') +
  ylab('Activation-Induced Differentiaton Index\n(AMDI)')+
  scale_y_continuous(limits = c(-30, 30))+
  scale_x_continuous(limits = c(-10,10))+
  scale_color_brewer(palette = 'Set1', direction = -1)


ggsave(paste0(macspectrumdir, '/macgg.jpg'), macgg, 
       height = 4, width = 5)


mac_onlyTKO <- read.csv(paste0(macspectrumdir, '/onlyTKO_MPI_AMDI_table.csv'))

macgg_onlyTKO <- ggplot(mac_onlyTKO, aes(MPI, AMDI, col = Feature, label=Samples))+
  geom_point()+
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0)+
  ggrepel::geom_text_repel()+
  theme_classic()+
  xlab('Macrophage Polarization Index\n(MPI)') +
  ylab('Activation-Induced Differentiaton Index\n(AMDI)')+
  scale_y_continuous(limits = c(-30, 30))+
  scale_x_continuous(limits = c(-10,10))+
  scale_color_brewer(palette = 'Set1', direction = -1)

ggsave(paste0(macspectrumdir, '/macgg_onlyTKO.jpg'), macgg_onlyTKO, 
       height = 4, width = 5)




### remake whole msigdb GSEA

c1 = "TKO"
c2 = "DKO"
comptitle <- paste0(c1, '-vs-', c2)

#make whole results dir
compoutdir <- paste0('results/comparative-de/', comptitle)
#whole msigdb gsea
cat = 'Whole_MSIGDB'
catdir <-  paste0(compoutdir, '/3.gsea/', cat)


gseares <- read.csv(paste0(catdir, '/fgsea-results.csv'))
gseadotplot <- FerrenaBulkRNAseq::gsea.dotplot.onecol(gseares = gseares,
                                                      pathwayfontsize = 11, ntop = 15) +
  ggtitle(cat)

#save dotplot
ggsave(plot=gseadotplot, 
       filename = '/Users/ferrenaa/Dropbox/Result_from_Alex/2021october-TKOvsDKO_bulk/results/comparative-de/TKO-vs-DKO/3.gsea/Whole_MSIGDB/gseadotplot_zoom.jpg',
       height = 10, width = 10, bg = 'white')













#remake volcanoplot

# yaxis, padj
# also make bigger, height lower...
res <- read.csv('results/comparative-de/TKO-vs-DKO/2.deresults/DEresults-normalized_count_matrix.csv')
lfc_thres=1

gencode <- read.csv('data/gencode.vM23.ids_names_types.csv')

rownames(res) <- res$ensembl_gene_id

volcano <- volcanoplot(results = res, lfc_thres = lfc_thres,
                                          colors = c('firebrick', 'steelblue'),
                                          change_gene_label=T, use_padj=T,
                                          gene_label_equivalency=gencode)



volcano



ggsave('results/comparative-de/TKO-vs-DKO/2.deresults/volcano_bigger.jpg', width = 5, height = 4.5, dpi=300)









##### DE of some key genes related to E2F1 / Apoptosis, qPCR style #####

res <- read.csv('results/comparative-de/TKO-vs-DKO/2.deresults/DEresults-normalized_count_matrix.csv')




# get key genes
keygenes <- c('Skp2','E2f1', 'Bbc3', 'Bid', 'Bcl2l11', 'Casp3', 'Apoe', 'Ctss', 'Sh2d6', 'Lcn2')

# set ref gene
ref <- c('Gapdh')


#make sure they are in the gene list
keygenes %in% res$mgi_symbol
ref %in% res$mgi_symbol




#qpcr style DE:
# ratio of each gene over ref gene, gapdh.
# compare ratios in TKO vs DKO via wilcox...


#get the matrix
smallres <- res[match(c(keygenes,ref), res$mgi_symbol),]
smallmat <- smallres[,c(2,9:ncol(smallres))]
rownames(smallmat) <- smallmat[,1] ; smallmat <- smallmat[,-1]

#get test vs ref
testmat <- smallmat[keygenes,]
refmat <- smallmat[ref,]

#divide each column; do it very explicity and carefully...
outlist <- list()
for(i in 1:ncol(testmat)){
  
  samp <- colnames(testmat)[i]
  testcol <- testmat[,i]
  refscalar <- refmat[,i]
  
  out <- data.frame(testcol / refscalar)
  rownames(out) <- rownames(testmat)
  colnames(out) <- samp
  
  outlist[[i]] <- out
}

ratiomat <- dplyr::bind_cols(outlist)


# ANALYZE:
# wilcox test
c1ratios <- ratiomat[,1:3]
c2ratios <- ratiomat[,4:6]

compres <- lapply(1:length(keygenes), function(i) {
  gene <- keygenes[i]
  
  xvec <- t(c1ratios[i,])[,1]
  yvec <- t(c2ratios[i,])[,1]
  
  wilcoxout <- wilcox.test(xvec , yvec)
  tout <- t.test( xvec,  yvec )
  
  data.frame(gene = gene,
             wilcoxp = wilcoxout$p.value,
             ttestp = tout$p.value)
}  )



compres <- dplyr::bind_rows(compres)



#make some plots
plotmat <- as.data.frame(t(ratiomat))

plotmat$Genotype <- c(rep(c('TKO', 'DKO'), each=3))

plotmatwide <- reshape2::melt(plotmat)
plotmatwide$value <- log2(plotmatwide$value)

ggplot(plotmatwide, aes(Genotype, value, col = Genotype))+
  geom_point()+
  facet_wrap(~variable, nrow = 2) +
  ylab('Log2(gene / Gapdh)')

compres
