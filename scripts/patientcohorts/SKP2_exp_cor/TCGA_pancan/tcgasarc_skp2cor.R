### parse and read in...?
library(TCGAbiolinks)
library(tidyverse)
library(SummarizedExperiment)
library(maftools) # v‘2.2.10’

library(survival)
library(survminer)


setwd('~/Dropbox/data/bangdata/target/pancan/')
# setwd('/home/Alex/data/pancan')
set.seed(2022)


#2022.may.24
packageVersion('TCGAbiolinks') #‘2.25.0’
packageVersion("tidyverse") #‘1.3.1'
packageVersion("SummarizedExperiment") #‘1.22.0’
packageVersion("maftools") #‘2.8.5’








### check sarc ###
proj = "TCGA-SARC"

samp <- list.files('parsed', pattern = proj, full.names = T)

message('Reading ', proj)

#read in data list
dl <- readRDS(samp)

#get the data
gem <- dl[[2]]
md <- dl[[3]]



### ignore subbtypes for now...?
#subtypes are a mess between two columns, 'primary_diagnosis' and 'paper_histology'
# lets also add in ICD10 code...?
# st <- as.data.frame(md[,c('barcode', 'primary_diagnosis', 'paper_histology' ,'tissue_or_organ_of_origin')])
# # st[is.na(st$paper_histology),"paper_histology"] <- st[is.na(st$paper_histology),"primary_diagnosis"]
# 
# write.csv('~/Desktop/tcgasarchistology.csv', x=st, quote = T, row.names = F)






### calculate MCP

# log counts...
mcp <- MCPcounter::MCPcounter.estimate(log(gem+1), featuresType = 'HUGO_symbols')


#prep to calculate cor with skp2 (TPM)
gem <- dl[[1]]
mcp <- as.data.frame(t(mcp))

#add genes
genes = c('SKP2')
genes <- genes[genes %in% rownames(gem)]
genedf <- data.frame( gem[match(genes, rownames(gem)),])
colnames(genedf) <- genes
mcp <- cbind(mcp, genedf)

#correlate
cr <- cor.test(mcp$SKP2, mcp$`Monocytic lineage`)

cr <- data.frame(proj=proj, r=cr$estimate, p=cr$p.value, meanMac=mean(mcp$`Monocytic lineage`), meanSKP2=mean(mcp$SKP2))


mm <- mcp

pdf('../finalized_results/skp2_correlations/tcgasarc_skp2_monocyte_cor.pdf', 
    height = 4, width = 4)

x = cor.test(mcp$SKP2, mcp$`Monocytic lineage`)
ggplot(mm, aes(SKP2, `Monocytic lineage`))+
  geom_point()+
  labs(subtitle = paste0('Cor = ', formatC(x$estimate, 3, format = "f"), '\nP = ', formatC(x$p.value, 3, format = "f")))+
  ylab('MCP-Counter Monocyte Lineage Score') + xlab('SKP2 (TPM)')+
  theme_classic2()

dev.off()




### find all skp2 cor genes
skp2 <- t(gem['SKP2',])

genes <- rownames(gem)
genes <- genes[!(genes %in% c('SKP2'))]

reslist <- list()
for(gene in genes){
  genevec <- t(gem[gene,])
  y <- cor.test(genevec, skp2)
  reslist[[gene]] <- data.frame(gene = gene, cor = y$estimate, p = y$p.value)
}

res <- dplyr::bind_rows(reslist)

#remove NAs, caaused by all 0 genes
res <- res[complete.cases(res),]

res$fdr <- p.adjust(res$p, method = 'fdr')

res <- res[order(res$cor, decreasing = T),]



sigres <- res[res$fdr<0.05,]



