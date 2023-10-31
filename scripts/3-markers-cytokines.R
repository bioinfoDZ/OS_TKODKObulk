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

### check markers, cytokines in TKO vs DKO de


deres <- 'results/comparative-de/TKO-vs-DKO/2.deresults/DEresults-normalized_count_matrix.csv'
res <- read.csv(deres)

md <- readxl::read_excel('data/metadata.xlsx')

#remove point / version info from ensembl IDs
res$ensembl_gene_id <- str_split_fixed(res$ensembl_gene_id, '\\.', 2)[,1]

#get normalized gene exp matrix
gem <- res[,match(md$Sample, colnames(res)) ]

#reset names
# if duplicate mgi symbol, always use ensebmlID
dups <- names( table(res$mgi_symbol)[table(res$mgi_symbol)>1] )
res[res$mgi_symbol %in% dups,"mgi_symbol"] <- res[res$mgi_symbol %in% dups,"ensembl_gene_id"] 

rownames(gem) <- res$mgi_symbol

#log transform the gene exp matrix...
gem <- log2(gem+1)


#marker list
# start with mMCP-Counter markers
library(mMCPcounter)
data("mMCPcounter_signatures")

mMCPcounter_signatures <-mMCPcounter_signatures[mMCPcounter_signatures$ENSEMBL.ID %in% res$ensembl_gene_id,]

table( mMCPcounter_signatures$Denomination )


# get mmcp res
mmcpres <- 'results/comparative-de/TKO-vs-DKO/4.downstream/mMCP_counter/mmcpres.csv'
mmcpres <- read.csv(mmcpres)


#get color pal for later
ctag <- aggregate(mMCP_estimate ~ celltype, mmcpres, mean)
ctag <- ctag[order(ctag$mMCP_estimate, decreasing = T),]

pal <- c( RColorBrewer::brewer.pal(Inf, "Set1"),
          RColorBrewer::brewer.pal(Inf, "Set2")
)

pal <- sample(pal, size = nrow(ctag), replace = F)
ctag$color <- pal 

#select significant only
mmcpres <- mmcpres[mmcpres$pval < 0.05,]

#sort celltypes pval by lowest pvalue...
ctpval <- data.frame(celltype = mmcpres$celltype, 
                     pval = mmcpres$pval)
ctpval <- ctpval[!duplicated(ctpval$celltype),]
ctpval <- ctpval[order(ctpval$pval),]

ctpval$color <- ctag[match(ctpval$celltype, ctag$celltype), "color"]


#get signatures of significant mmcp difference results
data("mMCPcounter_signatures")

#get color pal

mMCPcounter_signatures <- mMCPcounter_signatures[mMCPcounter_signatures$Denomination %in% ctpval$celltype,]

mMCPcounter_signatures <- mMCPcounter_signatures[order(mMCPcounter_signatures$Denomination),]

#keep genes in the dataset
mMCPcounter_signatures <- mMCPcounter_signatures[mMCPcounter_signatures$ENSEMBL.ID %in% res$ensembl_gene_id,]

#set up genes with annots
genes <- mMCPcounter_signatures$Gene.Symbol
geneannots <- data.frame(genes = genes,
                         celltype = mMCPcounter_signatures$Denomination)


geneannots$color <- ctpval[match(geneannots$celltype, ctpval$celltype), "color"]



#reorder genannots...
# this makes the heatmap the right order...
geneannotslist <- list()
for(ct in ctpval$celltype){
  ga <- geneannots[geneannots$celltype==ct,]
  geneannotslist[[ct]] <- ga
}
geneannots <- dplyr::bind_rows(geneannotslist)

#factorize also
# this makes the legend follow the right order
geneannots$celltype <- factor(geneannots$celltype, levels = ctpval$celltype)



#heatmap of cell types...
# want to plot cell types
# annotate genotype
# annotate marker class / cell type

markerheatmap <- function(gem, genes, geneannots, metadata, res, lfc_thres, pval_thres, scale, ...){
  
  require(ComplexHeatmap)
  
  
  
  # add lfc;pval info to gene
  if( missing(lfc_thres) ){lfc_thres <- 1}
  if( missing(pval_thres) ){pval_thres <- 0.05}
  if( missing(scale) ){scale <- T}
  
  
  # use res to get sig de genes
  res$symbol <- res$mgi_symbol
  res[abs(res$log2FoldChange) >  lfc_thres & res$padj < pval_thres, "symbol"] <- paste0('* ', res[abs(res$log2FoldChange) >  lfc_thres & res$padj < pval_thres, "symbol"], ' *' )
  
  sigres <- res[abs(res$log2FoldChange) >  lfc_thres & res$padj < pval_thres,]
  
  genes[genes %in% sigres$mgi_symbol] <-  na.omit( sigres[match(genes, sigres$mgi_symbol),"symbol"] )
  
  #changen ames...
  tmpgem <- gem; rownames(tmpgem) <- res$symbol
  
  #get gem, scale it
  # works best if you log transform before too
  tmpgem <- tmpgem[match(genes, rownames(tmpgem)),]
  
  if(scale==T){
    tmpgem <- t(scale(t(tmpgem)))
  } else{
    tmpgem <- as.matrix(tmpgem)
  }
  
  #set up the column annotation
  #annotation df, has sample (column) name and color
  annotdf <- data.frame(sample = colnames(tmpgem),
                        sampcheck = metadata[match(colnames(tmpgem), metadata$Sample), "Sample"],
                        cond = metadata[match(colnames(tmpgem), metadata$Sample), "Condition"],
                        color = metadata[match(colnames(tmpgem), metadata$Sample), "Color"])
  
  #define colors, need to use this list thing for complexheatmap
  hacol <- list(Condition = annotdf$Color) ; names(hacol[[1]]) <- annotdf$Condition
  
  #create the annotation object
  ha <- ComplexHeatmap::HeatmapAnnotation(Condition = annotdf[,3], col = hacol, name = 'Condition')
  
  
  
  #set up the marker annotation
  
  #define colors, need to use this list thing for complexheatmap
  hacol <- list(Celltype = geneannots$color) ; names(hacol[[1]]) <- geneannots$Celltype
  
  #create the annotation object
  ha_gene <- ComplexHeatmap::rowAnnotation(Celltype = geneannots$Celltype, col = hacol)
  
  
  ComplexHeatmap::Heatmap(tmpgem,
                          top_annotation = ha,
                          right_annotation = ha_gene,
                          ...
  )
  
  
}

geneannots <- geneannots[geneannots$genes %in% res$mgi_symbol,]
markerheatmap(gem, geneannots$genes, metadata = md, geneannots = geneannots, res = res, cluster_rows = F,
              name = 'scaled\nlog2\nnormalized\ncounts')


markerheatmap(gem, geneannots$genes, metadata = md, geneannots = geneannots, res = res, cluster_rows = F,
              scale=F, name = 'log2\nnormalized\ncounts')




### annotated markers ###


geneannots <- data.frame(genes = c('Ptprc'),
                         celltype = 'Immune',
                         color = 'black')

geneannots <- rbind(geneannots,
                    data.frame(genes = c('Col1a1', 'Sp7', 'Col1a2', 'Sox9', 'Ibsp'),
                               celltype = 'Malignant OS',
                               color = 'purple')
)

geneannots <- rbind(geneannots,
                    data.frame(genes = c('Cd68', 'Csf1r', 'Aif1', 'H2-Eb1', 'Adgre1'),
                               celltype = 'Macrophage',
                               color = '#FFFF33')
)



geneannots <- rbind(geneannots,
                    data.frame(genes = c('Cd14', 'Itgam'),
                               celltype = 'Monocyte',
                               color = 'yellow4')
)


geneannots <- rbind(geneannots,
                    data.frame(genes = c('Nfatc1', 'Acp5'),
                               celltype = 'Osteoclast',
                               color = 'orange')
)



geneannots <- rbind(geneannots,
                    data.frame(genes = c('Cd19', 'Ms4a1', 'Jchain'),
                               celltype = 'Bcell',
                               color = '#E5C494')
)

geneannots <- rbind(geneannots,
                    data.frame(genes = c('Cd3e', 'Cd3d', 'Cd8a', 'Cd4'),
                               celltype = 'Tcell',
                               color = '#4DAF4A')
)


geneannots <- rbind(geneannots,
                    data.frame(genes = c('Pecam1', 'Acta2'),
                               celltype = 'Endothelial',
                               color = 'red')
)


geneannots <- geneannots[geneannots$genes %in% rownames(gem),]
geneannots$celltype <- factor(geneannots$celltype, levels = unique(geneannots$celltype))



markerheatmap(gem, geneannots$genes, metadata = md, geneannots = geneannots, res=res, cluster_rows = F,
              name = 'scaled\nlog2\nnormalized\ncounts')

markerheatmap(gem, geneannots$genes, metadata = md, geneannots = geneannots, res=res, cluster_rows = F,
              scale=F, name = 'log2\nnormalized\ncounts')


mhm <- markerheatmap(gem, geneannots$genes, metadata = md, geneannots = geneannots, res=res, cluster_rows = F,
                     name = 'scaled\nlog2\nnormalized\ncounts')


jpeg(paste0('/Users/ferrenaa/Dropbox/data/bangdata/2021october-TKOvsDKO/results/comparative-de/TKO-vs-DKO/4.downstream/mMCP_counter/markerheatmap.jpg'),
     height = 6, width = 6, units = 'in', res = 300)


mhm

dev.off()


#M1 and M2
# https://www.sciencedirect.com/science/article/pii/S1074761314002283?via%3Dihub

#Csf1 vs GMCSF
# csf1, induces mature macrophage diff? not M1 vs M2 marker.
# cd14, immature, circulating monocyte marker?


#M1
#key cytokine, IFNG
#key transudcer, STAT1

#expression:
# TFs: IRF5, SOCS1
# cytokines, mouse and human: IL6, IL23a, Il12a, TNF
# chemokines (human only?): Cxcl10, il8, ccl5, cxcl9, cxcl10, cxcl11
# others: Marco, Mmp9, Nos2 hi/Arg1 low?

# genetics: Akt2? Klf6?

# what about Ifng receptor?

m1genes <- c('Ifng', 'Stat1', 'Il23a', 'Il12a', 'Tnf', 'Marco', 'Mmp9', 'Nos2', 'Akt2', 'Klf6')

m1genes <- data.frame(genes = m1genes,
                      celltype = 'M1_macrophage',
                      color = 'orange')



#M2
#key cytokine, Il4
#Key transducer, STAT6

#expression:
# TFs: IRF4, SOCS2
# cytokines: IL10, IL6
# chemokines: (Mouse only) CCL17, CCl24, CCL22
# chemokines; human only: CCL4, CCL13, CCL17, CCL18
# metabolic: Arg1? Nos2 lo? not really best marker
# others: TGFBR2

# genetics: Akt1? Klf4?


m2genes <- c('Il4', 'Stat6', 'Il10', 'Tgfb1', 'Ccl17', 'Arg1', 'Akt1', 'Klf4')




m2genes <- data.frame(genes = m2genes,
                      celltype = 'M2_macrophage',
                      color = 'skyblue')

geneannots <- rbind(m1genes, m2genes)

geneannots = geneannots[geneannots$genes %in% rownames(gem),]


markerheatmap(gem, geneannots$genes, geneannots = geneannots, md, res=res, km = 2)

markerheatmap(gem, geneannots$genes, metadata = md, geneannots = geneannots, res=res, km = 2,
              name = 'scaled\nlog2\nnormalized\ncounts')

markerheatmap(gem, geneannots$genes, metadata = md, geneannots = geneannots, res=res, km = 2,
              scale=F, name = 'log2\nnormalized\ncounts')





#pro and anti inflam cytokines


inflam <- c('Tnf', 'Ifng', 'Il6')
anti <- c('Il10', 'Tgfb1', 'Tgfb2', 'Tgfb3')


m1genes <- data.frame(genes = inflam,
                      celltype = 'Pro-Inflammatory',
                      color = 'red')
m2genes <- data.frame(genes = anti,
                      celltype = 'Anti-Inflammatory',
                      color = 'skyblue')

geneannots <- rbind(m1genes, m2genes)

geneannots = geneannots[geneannots$genes %in% rownames(gem),]

geneannots$celltype <- factor(geneannots$celltype, levels = unique(geneannots$celltype))



markerheatmap(gem, geneannots$genes, metadata = md, geneannots = geneannots, res=res, cluster_rows = F)





### t cell cytotoxic vs exhaustion

gzm <- sort(res$mgi_symbol[grepl('Gzm', res$mgi_symbol)])

cytotoxic <- c('Tnf', 'Ifng', 'Prf1', gzm)

exhaustion <- c('Pdcd1', 'Ctla4', 'Tox', 'Foxp3')



m1genes <- data.frame(genes = cytotoxic,
                      celltype = 'Effector T cell',
                      color = 'red')
m2genes <- data.frame(genes = exhaustion,
                      celltype = 'Exhausted T cell',
                      color = 'skyblue')

geneannots <- rbind(m1genes, m2genes)

geneannots = geneannots[geneannots$genes %in% rownames(gem),]

geneannots$celltype <- factor(geneannots$celltype, levels = unique(geneannots$celltype))



markerheatmap(gem, geneannots$genes, metadata = md, geneannots = geneannots, res=res, cluster_rows = F)








#read in panglaodb markers...

# pdb <- read.table('~/Dropbox/data/PanglaoDB_markers_27_Mar_2020.tsv', sep='\t', header = T)
# pdb <- read.table('~/Dropbox/data/general_data_utilities/panglaodb/PanglaoDB_markers_27_Mar_2020.tsv', sep='\t', header = T)
# 
# 
# 
# #select only mouse or mouse/human both genes
# pdb <- pdb[grepl('Mm', pdb$species),]
# 
# #select only fulll rows
# pdb <- pdb[complete.cases(pdb),]
# 
# 
# 
# #select few cell types
# cts <- c('Macrophages', 'Monocytes', 'Myeloid-derived suppressor cells', )


# the download is buggest....


# just pick some cell types, and get top 3 markers....






geneannots <- data.frame(genes = c('Ptprc'),
                         celltype = 'Immune',
                         color = 'black')

geneannots <- rbind(geneannots,
                    data.frame(genes = c('Sp7', 'Sox9', 'Ibsp', 'Col1a1', 'Col1a2'),
                               celltype = 'Malignant OS',
                               color = 'purple')
)


# 
# geneannots <- rbind(geneannots,
#                     data.frame(genes = c('Nlrp3', 'Fgfr1', 'Cx3cl1'),
#                                celltype = 'Osteoblast',
#                                color = 'pink')
# )
# 
# geneannots <- rbind(geneannots,
#                     data.frame(genes = c('Sost', 'Pfn1', 'Efnb2'),
#                                celltype = 'Osteocyte',
#                                color = 'cyan')
# )
# 
# geneannots <- rbind(geneannots,
#                     data.frame(genes = c('Emp1', 'Cytl1', 'Lcn2'),
#                                celltype = 'Chondrocyte',
#                                color = 'blue')
# )


geneannots <- rbind(geneannots,
                    data.frame(genes = c('Cd68', 'H2-Eb1', 'Aif1', 'Adgre1'),
                               celltype = 'Macrophage',
                               color = '#FFFF33')
)



geneannots <- rbind(geneannots,
                    data.frame(genes = c('Cd14', 'Itgam', 'Csf1r', 'Ccr2', 'Cx3cr1'),
                               celltype = 'Monocyte',
                               color = 'yellow4')
)


geneannots <- rbind(geneannots,
                    data.frame(genes = c('Nfatc1', 'Mmp9', 'Ctsk'),
                               celltype = 'Osteoclast',
                               color = 'orange')
)



geneannots <- rbind(geneannots,
                    data.frame(genes = c('Cd19', 'Ms4a1', 'Jchain'),
                               celltype = 'Bcell',
                               color = '#E5C494')
)

geneannots <- rbind(geneannots,
                    data.frame(genes = c('Trbc2','Cd3e', 'Cd3d', 'Cd8a', 'Cd4'),
                               celltype = 'Tcell',
                               color = '#4DAF4A')
)


geneannots <- rbind(geneannots,
                    data.frame(genes = c('Acta2', 'Mylk', 'Myl9'),
                               celltype = 'Smooth Muscle',
                               color = 'brown')
)


geneannots <- rbind(geneannots,
                    data.frame(genes = c('Pecam1', 'Egfl7', 'Id3'),
                               celltype = 'Endothelial',
                               color = 'red')
)


geneannots <- geneannots[geneannots$genes %in% rownames(gem),]
geneannots$celltype <- factor(geneannots$celltype, levels = unique(geneannots$celltype))


colnames( geneannots )[2] <- 'Celltype'




mhm <- markerheatmap(gem, geneannots$genes, metadata = md, geneannots = geneannots, res=res, cluster_rows = F,
              name = 'scaled\nlog2\nnormalized\ncounts')




pdf('results/comparative-de/TKO-vs-DKO/4.downstream/mMCP_counter/markerheatmap_updated.pdf',
    height = 6, width = 6)


mhm

dev.off()







###### CCL genes and others



geneannots <- data.frame(genes = str_sort(grep(rownames(gem), pattern = 'Ccl', ignore.case = F, value = T), numeric = T),
                         celltype = 'CCL genes',
                         color = 'orange')


geneannots <- geneannots[geneannots$genes %in% rownames(gem),]
geneannots$celltype <- factor(geneannots$celltype, levels = unique(geneannots$celltype))



markerheatmap(gem, geneannots$genes, metadata = md, geneannots = geneannots, res=res, cluster_rows = F,
              name = 'scaled\nlog2\nnormalized\ncounts')



### il genes excluding ilrs...


genes <-str_sort(grep(rownames(gem), pattern = '^Il', ignore.case = F, value = T), numeric = T)

genes <- c(genes, 'Tgfb1', 'Tgfb2', 'Tnf')

geneannots <- data.frame(genes = genes,
                         celltype = 'IL genes and others',
                         color = 'orange')
geneannots <- geneannots[!(grepl(geneannots$genes, pattern = 'r')),]


geneannots <- geneannots[geneannots$genes %in% rownames(gem),]
geneannots$celltype <- factor(geneannots$celltype, levels = unique(geneannots$celltype))



markerheatmap(gem, geneannots$genes, metadata = md, geneannots = geneannots, res=res, cluster_rows = F,
              name = 'scaled\nlog2\nnormalized\ncounts')



### ILRs excluding 

genes <-str_sort(grep(rownames(gem), pattern = '^Il', ignore.case = F, value = T), numeric = T)

geneannots <- data.frame(genes = genes,
                         celltype = 'ILR genes',
                         color = 'orange')
geneannots <- geneannots[grepl(geneannots$genes, pattern = 'r'),]


geneannots <- geneannots[geneannots$genes %in% rownames(gem),]
geneannots$celltype <- factor(geneannots$celltype, levels = unique(geneannots$celltype))



markerheatmap(gem, geneannots$genes, metadata = md, geneannots = geneannots, res=res, cluster_rows = F,
              name = 'scaled\nlog2\nnormalized\ncounts')




### CXCL

geneannots <- data.frame(genes = str_sort(grep(rownames(gem), pattern = 'Cxcl', ignore.case = F, value = T), numeric = T),
                         celltype = 'CXCL genes',
                         color = 'orange')


geneannots <- geneannots[geneannots$genes %in% rownames(gem),]
geneannots$celltype <- factor(geneannots$celltype, levels = unique(geneannots$celltype))



markerheatmap(gem, geneannots$genes, metadata = md, geneannots = geneannots, res=res, cluster_rows = F,
              name = 'scaled\nlog2\nnormalized\ncounts')






### CXCR

geneannots <- data.frame(genes = str_sort(grep(rownames(gem), pattern = 'Cxcr', ignore.case = F, value = T), numeric = T),
                         celltype = 'CXCR genes',
                         color = 'orange')


geneannots <- geneannots[geneannots$genes %in% rownames(gem),]
geneannots$celltype <- factor(geneannots$celltype, levels = unique(geneannots$celltype))



markerheatmap(gem, geneannots$genes, metadata = md, geneannots = geneannots, res=res, cluster_rows = F,
              name = 'scaled\nlog2\nnormalized\ncounts')









##### get a list of all significantly DE cytokines/chemokines

#get KEGG cytokin e pathway... seems to be human specific... so add to grep search
pathways <- readRDS('~/Dropbox/data/general_data_utilities/msigdb/msigdb_mouse.2022.may.12.rds')

cytokegg <- pathways[pathways$gs_name=='KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION',]$gene_symbol

#exclude receptors and CD molecules 
cytokegg <- cytokegg[!grepl('r', cytokegg, ignore.case = T)]
#add back in Prl, 
cytokegg <- c(cytokegg, 'Prl')
cytokegg <- str_sort(cytokegg, numeric = T)

#exclude CD mollecules
cytokegg <- cytokegg[!grepl('CD', cytokegg, ignore.case = T)]



### my own list...

cytokines <- c()

# ccls
cytokines <- c(cytokines, str_sort(grep(rownames(gem), pattern = 'Ccl', ignore.case = F, value = T), numeric = T) )

#ILs, excluding ILrs
genes <- str_sort(grep(rownames(gem), pattern = '^Il', ignore.case = F, value = T), numeric = T)
genes <- genes[!(grepl(genes, pattern = 'r'))]
genes <- genes[!(grepl(genes, pattern = 't'))]
genes <- genes[!(grepl(genes, pattern = 'bp'))]
genes <- genes[!(grepl(genes, pattern = 'v'))]
genes <- genes[!(grepl(genes, pattern = 'kap'))]
genes <- genes[!(grepl(genes, pattern = 'Ilf'))]

cytokines <- c(cytokines, genes)



# cxcls
cytokines <- c(cytokines, str_sort(grep(rownames(gem), pattern = 'Cxcl', ignore.case = F, value = T), numeric = T) )


#combine the lists

cytokines <- unique( c(cytokegg, cytokines) )
cytokines <- str_sort(cytokines, numeric = T)


### final filtering ###
#get rid of "FLT1,2,4...
cytokines <- cytokines[!grepl('Flt[[:digit:]]$', cytokines)]

#get rid of fas (receptor), but keep fas Ligand
cytokines <- cytokines[!grepl('Fas$', cytokines)]

#get only significant
sig <- res[res$padj < 0.05,]
sig <- sig[abs(sig$log2FoldChange) > 1,]
siggenes <- sig$mgi_symbol
cytokines <- cytokines[cytokines %in% siggenes]



# markerheatmap(gem, cytokines, metadata = md,  res=res, cluster_rows = F,
#               name = 'scaled\nlog2\nnormalized\ncounts')


FerrenaBulkRNAseq::heatmapplot(gem, cytokines, metadata = md)



geneannots <- data.frame(genes = c('Cxcl12', 'Pf4', 'Cxcl13' ,'Ccl6', 'Cxcl16', 'Cxcl14', 'Ccl2', 'Ccl7'),
                         celltype = 'Monocyte migration',
                         color = 'orange')

geneannots <- rbind(geneannots, 
                    data.frame(genes = c('Ccl4', 'Ccl3', 'Il1b', 'Ccl9', 'Osm', 'Il6', 'Il12b', 'Tslp', 'Il33', 'Cxcl5'),
                         celltype = 'Proinflammatory',
                         color = 'firebrick') )

geneannots <- rbind(geneannots, 
                    data.frame(genes = c('Ccl22', 'Tnfsf14', 'Tnfsf9', 'Il7'),
                         celltype = 'T cell migration/stimulation',
                         color = 'yellow') )

geneannots <- rbind(geneannots, 
                    data.frame(genes = c('Ccl11', 'Csf3'),
                         celltype = 'Other myeloid function',
                         color = 'purple'))


geneannots <- rbind(geneannots, 
                    data.frame(genes = c('Ltb', 'Tgfb3', 'Lif', 'Cxcl1'),
                         celltype = 'Unclear/pleiotropic',
                         color = 'grey'))

geneannots <- rbind(geneannots, 
                    data.frame(genes = c('Tnfsf10'),
                         celltype = 'Apoptosis',
                         color = 'black'))


geneannots <- rbind(geneannots, 
                    data.frame(genes = c('Pdgfb','Pdgfa', 'Flt3l'),
                         celltype = 'Mitogen',
                         color = 'steelblue'))


geneannots <- rbind(geneannots, 
                    data.frame(genes = c('Vegfa','Vegfd'),
                         celltype = 'Angiogenesis',
                         color = 'green'))


colnames(geneannots)[2] <- 'Function'

markerheatmap(gem, geneannots$genes, metadata = md, geneannots = geneannots,  res=res, 
              cluster_rows = T, row_km = 2,
              name = 'scaled\nlog2\nnormalized\ncounts')




# heatmap of cell types...
# want to plot cell types
# annotate genotype
# annotate marker class / cell type

cytokineheatmap <- function(gem, genes, geneannots, metadata, res, lfc_thres, pval_thres, scale, ...){
  
  require(ComplexHeatmap)
  
  
  
  # add lfc;pval info to gene
  if( missing(lfc_thres) ){lfc_thres <- 1}
  if( missing(pval_thres) ){pval_thres <- 0.05}
  if( missing(scale) ){scale <- T}
  
  
  # use res to get sig de genes
  res$symbol <- res$mgi_symbol
  #res[abs(res$log2FoldChange) >  lfc_thres & res$padj < pval_thres, "symbol"] <- paste0('* ', res[abs(res$log2FoldChange) >  lfc_thres & res$padj < pval_thres, "symbol"], ' *' )
  
  sigres <- res[abs(res$log2FoldChange) >  lfc_thres & res$padj < pval_thres,]
  
  genes[genes %in% sigres$mgi_symbol] <-  na.omit( sigres[match(genes, sigres$mgi_symbol),"symbol"] )
  
  #changen ames...
  tmpgem <- gem; rownames(tmpgem) <- res$symbol
  
  #get gem, scale it
  # works best if you log transform before too
  tmpgem <- tmpgem[match(genes, rownames(tmpgem)),]
  
  if(scale==T){
    tmpgem <- t(scale(t(tmpgem)))
  } else{
    tmpgem <- as.matrix(tmpgem)
  }
  
  #set up the column annotation
  #annotation df, has sample (column) name and color
  annotdf <- data.frame(sample = colnames(tmpgem),
                        sampcheck = metadata[match(colnames(tmpgem), metadata$Sample), "Sample"],
                        cond = metadata[match(colnames(tmpgem), metadata$Sample), "Condition"],
                        color = metadata[match(colnames(tmpgem), metadata$Sample), "Color"])
  
  #define colors, need to use this list thing for complexheatmap
  hacol <- list(Condition = annotdf$Color) ; names(hacol[[1]]) <- annotdf$Condition
  
  #create the annotation object
  ha <- ComplexHeatmap::HeatmapAnnotation(Condition = annotdf[,3], col = hacol, name = 'Condition')
  
  
  
  #set up the marker annotation
  
  #define colors, need to use this list thing for complexheatmap
  hacol <- list(Function = geneannots$color) ; names(hacol[[1]]) <- geneannots$Function
  
  #create the annotation object
  ha_gene <- ComplexHeatmap::rowAnnotation(Function = geneannots$Function, col = hacol)
  
  
  ComplexHeatmap::Heatmap(tmpgem,
                          top_annotation = ha,
                          right_annotation = ha_gene,
                          ...
  )
  
  
}


geneannots <- geneannots[order(geneannots$Function),]


#give it a label...
md$Sample
code <- data.frame(Sample = c("TL724M0T", "TL30F1T",  "TL744F2T", "DL761M2R", "DL579FOR", "DL608F1R"),
           Code = c('TKO1', 'TKO2', 'TKO3', 'DKO1', 'DKO2', 'DKO3'))

md$Code <- code$Code
colnames(gem) <- code$Code
md$Sample <- code$Code


cyhm <- cytokineheatmap(gem, geneannots$genes, metadata = md, geneannots = geneannots,  res=res, 
                        cluster_rows = F, row_km = 2, cluster_columns=F,
                        name = 'scaled\nlog2\nnormalized\ncounts')

pdf('results/comparative-de/TKO-vs-DKO/4.downstream/mMCP_counter/cytokineheatmap.pdf',
    height = 6, width = 6)



cyhm
dev.off()
