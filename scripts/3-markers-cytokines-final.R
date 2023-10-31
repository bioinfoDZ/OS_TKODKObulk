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




#give it a label...
md$Sample
code <- data.frame(Sample = c("TL724M0T", "TL30F1T",  "TL744F2T", "DL761M2R", "DL579FOR", "DL608F1R"),
                   Code = c('TKO1', 'TKO2', 'TKO3', 'DKO1', 'DKO2', 'DKO3'))

md$Code <- code$Code
colnames(gem) <- code$Code
md$Sample <- code$Code







mhm <- markerheatmap(gem, geneannots$genes, metadata = md, geneannots = geneannots, res=res, cluster_rows = F,
                     name = 'scaled\nlog2\nnormalized\ncounts')





mhm



pdf('results/comparative-de/TKO-vs-DKO/4.downstream/mMCP_counter/markerheatmap_updated.pdf',
    height = 6, width = 6)


mhm

dev.off()





### update; 
# put TKO last
# put malignant OS on bottom



geneannots <- data.frame(genes = c('Ptprc'),
                         celltype = 'Immune',
                         color = 'black')


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


geneannots <- rbind(geneannots,
                    data.frame(genes = c('Sp7', 'Sox9', 'Ibsp', 'Col1a1', 'Col1a2'),
                               celltype = 'Malignant OS',
                               color = 'purple')
)



geneannots <- geneannots[geneannots$genes %in% rownames(gem),]
geneannots$celltype <- factor(geneannots$celltype, levels = unique(geneannots$celltype))


colnames( geneannots )[2] <- 'Celltype'




#give it a label...
md$Sample
code <- data.frame(Sample = c("TL724M0T", "TL30F1T",  "TL744F2T", "DL761M2R", "DL579FOR", "DL608F1R"),
                   Code = c('TKO1', 'TKO2', 'TKO3', 'DKO1', 'DKO2', 'DKO3'))

md$Code <- code$Code
colnames(gem) <- code$Code
md$Sample <- code$Code




gem2 <- gem[,c(4:6, 1:3)]

md2 <- md[match(colnames(gem2), md$Sample),]


mhm2 <- markerheatmap(gem2, geneannots$genes, metadata = md2, geneannots = geneannots, 
                      res=res, cluster_rows = F, cluster_columns=F,
                      #column_dend_reorder=F,
                     name = 'scaled\nlog2\nnormalized\ncounts')





mhm2



pdf('results/comparative-de/TKO-vs-DKO/4.downstream/mMCP_counter/markerheatmap_updated_FINALIZED_TKORIGHT_OSBOTTOM.pdf',
    height = 6, width = 6)


mhm2

dev.off()










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


gem2 <- gem[,c(4:6, 1:3)]

md2 <- md[match(colnames(gem2), md$Sample),]

FerrenaBulkRNAseq::heatmapplot(gem2, cytokines, metadata = md2)



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




markerheatmap(gem2, geneannots$genes, metadata = md2, geneannots = geneannots,  res=res, 
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




gem2 <- gem[,c(4:6, 1:3)]
md2 <- md[match(colnames(gem2), md$Sample),]



cyhm <- cytokineheatmap(gem2, geneannots$genes, metadata = md2, geneannots = geneannots,  res=res, 
                        cluster_rows = F, row_km = 2, cluster_columns=F,
                        name = 'scaled\nlog2\nnormalized\ncounts')

pdf('results/comparative-de/TKO-vs-DKO/4.downstream/mMCP_counter/cytokineheatmap_FINAL_TKOLAST.pdf',
    height = 6, width = 6)



cyhm
dev.off()

