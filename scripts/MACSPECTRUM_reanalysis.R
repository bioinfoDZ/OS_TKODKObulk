library(tidyverse)
library(Seurat)
library(FerrenaSCRNAseq)


set.seed(2022)


setwd('~/Dropbox/data/deyou/macspectrum/')

samps <- list.files('rawdata_geo/')



sobjlist <- list()
autofilterlist <- list()

for(samp in samps){
  message(samp)
  
  fp <- paste0('rawdata_geo/', samp)
  sobj <- CreateSeuratObject( Read10X(fp), project = samp )
  
  #mito content, add to metadata
  mito.features <- grep(pattern = "^mt-", x = rownames(x = sobj), value = TRUE, ignore.case = T)
  sobj[["percent.mito"]] <- Seurat::PercentageFeatureSet(sobj, features = mito.features)
  
  #normalize and cluster
  suppressWarnings(sobj <- Seurat::SCTransform(sobj, verbose = T))
  
  sobj <- Seurat::RunPCA(object = sobj, verbose = F)
  
  sobj <- Seurat::FindNeighbors(object = sobj, dims = 1:20, verbose = F)
  sobj <- Seurat::FindClusters(object = sobj, resolution = 0.1, verbose = F, algorithm = 4)
  
  sobj <- RunUMAP(sobj, dims = 1:20)
  
  
  
  ### run auto filter ###
  reportlist <- FerrenaSCRNAseq::automatedfiltering(sobj, clusters = 'SCT_snn_res.0.1',
                                                    iterativefilter.mito = F)
  
  autofilterlist[[samp]] <- reportlist
  
  #add autofilter results to metadata
  autofilterres <- reportlist[[1]]
  sobj$filteredout <- autofilterres$filteredout
  sobj$filterreason <- autofilterres$filterreason
  
  #filter
  goodcells <- autofilterres[autofilterres$filteredout == 'No', 'barcodes']
  sobj <- sobj[,goodcells]
  
  #reprocess
  #normalize and cluster
  suppressWarnings(sobj <- Seurat::SCTransform(sobj, verbose = T))
  
  sobj <- Seurat::RunPCA(object = sobj, verbose = F)
  
  sobj <- Seurat::FindNeighbors(object = sobj, dims = 1:20, verbose = F)
  sobj <- Seurat::FindClusters(object = sobj, resolution = 0.1, verbose = F, algorithm = 4)
  
  sobj <- RunUMAP(sobj, dims = 1:20)
  
  #use doublet filtering
  dfdf <- FerrenaSCRNAseq::doubletfinderwrapper(sobj, clusters = 'SCT_snn_res.0.1')
  
  
  #add doublet info to seurat
  sobj$DoubletFinderClassification <- dfdf$DoubletFinderClassification
  
  
  #remove all doublets
  md <- sobj@meta.data
  md <- md[md$DoubletFinderClassification == 'Singlet' ,]
  goodcells <- rownames(md)
  
  sobj <- sobj[,goodcells]
  
  
  #reprocess without doublets
  suppressWarnings(sobj <- Seurat::SCTransform(sobj, verbose = T))
  
  sobj <- Seurat::RunPCA(object = sobj, verbose = F)
  
  sobj <- Seurat::FindNeighbors(object = sobj, dims = 1:20, verbose = F)
  sobj <- Seurat::FindClusters(object = sobj, resolution = 0.1, verbose = F, algorithm = 4)
  
  sobj <- RunUMAP(sobj, dims = 1:20)
  
  
  sobjlist[[samp]] <- sobj
}





beepr::beep()


#add cell cycle scoring...

sobjlist <- lapply( sobjlist,
                    function(x) { CellCycleScoring(x, s.features = Seurat::cc.genes.updated.2019$s.genes, g2m.features = Seurat::cc.genes.updated.2019$g2m.genes)}
)


#save sobjlist
saveRDS(sobjlist, 'sobjlist.rds')




#### get avgs....
# just need unnormalized UMIs
sobjlist <- readRDS('sobjlist.rds')

avglist <- lapply(sobjlist, function(x){  rowMeans(x@assays$RNA@counts)  })



rm(sobjlist, sobj, reportlist)

# saveRDS(avglist, 'avglist.rds')
# avglist <- readRDS('avglist.rds')

#keep only intersect of gene names
genes <- Reduce(intersect, lapply(avglist, FUN = names))

avglist <- lapply(avglist, function(x){x[match(genes, names(x))]})


#make to df
avgdf <- data.frame(avglist)



#get ensembl ids
anno <- as.data.frame(annotables::grcm38)

#filter...
anno <- anno[anno$symbol %in% genes,]
avgdf <- avgdf[rownames(avgdf) %in% anno$symbol,]




#apply a very small cutoff... does reduce a lot..
# avgdf <- avgdf[rowMeans(avgdf)>0.5,]

#get genes after cutoff
genes <- rownames(avgdf)
anno <- anno[match(genes, anno$symbol),]


#get rid of dups...
# no more
table(duplicated(rownames(avgdf)))
table(duplicated(anno$ensgene))


#chage rownames to ensembl
rownames(avgdf) <- anno$ensgene


#format properly...
mat <- as.data.frame(avgdf)
mat <- cbind(rownames(mat), mat)
colnames(mat)[1] <- 'Ensembl_ID'

#finally, make features df
features <- c('M0', 'M1', 'M2')


### save

write.csv(mat, paste0('mat.csv'), quote = F, row.names = F)
write.csv(features, paste0('features.csv'), quote = F, row.names = F)




### also, let's try to bind this df together with the results from others; TKO vs DKO

macspectrumdir <- 'results/comparative-de/TKO-vs-DKO/4.downstream/macspectrum'
macspectrumdir <- paste0('~/Dropbox/data/bangdata/2021october-TKOvsDKO/', macspectrumdir)




mat2 <- read.csv(paste0(macspectrumdir, '/tpm.csv'))
f2 <- read.csv(paste0(macspectrumdir, '/feature.csv'))


# read ub thge m1 m2 m0 

mat <- read.csv(paste0('mat.csv'))
f <- read.csv(paste0('features.csv'))



#the lib sizes are very differet, but watever let's just do it...

#intersect of genes...
intgenes <- intersect(mat[,1], mat2[,1])

mat <- mat[match(intgenes, mat[,1]) ,]
mat2 <- mat2[match(intgenes, mat2[,1]), ]



#maybe quantile norm?
# try after binding?
#gemt <- as.data.frame(preprocessCore::normalize.quantiles(as.matrix(gem)))

# mat[,-1] <- as.data.frame(preprocessCore::normalize.quantiles(as.matrix(mat[,-1])))
# 
# mat2[,-1] <- as.data.frame(preprocessCore::normalize.quantiles(as.matrix(mat2[,-1])))


colnames(f) <- 'x'
colnames(f2) <- 'x'

#bind mats
bigmat <- cbind(mat, mat2[,-1])
bigfs_in_chat <- rbind(f,f2)

#quantile norm?
bigmat[,-1] <- as.data.frame(preprocessCore::normalize.quantiles(as.matrix(bigmat[,-1])))



#try to write it....
write.csv(bigmat, paste0('bigmat_bound_macspec_tkodko.csv'), quote = F, row.names = F)
write.csv(bigfs_in_chat, paste0('bigfeatures__bound_macspec_tkodko.csv'), quote = F, row.names = F)





### just cor...
ComplexHeatmap::Heatmap(cor(bigmat[,-1], method = 'spearman'))


### remove dko

bigmat_nodko <- bigmat[,1:7]
bigf_nodko <- bigfs_in_chat[1:6,]

ComplexHeatmap::Heatmap(cor(bigmat_nodko[,-1], method = 'spearman'))

#try to write it....
write.csv(bigmat_nodko, paste0('bigmat_bound_macspec_tkoOnly.csv'), quote = F, row.names = F)
write.csv(bigf_nodko, paste0('bigfeatures__bound_macspec_tkoOnly.csv'), quote = F, row.names = F)





### means, get cor???

cor( rowMeans(mat2[,2:4]), bigmat$m0, method = 'spearman')
cor( rowMeans(mat2[,2:4]), bigmat$m1, method = 'spearman')
cor( rowMeans(mat2[,2:4]), bigmat$m2, method = 'spearman')





#### try RISC, marker analysis??


#first just concat...
sobjlist <- readRDS("sobjlist.rds")

sobjint <- merge(sobjlist[[1]], sobjlist[[2]])
sobjint <- merge(sobjint, sobjlist[[3]])


DefaultAssay(sobjint) <- 'RNA'
sobjint$seurat_clusters_singlesample <- sobjint$seurat_clusters


#reprocess without doublets
suppressWarnings(sobjint <- Seurat::SCTransform(sobjint, verbose = T))

sobjint <- Seurat::RunPCA(object = sobjint, verbose = F)

sobjint <- Seurat::FindNeighbors(object = sobjint, dims = 1:20, verbose = F)
sobjint <- Seurat::FindClusters(object = sobjint, resolution = 0.1, verbose = F, algorithm = 4)

sobjint <- RunUMAP(sobjint, dims = 1:20)

beepr::beep()

DimPlot(sobjint, group.by = 'orig.ident') + DimPlot(sobjint, label = T) + 
  DimPlot(sobjint, group.by = 'Phase') + FeaturePlot(sobjint, 'nCount_RNA')



###  m2 and m0 similar, m1 is different
# batch or biology? try risc...




library(RISC)
ncore <- parallel::detectCores() - 4

#get object list
sobjlist <- readRDS("sobjlist.rds")


#prepare the RISC functions
process0 <- function(obj0){
  # Filter cells and genes
  obj0 = scFilter(obj0, min.UMI = 1000, max.UMI = 40000, min.gene = 200, min.cell = 0)
  # Normalize the raw counts
  obj0 = scNormalize(obj0, ncore = ncore)
  # Find highly variable genes
  obj0 = scDisperse(obj0)
  # print(length(obj0@vargene))
  return(obj0)
}


rm(sobjint, mat,mat2,bigmat, bigmat_nodko, intgenes)

#### integration with RISC #######


risclist <- list()

for(sampidx in c(1:length(sobjlist)) ){
  
  samp <- names(sobjlist)[sampidx]
  sobj <- sobjlist[[sampidx]]
  
  message(' - Prepping RISC object')
  #get RNA counts matrix
  mat0 <- sobj@assays$RNA@counts
  
  #get metadata, keep important columns
  coldata0 <- sobj@meta.data
  
  #add cell names and rename clustering/cellname columns
  
  #get barcodes, strip numeric suffix, and add samplename prefix
  barcodes <- stringr::str_split_fixed(rownames(coldata0), '-', 2)[,1]
  barcodes <- paste0(coldata0$orig.ident, '.',barcodes)
  
  coldata0 <- cbind(barcodes, coldata0)
  
  colnames(coldata0)[ncol(coldata0)] <- paste0('Seurat-individualsampleclust-', colnames(coldata0)[ncol(coldata0)])
  
  
  #make the rowdatadf...
  rowdata0 = data.frame(Symbol = rownames(mat0), row.names = rownames(mat0))
  
  #make the risc object
  dat0 = readsc(mat0, coldata0, rowdata0, is.filter = T)
  
  rm(sobj, mat0, rowdata0, coldata0, barcodes)
  
  #process the risc object
  message(' - Processing RISC object')
  
  dat0 <- process0(dat0)
  
  
  
  
  risclist[[samp]] <- dat0
  rm(dat0)
  
}

beepr::beep()

rm(sobjlist)



#get the intersect of gene names
var0 <- Reduce(intersect, lapply(risclist, FUN = function(x){x@rowdata$Symbol}))
# var0 <- Reduce(unique, lapply(risclist, FUN = function(x){x@rowdata$Symbol}))


#look at it and choose...
InPlot(risclist, var.gene = var0, Std.cut = 0.95, ncore = ncore) 



# set risc params
ref <- 1
eigens <- 15

#need to rearrange the list...

if(ref != 1){
  
  data0 <- list(risclist[[ref]])
  names(data0) <- names(risclist)[ref]
  
  for(i in 1:length(risclist)){
    if(i != ref){
      name = names(risclist)[i]
      data0[[name]] <- risclist[[i]]
    }
    
  }
  
} else{
  data0 <- risclist
}

rm(risclist)


#actually integrate
data0 = scMultiIntegrate(
  objects = data0, eigens = eigens, add.Id = NULL, var.gene = var0,
  # method = "RPCI", 
  align = 'OLS', npc = 50, adjust = TRUE,
  ncore = ncore, 
  #do.fast = "AUTO"
)



data0 = scUMAP(data0, npc = eigens, use = "PLS")

saveRDS(data0, 'RISC_int_1.6.rds')




#plots
d1=RISC::DimPlot(data0, colFactor = 'orig.ident', Alpha = 0.3, size = 0.5)

d1


beepr::beep()

# saveRDS(data0, 'RISC_int.rds')

data0 <- readRDS("RISC_int.rds")



# get the metadata
md <- data0@coldata

# get the gene expression matrix
gem <- data0@assay$logcount
gem <- do.call(cbind, gem)

#get UMAP coords
umap <- data0@DimReduction$cell.umap

#make seurat object, for nice pre-made plotting
sobj <- CreateSeuratObject(counts = gem, meta.data = md)
sobj[['umap']] <- Seurat::CreateDimReducObject(umap, key = 'umap')


rm(data0, md, umap,gem)



#find markers of m0, m1, m2
sobj <- SetIdent(sobj, value = sobj$orig.ident)



m <- FindAllMarkers(sobj, only.pos = T)



#pairwise
### m0 ###
m0_tot <- FindMarkers(sobj, ident.1 = 'm0', only.pos = T)
m0_vs_m1 <- FindMarkers(sobj, ident.1 = 'm0', ident.2 = 'm1', only.pos = T)
m0_vs_m2 <- FindMarkers(sobj, ident.1 = 'm0', ident.2 = 'm2', only.pos = T)

m0_spec <- m0_tot[rownames(m0_tot) %in% rownames(m0_vs_m1),]
m0_spec <- m0_tot[rownames(m0_tot) %in% rownames(m0_vs_m2),]
m0_spec <- m0_spec[order(m0_spec$avg_log2FC, decreasing = T),]

### m1
m1_tot <- FindMarkers(sobj, ident.1 = 'm1', only.pos = T)
m1_vs_m0 <- FindMarkers(sobj, ident.1 = 'm1', ident.2 = 'm0', only.pos = T)
m1_vs_m2 <- FindMarkers(sobj, ident.1 = 'm1', ident.2 = 'm2', only.pos = T)

m1_spec <- m1_tot[rownames(m1_tot) %in% rownames(m1_vs_m0),]
m1_spec <- m1_tot[rownames(m1_tot) %in% rownames(m1_vs_m2),]
m1_spec <- m1_spec[order(m1_spec$avg_log2FC, decreasing = T),]

### m2
m2_tot <- FindMarkers(sobj, ident.1 = 'm2', only.pos = T)
m2_vs_m0 <- FindMarkers(sobj, ident.1 = 'm2', ident.2 = 'm0', only.pos = T)
m2_vs_m1 <- FindMarkers(sobj, ident.1 = 'm2', ident.2 = 'm1', only.pos = T)

m2_spec <- m2_tot[rownames(m2_tot) %in% rownames(m2_vs_m0),]
m2_spec <- m2_tot[rownames(m2_tot) %in% rownames(m2_vs_m1),]
m2_spec <- m2_spec[order(m2_spec$avg_log2FC, decreasing = T),]

beepr::beep()

#exclude overlapping genes; there are not many of these
m0_spec <- m0_spec[!(rownames(m0_spec) %in% rownames(m1_spec)),]
m0_spec <- m0_spec[!(rownames(m0_spec) %in% rownames(m2_spec)),]

m1_spec <- m1_spec[!(rownames(m1_spec) %in% rownames(m0_spec)),]
m1_spec <- m1_spec[!(rownames(m1_spec) %in% rownames(m2_spec)),]

m2_spec <- m2_spec[!(rownames(m2_spec) %in% rownames(m0_spec)),]
m2_spec <- m2_spec[!(rownames(m2_spec) %in% rownames(m1_spec)),]


### add stringency: more conservative lfc cutooff
m0_spec <- m0_spec[m0_spec$avg_log2FC>1,]
m1_spec <- m1_spec[m1_spec$avg_log2FC>1,]
m2_spec <- m2_spec[m2_spec$avg_log2FC>1,]

#save it all
mlist <- list(m0_tot=m0_tot, m0_vs_m1=m0_vs_m1, m0_vs_m2=m0_vs_m2,
              m1_tot=m1_tot, m1_vs_m0=m1_vs_m0, m1_vs_m2=m1_vs_m2,
              m2_tot=m2_tot, m2_vs_m0=m2_vs_m0, m2_vs_m1=m2_vs_m1,
              m0_spec=m0_spec,
              m1_spec=m1_spec,
              m2_spec=m2_spec)

saveRDS(mlist, 'mlist_1.6.rds')

mlist <- readRDS('mlist.rds')

m0_spec <- mlist[['m0_spec']]
m1_spec <- mlist[['m1_spec']]
m2_spec <- mlist[['m2_spec']]



#module scores
md <- sobj@meta.data[,1:23]
sobj@meta.data <- md
sobj <- AddModuleScore(sobj, features = list(m0 = rownames(m0_spec) ), name = 'm0_' )
sobj <- AddModuleScore(sobj, features = list(m1 = rownames(m1_spec) ), name = 'm1_' )
sobj <- AddModuleScore(sobj, features = list(m2 = rownames(m2_spec) ), name = 'm2_' )




DimPlot(sobj, label = T) + FeaturePlot(sobj,  c('m0_1', 'm1_1', 'm2_1'), label = T, ncol = 1)


#bind them together
m0_spec$gene <- rownames(m0_spec); m0_spec$cluster <- 'm0'
m1_spec$gene <- rownames(m1_spec); m1_spec$cluster <- 'm1'
m2_spec$gene <- rownames(m2_spec); m2_spec$cluster <- 'm2'

m <- rbind(m0_spec, m1_spec, m2_spec)


n=10
top <- m %>% group_by(cluster) %>% top_n(n = n, wt = avg_log2FC)

sobj <- ScaleData(sobj, features = m$gene)
DoHeatmap(sobj, features = top$gene, raster = T)




### check all CCL genes
sobj <- ScaleData(sobj)
allgenes <- rownames(sobj)

#CCL
genes <- str_sort( grep(allgenes, pattern = 'CCL', value = T, ignore.case = T), numeric = T) 
DoHeatmap(sobj, genes)

#CXCL
genes <- str_sort( grep(allgenes, pattern = 'CXCL', value = T, ignore.case = T), numeric = T) 
DoHeatmap(sobj, genes)

# IL
genes <- str_sort( grep(allgenes, pattern = '^IL', value = T, ignore.case = T), numeric = T) 
genes <- genes[!grepl(genes, pattern = '^ILR', ignore.case = T)]
DoHeatmap(sobj, genes)


### cytokine list from cytreg
cts <- read.csv('~/Dropbox/data/general_data_utilities/cytreg/cytreg_v2.csv')
hom <- read.csv('~/Dropbox/data/general_data_utilities/biomart/biomart_nodups_may05-2022_feb2021archive.csv')

#keep only mouse and human....
cts <- cts[cts$species == 'H' | cts$species == 'M', ]

#just get cytokines...
cyto <- str_sort(unique(cts$cytokine), numeric = T)
cyto[!(cyto %in% hom$HGNC.symbol)] #there are some important ones missing... esp IL17, others...

# some of them may be exxcluded because there are PARALOGS in mouse (ie il17 )

#get only ones with homologs...
cyto <- cyto[cyto %in% hom$HGNC.symbol]
hom <- hom[match(cyto, hom$HGNC.symbol),]
cyto <- hom$MGI.symbol

DoHeatmap(sobj, cyto)


write.csv('finalized_m0m1m2_markers.csv', x = m, row.names = F, quote = T)



rm(list=ls())






















### try to see TKO module scores of m0, m1, m2



setwd('~/Dropbox/data/deyou/macspectrum/')


m <- read.csv('finalized_m0m1m2_markers.csv')








### read in tko dko tpm

macspectrumdir <- 'results/comparative-de/TKO-vs-DKO/4.downstream/macspectrum'
macspectrumdir <- paste0('~/Dropbox/data/bangdata/2021october-TKOvsDKO/', macspectrumdir)




mat2 <- read.csv(paste0(macspectrumdir, '/tpm.csv'))
f2 <- read.csv(paste0(macspectrumdir, '/feature.csv'))




# set up module score function
set.seed(2022)
modulescore <- function(gem, genelist, numbins=NULL, numcontrolgenesperbin=NULL){
  if(is.null(numbins)){numbins=24}
  if(is.null(numcontrolgenesperbin)){numcontrolgenesperbin=100}
  
  testgenes <- genelist
  controlgenes <- rownames(gem)
  controlgenes <- controlgenes[!(controlgenes %in% testgenes)]
  
  #testmeans <- Matrix::rowMeans(gem[rownames(gem) %in% genes,])
  
  ### bin the genes by avg expression ###
  
  #get average gene expression
  avgs <- Matrix::rowMeans(gem)
  avgs <- avgs[order(avgs)]
  
  #cut; get the gene bins
  bins <- cut_number(avgs, n=numbins)
  
  bindf <- data.frame(genenames = names(avgs), 
                      bins = bins)
  
  #select control genes from same expression bins
  controlgenes <- c()
  for(genedex in 1:length(testgenes)){
    
    #get gene and bin
    gene <- testgenes[genedex]
    genebin <- bindf[bindf$genenames == gene,2]
    
    #select all genes from that bin, to sample from
    tmpcontrols <- bindf[bindf$bins == genebin,]
    
    #exclude the actual test genes
    tmpcontrols <- tmpcontrols[!(tmpcontrols$genenames %in% testgenes),]
    
    #if num controls exceeds bin size, take all genes in bin...
    numtotake <- ifelse(nrow(tmpcontrols) < numcontrolgenesperbin,
                        yes = nrow(tmpcontrols),
                        no = numcontrolgenesperbin)
    
    #sample the genes
    tmpcontrols <- sample(tmpcontrols$genenames, size = numtotake, replace = F)
    
    controlgenes <- unique(c(controlgenes, 
                             tmpcontrols
    ))
  }
  
  
  
  #get control gene mean expression for each sample
  samplemeans_controls <- Matrix::colMeans( gem[rownames(gem) %in% controlgenes,] )
  
  #get test genes mean expression for each samoke
  samplemeans_test <- Matrix::colMeans( gem[rownames(gem) %in% testgenes,] )
  
  #subtract them to get module score
  modulescore <- samplemeans_test - samplemeans_controls
  
  return(modulescore)
}




### m0, m1, m2 modulle scores
# translate ensembl back to gencode...
gencode <- read.csv('~/Dropbox/data/bangdata/2021october-TKOvsDKO/data/gencode.vM23.ids_names_types.csv')

#fix ensembl...
gencode$gene_id <- str_split_fixed(gencode$gene_id, '\\.', 2)[,1]

#keep same
# remove duplicate genes, so annoying
dupgenes <- names(table(gencode$gene_name)[table(gencode$gene_name)>1])
gencode <- gencode[!(gencode$gene_name %in% dupgenes),]
# gencode <- gencode[gencode$gene_name %in% m$gene,]
gencode <- gencode[gencode$gene_id %in% mat2$Ensembl_ID,]
mat2 <- mat2[mat2$Ensembl_ID %in% gencode$gene_id,]

#match order
gencode <- gencode[match(mat2$Ensembl_ID, gencode$gene_id),]

#reformat, gene symbols as rownames
gem <- mat2[,-1]
rownames(gem) <- gencode$gene_name


mod <- data.frame( m0 = modulescore(gem, m[m$cluster=='m0', "gene"]),
                   m1 = modulescore(gem, m[m$cluster=='m1', "gene"]),
                   m2 = modulescore(gem, m[m$cluster=='m2', "gene"])
)


mod$genotype <- c(rep('TKO',3), rep('DKO',3))
modx <- reshape2::melt(mod)
ggplot(modx, aes(variable, value))+
  geom_point() + geom_boxplot() +
  facet_wrap(~genotype)


#eclude DKO
mod <- mod[mod$genotype=='TKO',]
modx <- modx[modx$genotype=='TKO',]
ggplot(modx, aes(variable, value))+
  geom_boxplot() + geom_point() +
  theme_light()





### try pathway analysis wwith ORA test??

res <- read.csv('~/Dropbox/data/bangdata/2021october-TKOvsDKO/results/comparative-de/TKO-vs-DKO/2.deresults/DEresults-normalized_count_matrix.csv')

#remove upds...
dupgenes <- names(table(res$mgi_symbol)[table(res$mgi_symbol)>1])
res <- res[!(res$mgi_symbol %in% dupgenes),]

#select significant genes...
sig <- res[res$padj < 0.05,]
sig <- sig[sig$log2FoldChange > 1,]


#rank by FC
genes <- sig$log2FoldChange
names(genes) <- sig$mgi_symbol

#sort them by decreasing magnitude
genes <- genes[order(abs(genes), decreasing = T)]

# just in case, remove any NAs...
genes <- na.omit(genes)


#just give the actual gene names as input...
genenames <- names(genes)



#prep term2gene...
term2gene <- m[,c('cluster', 'gene')]



## run the enrichment analysis ##
library(clusterProfiler)
up_msigdb_ora <- enricher(genenames,
                          TERM2GENE = term2gene, 
                          # TERM2NAME = term2name,
                          # OrgDb = org.Mm.eg.db,
                          # keyType = 'ENSEMBL',
                          pvalueCutoff = 1, qvalueCutoff = 1, maxGSSize = 1000
)


#run termsimilarity
up_msigdb_ora <- enrichplot::pairwise_termsim(up_msigdb_ora)

#emaapplot, with just the top 30...
up_emap_pval <- emapplot(up_msigdb_ora, layout='graphopt', repel=T)


#try to plot color = num genes in overexp gene list?
up_msigdb_ora@result$Percent_of_DEGs <- (up_msigdb_ora@result$Count / length(genenames))*100



up_emap_perc <- emapplot(up_msigdb_ora, color='Percent_of_DEGs', layout='graphopt', repel=T,
                         coords = up_emap_pval$data[,1:2])+
  scale_fill_distiller(palette = 'Reds', direction = 1, name='Percent of DEGs')



### onlly m0 is significantly enriched in TKO upreg genes
dotplot(up_msigdb_ora, x = 'Count')


# TKO is high in m0, m2
ComplexHeatmap::Heatmap(mod[1:3,1:3])





### just do it myself, thsi package sucks
term2gene <- m[,c('cluster', 'gene')]

#get universe genes
univ <- res$mgi_symbol

#keep only pway gnees in the universe?
term2gene <- term2gene[term2gene$gene %in% res$mgi_symbol,]
pwaylens <- table(term2gene$cluster)

#get my degs
mydegs <- sig$mgi_symbol


freslist <- list()
for(pway in names(pwaylens)){
  
  
  #set up contingency table
  
  #         # DE genes    # non DE genes
  
  # in Pway
  # not in Pway
  
  # it should all add up to nrow(res)
  
  #get pway genes
  pwaygenes <- term2gene[term2gene$cluster==pway,"gene"]
  
  # find DE genes in pathway / not in pathway
  de_in_pway <- table(mydegs %in% pwaygenes)
  
  #remove DEGs from universe
  univ <- univ[!(univ %in% mydegs)]
  
  #find nonDE genes in pathwya / not in pathway
  pwaygenes_nonDE <- pwaygenes[!(pwaygenes %in% mydegs)]
  nonde_in_pway <- table(univ %in% pwaygenes_nonDE)
  
  
  #make sure the vectors defined above via table have both false and true...
  # if not, need to set to 0...
  if( !('TRUE' %in% names(de_in_pway)) ){ de_in_pway['TRUE'] <- 0 }
  if( !('FALSE' %in% names(de_in_pway)) ){ de_in_pway['FALSE'] <- 0 }
  
  if( !('TRUE' %in% names(nonde_in_pway)) ){ nonde_in_pway['TRUE'] <- 0 }
  if( !('FALSE' %in% names(nonde_in_pway)) ){ nonde_in_pway['FALSE'] <- 0 }
  
  
  
  #contingency table set up...
  ct <- data.frame(DE =     c(de_in_pway['TRUE'], de_in_pway['FALSE']),
                   Not_DE = c(nonde_in_pway['TRUE'], nonde_in_pway['FALSE']),
                   row.names = c('In_Pway', 'Not_in_pway')) 
  
  
  f <- fisher.test(ct, alternative = 'greater')
  
  freslist[[pway]] <- data.frame(ID = pway, OR = f$estimate, P = f$p.value, row.names = NULL, pwaysize = length(pwaygenes), overlapsize = ct[1,1])
  
  
  
}

fres <- dplyr::bind_rows(freslist)


fres$ID <- c('M0', 'M1', 'M2')


enrichment <- ggplot(fres, aes(ID, OR, col = -log10(P), size = overlapsize))+ 
  geom_point() + scale_color_gradient(low = 'grey', high = 'purple')+
  theme_light()


pdf("/Users/ferrenaa/Dropbox/data/bangdata/2021october-TKOvsDKO/results/comparative-de/TKO-vs-DKO/4.downstream/mMCP_counter/tkoenrichment_m0m1m2.pdf",
    height = 4, width = 4)

enrichment

dev.off()



#### try GSVA
library(GSVA)


#input is logTPM, need to specify gaussian kernel






setwd('~/Dropbox/data/deyou/macspectrum/')


m <- read.csv('finalized_m0m1m2_markers.csv')








### read in tko dko tpm

macspectrumdir <- 'results/comparative-de/TKO-vs-DKO/4.downstream/macspectrum'
macspectrumdir <- paste0('~/Dropbox/data/bangdata/2021october-TKOvsDKO/', macspectrumdir)




mat2 <- read.csv(paste0(macspectrumdir, '/tpm.csv'))
f2 <- read.csv(paste0(macspectrumdir, '/feature.csv'))


### m0, m1, m2 modulle scores
# translate ensembl back to gencode...
gencode <- read.csv('~/Dropbox/data/bangdata/2021october-TKOvsDKO/data/gencode.vM23.ids_names_types.csv')

#fix ensembl...
gencode$gene_id <- str_split_fixed(gencode$gene_id, '\\.', 2)[,1]

#keep same
# remove duplicate genes, so annoying
dupgenes <- names(table(gencode$gene_name)[table(gencode$gene_name)>1])
gencode <- gencode[!(gencode$gene_name %in% dupgenes),]
# gencode <- gencode[gencode$gene_name %in% m$gene,]
gencode <- gencode[gencode$gene_id %in% mat2$Ensembl_ID,]
mat2 <- mat2[mat2$Ensembl_ID %in% gencode$gene_id,]

#match order
gencode <- gencode[match(mat2$Ensembl_ID, gencode$gene_id),]

#reformat, gene symbols as rownames
gem <- as.matrix(mat2[,-1])
rownames(gem) <- gencode$gene_name


#exclude dko...
gem <- gem[,1:3]

#log it...
gem <- log(gem+1)


#make gene lists
gs <- list(m0 = m[m$cluster=='m0','gene'],
           m1 = m[m$cluster=='m1','gene'],
           m2 = m[m$cluster=='m2','gene'])

gsvares <- gsva(as.matrix(gem), gs, kcdf='Gaussian' )


modx <- reshape2::melt(gsvares)
ggplot(modx, aes(Var1, value))+
  geom_boxplot() + geom_point() +
  theme_light()







#### update Dec 9 for supp figure reanalysis of macspectrum

library(tidyverse)
library(RISC)
library(Seurat)
library(FerrenaSCRNAseq)


set.seed(2022)


setwd('~/Dropbox/data/deyou/macspectrum/')

risc <- readRDS('RISC_int.rds')

markers <- read.csv('finalized_m0m1m2_markers.csv')

n <- 20
top <- markers %>% group_by(cluster) %>% top_n(n = n, wt = avg_log2FC)



dir.create('supp_plots/')

#pdf()
drisc <- RISC::DimPlot(risc, colFactor = 'orig.ident')



hm <- RISC::Heat(risc, 
                 colFactor = 'orig.ident', 
                 genes = top$gene, 
                 gene.lab = T, 
                 smooth = F,
                 cell.lab.size = 8
)

hm_all <- RISC::Heat(risc, 
                     colFactor = 'orig.ident', 
                     genes = markers$gene, 
                     gene.lab = F, 
                     smooth = F,
                     gene.cluster = 3
)

data0 <- risc
# get the metadata
md <- data0@coldata

# get the gene expression matrix
gem <- data0@assay$logcount
gem <- do.call(cbind, gem)

#get UMAP coords
umap <- data0@DimReduction$cell.umap

#make seurat object, for nice pre-made plotting
sobj <- CreateSeuratObject(counts = gem, meta.data = md)
sobj[['umap']] <- Seurat::CreateDimReducObject(umap, key = 'umap')


rm(data0, md, umap,gem)


for(module in unique(markers$cluster) ){
  message('\n\n', module)
  genes <- list(module = markers[markers$cluster == module,'gene'])
  sobj <- AddModuleScore(sobj, features = genes, name = module)
  colnames(sobj@meta.data)[ncol(sobj@meta.data)] <- module
}

vlns <- VlnPlot(sobj, c('m0', 'm1', 'm2'), group.by = 'orig.ident', combine = F, pt.size = 0.1)
vlns <- lapply(vlns, function(gg){gg+
    scale_y_continuous(limits = c(-1,3))+
    NoLegend()
})

vlns <- patchwork::wrap_plots(vlns, ncol = 3)




pdf('supp_plots/dimplot.pdf', 4,4)
drisc
dev.off()

pdf('supp_plots/heatmaps.pdf', 6,6)
hm
grid::grid.newpage()
hm_all
dev.off()

pdf('supp_plots/vlns.pdf', height = 4,6)
vlns
dev.off()

