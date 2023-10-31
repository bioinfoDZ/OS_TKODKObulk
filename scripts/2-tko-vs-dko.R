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


### input parameters ###

#set conditions to compare
# c1 = factor of interest
# c2 = control to be compared against

c1 = "TKO"
c2 = "DKO"


#absolute l2fc and p-val thtesholds, set a priori
lfc_thres <- 1
pval_thres <- 0.05


### make the output folders ###

#get title of comparison
comptitle <- paste0(c1, '-vs-', c2)

#make whole results dir
compoutdir <- paste0('results/comparative-de/', comptitle)
dir.create(compoutdir)


#make qc outdir
dir.create(   paste0(compoutdir, '/0.qc/')  )
dir.create(   paste0(compoutdir, '/1.pca/')   )
dir.create(   paste0(compoutdir, '/2.deresults/')   )
dir.create(   paste0(compoutdir, '/3.gsea')   )





### QC & PCA ###

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



#get minimal md and sample names
condmd <- md[md$Condition %in% c(c1, c2),]
condsamps <- condmd$Sample



### filtering ###

#remove empty genes
gem <- gem[rowSums(gem) > 0,]

#filter lowly exp (under 10 counts) genes
#dim(gem[rowSums(gem) >= 10,])
gem <- gem[rowSums(gem) >= 10,]
#dim(gem)




# 
# 
# ### filtering outlier genes ###
# 
# # we have observed a strange phenomenon in which each individual replicate
# # expresses certain outlier genes at a very high level.
# # Thus we will try to define and filter out outlier genes.
# # Outlier genes will be defined as genes with 5x higher expression
# # than ALL other samples in the same genotype.
# # if the gene has this characteristic, it is removed from both genotypes
# 
# 
# 
# 
# #normalize with DESeq2
# gemx <- t( t(gem) / DESeq2::estimateSizeFactorsForMatrix(gem) )
# 
# 
# 
# badgenes <- c()
# for(genotype in unique(condmd$Condition) ){
#   
#   #get the gem for each genotype
#   mdtmp <- condmd[condmd$Condition == genotype,]
#   tmpgem <- gemx[,colnames(gemx) %in% mdtmp$Sample, drop=F]
#   
#   
#   message('\nStarting ', genotype)
#   
#   
#   #progress bar
#   total = nrow(tmpgem) - 1
#   pb <- txtProgressBar(min = 0, max = total, style = 3)
#   
#   
#   genotype_badgenes <- c()
#   for(genedex in c(1:nrow(tmpgem)) ){
#     
#     #get genename
#     genename <- rownames(tmpgem)[genedex]
#     
#     #get gem for genotype
#     generow <- tmpgem[genedex,]
#     
#     ## test if the max is higher than 5x any of the others... ##
#     #get the max gene
#     maxsamp <- max(generow)
#     
#     #test if ALL of the others are hgiher than maxsamp / 5
#     maxdivfive <- maxsamp / 5
#     
#     # actual test:
#     outliergene <- ifelse(test = all(test = generow[-which.max(generow)] < maxdivfive),
#                           yes = T, no = F)
#     
#     
#     if(outliergene){
#       genotype_badgenes <- c(genotype_badgenes, genename)
#     }
#     
#     
#     setTxtProgressBar(pb, genedex)
#     
#   } #end gem sweep
#   
#   
#   
#   badgenes <- c(badgenes, genotype_badgenes)
#   
# }
# 
# badgenes <- unique(badgenes)
# 
# #remove these from the gem and test...
# 
# gem <- gem[!(rownames(gem) %in% badgenes),]
# 
# rm(gemx, badgenes, tmpgem, mdtmp, pb, genotype_badgenes)
# rm(genedex,genename,generow,maxdivfive,maxsamp, outliergene, total, genotype)
# 
# 
# 
# 



##### PCA #####
#Novogene PCA showed one sample (DJ208F0) as a major outlier...

##use DESEQ2: size factor norm, and rlog ##

#create dds obj
dds <- DESeqDataSetFromMatrix(gem, md,
                              design = ~ Condition)


#run DESeq2 
dds <- DESeq(dds)





#for plotting only, see how size factor norm affect lib size


gemnorm <- counts(dds, normalized=T)


#### plot the lib size #####
pdf <- data.frame(samp = colnames(gemnorm), 
                  numreadsaligned = colSums(gemnorm),
                  condition = md$Condition,
                  batch = md$Batch,
                  color = md$Color,
                  stringsAsFactors = F
)


#subset to include only the samples of interest
pdf <- pdf[pdf$condition %in% c(c1,c2),]

#order them from hi to low
pdf$samp <- factor(pdf$samp, levels = pdf[order(pdf$numreadsaligned, decreasing = T),"samp"])
pdf$condition <- factor(pdf$condition, levels = unique(pdf[order(pdf$numreadsaligned, decreasing = T),"condition"]))
cols <- unique(pdf[order(pdf$numreadsaligned, decreasing = T),"color"])


libsize_norm <- ggplot(pdf, aes(x = samp, y = numreadsaligned, fill = condition))+
  geom_bar(stat = 'identity')+
  scale_fill_brewer(palette = 'Set1',direction = -1)+
  scale_y_continuous(labels = scales::comma)+
  theme_light()+
  #scale_fill_brewer(palette = 'Set2')+
  scale_fill_manual(values = cols)+
  scale_y_continuous(limits = c(0,25000000), labels = scales::comma)+
  coord_flip()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  labs(title = 'Number of reads aligned to protein-coding genes', 
       subtitle = 'With DEseq2 size factor normalization',
       y = 'Number of reads aligned', x = 'Sample')

#libsize_norm

ggsave(libsize_norm, 
       filename = paste0(compoutdir, '/0.qc/', 'libsize-postnorm.jpg'), 
       height = 5, width = 5, dpi = 300)







#rlog
rlog <- rlog(dds, blind = F)

#run PCA




## recode pca myself... ##

#get the regularized log transformed data
rlmat <- assay(rlog)

#subet to only get the key samples
rlmat <- rlmat[,colnames(rlmat) %in% condsamps]


#get the top row variance genes, as performed by DESeq2 plotPCA function
rvs <- rowVars(rlmat) ; names(rvs) <- rownames(rlmat)
rlmat <- rlmat[order(rvs, decreasing = T),]

#set the number of variable genes
ntop = 2000



jpeg(filename = paste0(compoutdir, '/1.pca/', 'gene-variance.jpg'), 
     height = 7, width = 5, res = 300, units = 'in')

#plot num genes
par(mfrow=c(2,1))

#see how ntop selects genes
plot(sort(rvs, decreasing = T), ylab ='Gene Variance', main = paste0('Selecting top ', ntop, ' genes')) ; abline(v = ntop)

#see how ntop selects genes, with log y
plot(sort(log(rvs), decreasing = T), ylab = ('Log Gene Variance')) ; abline(v = ntop)

dev.off()

#subset rlmat by ntop
rlmat <- rlmat[1:ntop,]

#instead of scaling, subtract by mean
# https://www.biostars.org/p/387863/
rlmat <- rlmat - rowMeans(rlmat)

#transpose the mat
rlmat <- t(rlmat)


# perform a PCA on the data in assay(x) for the selected genes
pcaobj <- prcomp( rlmat )

pdf <- as.data.frame(pcaobj$x[,1:2])
pdf$name <- condmd$Sample
pdf$batch <- condmd$Batch
pdf$age <- condmd$`Age_at_sack (week)`
pdf$tumor_location <- condmd$Tumor_location
pdf$condition <- factor(condmd$Condition)
pdf$color <- factor(condmd$Color, levels = unique(condmd$Color)[order(unique(pdf$condition))] )



#data.frame(levels(pdf$condition), levels(pdf$color))


pca_recoded <- ggplot(pdf, aes(PC1, PC2, col = condition, shape = batch))+
  geom_point(size = 3)+
  ggrepel::geom_text_repel(aes(label = name), box.padding = 0.4)+
  #scale_color_brewer(palette = 'Set2', direction = -1)+
  scale_color_manual(values = levels(pdf$color))+
  theme_light()

#pca_recoded

ggsave(pca_recoded, 
       filename = paste0(compoutdir, '/1.pca/', 'pca.jpg'), 
       height = 5, width = 5, dpi = 300)


#try to cluster with the HVGs
#HC on correlation matrix instead of distance

# Pairwise correlation between samples (columns)
cols.cor <- cor(t(rlmat), use = "pairwise.complete.obs", method = "pearson")
cols.cor <- as.dist(1-cols.cor)

col_dend <- dendsort(as.dendrogram(hclust(cols.cor)))

colors <- rev(colorRampPalette( brewer.pal(9, "RdYlBu"))(255))

colors <- circlize::colorRamp2(seq(min(-1), max(1), length = length(colors)), colors)

cormat <- cor(t(rlmat))



annotdf <- data.frame(sample = colnames(cormat),
                      sampcheck = md[match(colnames(cormat), md$Sample), "Sample"],
                      cond = md[match(colnames(cormat), md$Sample), "Condition"],
                      color = md[match(colnames(cormat), md$Sample), "Color"])

hacol <- list(Condition = annotdf$Color) ; names(hacol[[1]]) <- annotdf$Condition
ha <- HeatmapAnnotation(Condition = annotdf[,3], col = hacol, name = 'Condition')
hal <- rowAnnotation(Condition = annotdf[,3], col = hacol, name = 'Condition')

hm <- Heatmap(cormat, name = "Pearson r", 
              col=colors,
              cluster_rows = col_dend, cluster_columns = col_dend, 
              row_names_side = 'left',
              column_names_side = 'top',
              column_dend_height = unit(5, "cm"),
              row_dend_width =  unit(5, "cm"),
              #width = unit(12, "cm"), height = unit(12, "cm"),
              top_annotation = ha, left_annotation = hal
)


#hm

jpeg(paste0(compoutdir, '/1.pca/', 'hierarchicalclustering.jpg'),
     height = 8, width = 10, units = 'in', res = 300)

print(hm)

dev.off()





### comparative DE analysis ###

#get results:

res <- results(dds, contrast = c('Condition', c1, c2))


#set as DF and order by Pvalue
res <- as.data.frame(res)
res <- res[order(res$padj),]

#remove outlier genes with p = NA (cooks distance)
res <- res[!is.na(res$pvalue),]
res <- res[!is.na(res$padj),]





#make the volcanoplot

volcano <- FerrenaBulkRNAseq::volcanoplot(results = res, lfc_thres = lfc_thres,
                                          colors = c('firebrick', 'steelblue'),
                                          change_gene_label=T,
                                          gene_label_equivalency=gencode)







ggsave(volcano, 
       filename = paste0(compoutdir, '/2.deresults/', 'volcano.jpg'), 
       height = 5, width = 5, dpi = 300)








### heatmap of top DEGs ###

#get the norm counts matrix and subset to keep genes selected
gemnorm <- counts(dds, normalized=T)
gemnorm <- gemnorm[,colnames(gemnorm) %in% condmd$Sample]




#select top 30 genes by p-value in both directions
tmpres <- res[res$pvalue < 0.05,]
tmprespos <- rownames(tmpres[tmpres$log2FoldChange > 0, ])[1:30]
tmpresneg <- rownames(tmpres[tmpres$log2FoldChange < 0, ])[1:30]

genes <- c(tmprespos, rev(tmpresneg))





#set title, will be used as "column title" but goes on top so looks fine
hmtitle = paste0('top ', length(genes), ' DEGs by pvalue\n',
                 length(tmprespos), ' overexp and ', length(tmpresneg), ' underexp')


hm_top60 <- FerrenaBulkRNAseq::heatmapplot(gemnorm, genes, condmd, 
                                           heatmaptitle = hmtitle,
                                           change_gene_label=T,
                                           gene_label_equivalency=gencode
)


jpeg(paste0(compoutdir, '/2.deresults/', 'heatmap-top-60-degs-pvalue.jpg'),
     height = 8, width = 10, units = 'in', res = 300)

print(hm_top60)

dev.off()







### heatmap of top DEGs, excluding extreme LFC genes ###

#select top 30 genes by p-value in both directions, excluding extreme LFC genes
tmpres <- res[abs(res$log2FoldChange) <= 15,]
tmpres <- tmpres[tmpres$pvalue < 0.05,]
tmprespos <- rownames(tmpres[tmpres$log2FoldChange > 0, ])[1:30]
tmpresneg <- rownames(tmpres[tmpres$log2FoldChange < 0, ])[1:30]

genes <- c(tmprespos, rev(tmpresneg))




#set title, will be used as "column title" but goes on top so looks fine
hmtitle = paste0('top ', length(genes), ' DEGs by pvalue\n', 
                 'after excluding genes with LFC +/- ', 15)

hm_noext <- FerrenaBulkRNAseq::heatmapplot(gemnorm, genes, condmd, 
                                           heatmaptitle = hmtitle,
                                           change_gene_label=T,
                                           gene_label_equivalency=gencode)


jpeg(paste0(compoutdir, '/2.deresults/', 'heatmap-exclude-outlier-genes.jpg'),
     height = 8, width = 10, units = 'in', res = 300)

print(hm_noext)

dev.off()






#make heatmap of all DEGs, defined as pvalue < 0.05

#select significant genes by p-value in both directions
tmpres <- res[res$padj < pval_thres,]
tmpres <- tmpres[abs(tmpres$log2FoldChange) > lfc_thres,]

genes <- rownames(tmpres)



#set title, will be used as "column title" but goes on top so looks fine
hmtitle = paste0('all DEGs with padj < ', pval_thres,' and LFC > +/- ', lfc_thres, '\n',
                 table(sign(tmpres$log2FoldChange))[2], ' overexp and ', table(sign(tmpres$log2FoldChange))[1], ' underexp')




hm_alldegs_pval <- FerrenaBulkRNAseq::heatmapplot(gemnorm, genes, condmd, 
                                                  heatmaptitle = hmtitle,
                                                  change_gene_label=T,
                                                  gene_label_equivalency=gencode)


jpeg(paste0(compoutdir, '/2.deresults/', 'heatmap-alldegs-pval.jpg'),
     height = 8, width = 10, units = 'in', res = 300)

print(hm_alldegs_pval)

dev.off()











#make heatmap of DKOAA DEGs

#select DEGs from DKOAA results
dkoaares <- read.csv('/Users/ferrenaa/Dropbox/data/bangdata/dko_dkoaa_october2020/start_to_end/results_figures/DEG_filtered_padj0.05.csv')
dkoaadegs <- dkoaares[dkoaares$padj < 0.05,]

genes <- dkoaadegs$ensembl_gene_id



#set title, will be used as "column title" but goes on top so looks fine
hmtitle = paste0('DKOAA vs DKO genes with padj < 0.05')




hm_dkoaa <- FerrenaBulkRNAseq::heatmapplot(gemnorm, genes, condmd, 
                                           heatmaptitle = hmtitle,
                                           change_gene_label=T,
                                           gene_label_equivalency=gencode)




jpeg(paste0(compoutdir, '/2.deresults/', 'heatmap-dkoaa-degs.jpg'),
     height = 8, width = 10, units = 'in', res = 300)

print(hm_dkoaa)

dev.off()








### write out the results ###

#get normalized counts matrix
save <- as.data.frame(gemnorm)
save <- save[match(rownames(res), rownames(save)),]


#add ensemble ID and DEseq2 results
gcsave <- gencode[match(rownames(save), gencode$gene_id),]
save <- cbind(gcsave$id_with_version, gcsave$gene_name, res, save)
colnames(save)[1:2] <- c('ensembl_gene_id', "mgi_symbol")

write.csv(file = paste0(compoutdir, '/2.deresults/', 'DEresults-normalized_count_matrix.csv'),
          x = save, row.names = F, quote = F)


rm(gem, pcaobj, rlmat, rlog, dkoaares, dkoaadegs)








### GSEA ###


pathways <- msigdbr(species = 'Mus musculus')

cats <- unique(pathways$gs_cat)
cats <- str_sort(cats, numeric = T)

#for each msigdb cateogory, check if there are subcategories
# if not, analyze the whole category
# if yes, make a sub-dir and analyse them separately


for(cat in cats){
  cp <- pathways[pathways$gs_cat == cat,]
  
  message('\nRunning ', cat)
  
  catdir <-  paste0(compoutdir, '/3.gsea/', cat)
  dir.create(catdir)
  
  # if the category has no sub-categories:
  if( all(unique(cp$gs_subcat) == '') ){
    
    cp  <- cp %>% split(x = .$ensembl_gene, f = .$gs_name)
    
    
    #run GSEA and make main dotplot
    gseares <- FerrenaBulkRNAseq::gsea.results(save, pathways = cp)
    gseadotplot <- FerrenaBulkRNAseq::gsea.dotplot.onecol(gseares = gseares[[1]], pathwayfontsize = 7) +
      ggtitle(cat)
    
    #save gsearesults
    write.csv(file = paste0(catdir, '/fgsea-results.csv'),
              x = gseares[[1]], row.names = F, quote = F)
    
    #save dotplot
    ggsave(gseadotplot, 
           filename = paste0(catdir, '/gseadotpot.jpg'), 
           height = 10, width = 10, bg = 'white')
    
    
    
    
    #save leading edge
    le <- gseares[[2]]
    
    #make a padded DF from a list of char vectors of different lengths:
    # https://stackoverflow.com/questions/43159384/convert-a-list-of-character-vectors-into-dataframe
    le <- data.frame(lapply(le, `length<-`, max(lengths(le))))
    
    
    write.csv(file = paste0(catdir, '/leadingedge.csv'),
              x = le, row.names = F, quote = F)
    
    
    
    #loop through significant pathways, make leading edge plots
    gsearesdf <- gseares[[1]]
    gsearesdf <- gsearesdf[!is.na(gsearesdf$pval),]
    sigpathways <- gsearesdf[gsearesdf$padj < 0.05,]
    
    
    dir.create(paste0(catdir, '/heatmap-leadingedge'))
    dir.create(paste0(catdir, '/heatmap-wholepathway'))
    
    
    #for each significant pathway, 
    # make a heatmap using both leading edge and whole pathway genes
    for(spname in sigpathways$pathway){
      
      #get the pathway
      sp <- sigpathways[sigpathways$pathway == spname,]
      
      #get the leading edge genes
      sple <- le[[spname]]
      
      #get rid of NAs
      sple <- sple[!is.na(sple)]
      
      #if not enough genes causes problems...
      if(length(sple) < 3){
        next
      }
      
      
      #take only top 75 genes, do we really need to show more? lol
      sple <- sple[1:75]
      
      #get rid of NAs
      sple <- sple[!is.na(sple)]
      
      #use MGI names
      sple_mgi <- gencode[match(sple, gencode$gene_id),]
      
      # #if any are duplicated, use ensemble, else use MGI
      # if(!any(duplicated(sple_mgi$mgi_symbol))){
      #   sple <- sple_mgi$mgi_symbol
      # }
      
      #make a heatmap title
      hmtitle <- paste0(spname, '\n',
                        'NES=', round(sp$NES, 2), ', FDR=', round(sp$padj, 2))
      
      
      
      jpeg(paste0(catdir, '/heatmap-leadingedge/', spname, '.jpg'),
           height = 8, width = 10, units = 'in', res = 300)
      
      
      print(FerrenaBulkRNAseq::heatmapplot(gemnorm, sple, condmd[,1:10], 
                                           heatmaptitle = hmtitle,
                                           change_gene_label=T,
                                           gene_label_equivalency=gencode,
                                           column_km = 2))
      
      dev.off()
      
      #all genes in pathway
      
      allpathwaygenes <- cp[[spname]]
      allpathwaygenes <- allpathwaygenes[allpathwaygenes %in% rownames(gemnorm)]
      
      
      jpeg(paste0(catdir, '/heatmap-wholepathway/', spname, '.jpg'),
           height = 8, width = 10, units = 'in', res = 300)
      
      print(FerrenaBulkRNAseq::heatmapplot(gemnorm, allpathwaygenes, condmd[,1:10], 
                                           heatmaptitle = hmtitle,
                                           change_gene_label=T,
                                           gene_label_equivalency=gencode,
                                           column_km = 2))
      
      dev.off()
      
      
    } #finish heatmap loop
    
    
    
    
  } # finish no subcat if statement
  
  
  #if the category does have sub-categories
  else{
    #loop thru each subcategory and do gsea + heatmaps
    subcats <- unique(cp$gs_subcat)
    subcats <- subcats[!grepl('Legacy', subcats, ignore.case = T)]
    
    for(subcat in subcats){
      
      message('- ', subcat)
      #make outdir of this sub-category
      subcatdir <- paste0(catdir, '/', subcat)
      dir.create(subcatdir)
      
      #get the sub-cateogyr pathways and genes and format as list for fgsea
      cpsub <- cp[cp$gs_subcat == subcat,]
      cpsub  <- cpsub %>% split(x = .$ensembl_gene, f = .$gs_name)
      
      
      #run GSEA and make main dotplot
      gseares <- FerrenaBulkRNAseq::gsea.results(save, pathways = cpsub)
      gseadotplot <- FerrenaBulkRNAseq::gsea.dotplot.onecol(gseares = gseares[[1]], pathwayfontsize = 7) +
        ggtitle(subcat) 
      
      #save gsearesults
      write.csv(file = paste0(subcatdir, '/fgsea-results.csv'),
                x = gseares[[1]], row.names = F, quote = F)
      
      #save dotplot
      ggsave(gseadotplot, 
             filename = paste0(subcatdir, '/gseadotpot.jpg'), 
             height = 10, width = 10, bg = 'white')
      
      
      
      
      #save leading edge
      le <- gseares[[2]]
      
      #make a padded DF from a list of char vectors of different lengths:
      # https://stackoverflow.com/questions/43159384/convert-a-list-of-character-vectors-into-dataframe
      le <- data.frame(lapply(le, `length<-`, max(lengths(le))))
      
      
      write.csv(file = paste0(subcatdir, '/leadingedge.csv'),
                x = le, row.names = F, quote = F)
      
      
      
      #loop through significant pathways, make leading edge plots
      gsearesdf <- gseares[[1]]
      gsearesdf <- gsearesdf[!is.na(gsearesdf$pval),]
      sigpathways <- gsearesdf[gsearesdf$padj < 0.05,]
      
      
      dir.create(paste0(subcatdir, '/heatmap-leadingedge'))
      dir.create(paste0(subcatdir, '/heatmap-wholepathway'))
      
      
      #for each significant pathway, 
      # make a heatmap using both leading edge and whole pathway genes
      for(spname in sigpathways$pathway){
        
        #get the pathway
        sp <- sigpathways[sigpathways$pathway == spname,]
        
        #get the leading edge genes
        sple <- le[[spname]]
        
        #get rid of NAs
        sple <- sple[!is.na(sple)]
        
        #if not enough genes causes problems...
        if(length(sple) < 3){
          next
        }
        
        
        #take only top 75 genes, do we really need to show more? lol
        sple <- sple[1:75]
        
        #get rid of NAs
        sple <- sple[!is.na(sple)]
        
        #use MGI names
        sple_mgi <- gencode[match(sple, gencode$gene_id),]
        
        # #if any are duplicated, use ensemble, else use MGI
        # if(!any(duplicated(sple_mgi$mgi_symbol))){
        #   sple <- sple_mgi$mgi_symbol
        # }
        
        #make a heatmap title
        hmtitle <- paste0(spname, '\n',
                          'NES=', round(sp$NES, 2), ', FDR=', round(sp$padj, 2))
        
        
        
        jpeg(paste0(subcatdir, '/heatmap-leadingedge/', spname, '.jpg'),
             height = 8, width = 10, units = 'in', res = 300)
        
        
        print(FerrenaBulkRNAseq::heatmapplot(gemnorm, sple, condmd[,1:10], 
                                             heatmaptitle = hmtitle,
                                             change_gene_label=T,
                                             gene_label_equivalency=gencode,
                                             column_km = 2))
        
        dev.off()
        
        #all genes in pathway
        
        allpathwaygenes <- cpsub[[spname]]
        allpathwaygenes <- allpathwaygenes[allpathwaygenes %in% rownames(gemnorm)]
        
        
        jpeg(paste0(subcatdir, '/heatmap-wholepathway/', spname, '.jpg'),
             height = 8, width = 10, units = 'in', res = 300)
        
        print(FerrenaBulkRNAseq::heatmapplot(gemnorm, allpathwaygenes, condmd[,1:10], 
                                             heatmaptitle = hmtitle,
                                             change_gene_label=T,
                                             gene_label_equivalency=gencode,
                                             column_km = 2))
        
        dev.off()
        
        
      } #finish heatmap loop
      
    } #finish subcat loop
    
  } #finish category with sub-category loop
  
  
}



#whole msigdb gsea
cat = 'Whole_MSIGDB'
catdir <-  paste0(compoutdir, '/3.gsea/', cat)
dir.create(catdir)

cp  <- pathways %>% split(x = .$ensembl_gene, f = .$gs_name)


#run GSEA and make main dotplot
gseares <- FerrenaBulkRNAseq::gsea.results(save, pathways = cp)
gseadotplot <- FerrenaBulkRNAseq::gsea.dotplot.onecol(gseares = gseares[[1]], pathwayfontsize = 7) +
  ggtitle(cat)

#save gsearesults
write.csv(file = paste0(catdir, '/fgsea-results.csv'),
          x = gseares[[1]], row.names = F, quote = F)

#save dotplot
ggsave(gseadotplot, 
       filename = paste0(catdir, '/gseadotpot.jpg'), 
       height = 10, width = 10, bg = 'white')




#save leading edge
le <- gseares[[2]]

#make a padded DF from a list of char vectors of different lengths:
# https://stackoverflow.com/questions/43159384/convert-a-list-of-character-vectors-into-dataframe
le <- data.frame(lapply(le, `length<-`, max(lengths(le))))


write.csv(file = paste0(catdir, '/leadingedge.csv'),
          x = le, row.names = F, quote = F)



#loop through significant pathways, make leading edge plots
gsearesdf <- gseares[[1]]
gsearesdf <- gsearesdf[!is.na(gsearesdf$pval),]
sigpathways <- gsearesdf[gsearesdf$padj < 0.05,]


dir.create(paste0(catdir, '/heatmap-leadingedge'))
dir.create(paste0(catdir, '/heatmap-wholepathway'))


#for each significant pathway, 
# make a heatmap using both leading edge and whole pathway genes
for(spname in sigpathways$pathway){
  
  #get the pathway
  sp <- sigpathways[sigpathways$pathway == spname,]
  
  #get the leading edge genes
  sple <- le[[spname]]
  
  #get rid of NAs
  sple <- sple[!is.na(sple)]
  
  #if not enough genes causes problems...
  if(length(sple) < 3){
    next
  }
  
  
  #take only top 75 genes, do we really need to show more? lol
  sple <- sple[1:75]
  
  #get rid of NAs
  sple <- sple[!is.na(sple)]
  
  #use MGI names
  sple_mgi <- gencode[match(sple, gencode$gene_id),]
  
  # #if any are duplicated, use ensemble, else use MGI
  # if(!any(duplicated(sple_mgi$mgi_symbol))){
  #   sple <- sple_mgi$mgi_symbol
  # }
  
  #make a heatmap title
  hmtitle <- paste0(spname, '\n',
                    'NES=', round(sp$NES, 2), ', FDR=', round(sp$padj, 2))
  
  
  
  jpeg(paste0(catdir, '/heatmap-leadingedge/', spname, '.jpg'),
       height = 8, width = 10, units = 'in', res = 300)
  
  
  print(FerrenaBulkRNAseq::heatmapplot(gemnorm, sple, condmd[,1:10], 
                                       heatmaptitle = hmtitle,
                                       change_gene_label=T,
                                       gene_label_equivalency=gencode,
                                       column_km = 2))
  
  dev.off()
  
  #all genes in pathway
  
  allpathwaygenes <- cp[[spname]]
  allpathwaygenes <- allpathwaygenes[allpathwaygenes %in% rownames(gemnorm)]
  
  
  jpeg(paste0(catdir, '/heatmap-wholepathway/', spname, '.jpg'),
       height = 8, width = 10, units = 'in', res = 300)
  
  print(FerrenaBulkRNAseq::heatmapplot(gemnorm, allpathwaygenes, condmd[,1:10], 
                                       heatmaptitle = hmtitle,
                                       change_gene_label=T,
                                       gene_label_equivalency=gencode,
                                       column_km = 2))
  
  dev.off()
  
  
} #finish heatmap loop










writeLines(text = capture.output(sessionInfo()), 
           con = paste0(compoutdir, '/sessionInfo_', gsub(pattern = ' ', replacement = '_', Sys.time() ), '.txt' ) )

beepr::beep()


