library(tidyverse)
library(clusterProfiler)
library(msigdbr)
library(enrichplot)
library(org.Mm.eg.db)

set.seed(2021)

# > packageVersion('clusterProfiler')
# [1] ‘4.0.5’
# > packageVersion('enrichplot')
# [1] ‘1.12.3’

#make out dir
outdir <- 'results/comparative-de/TKO-vs-DKO/4.downstream/clusterprofiler/up'

dir.create(outdir)

### THINGS TO ADJUST ###

# INPUT GENE LIST

# COLOR PALETTE


# input gene list
#read in DE res

res <- read.csv('results/comparative-de/TKO-vs-DKO/2.deresults/DEresults-normalized_count_matrix.csv')

#fix ensembl...
res$ensembl_gene_id <- str_split_fixed(res$ensembl_gene_id, '\\.', 2)[,1]


res <- res[res$padj < 0.05,]
res <- res[res$log2FoldChange > 1,]
res <- res[order(abs(res$log2FoldChange), decreasing = T),]

genenames <- res$mgi_symbol



namedfc <- res$log2FoldChange
names(namedfc) <- res$mgi_symbol


# COLOR PALETTE
# up = Accent, down = 'Set3'
palname = 'Accent'














#get pathways:
# hallmarks
# C2 "CP" (biocarta,kegg,PID,reactome,wikipways)
# C3 TFT
# C5 GO BP/MF/CC
# C6 Oncogenic
# C8 cell type


# pathways <- msigdbr::msigdbr(species = 'Mus musculus')
pathways <- readRDS('~/Dropbox/data/general_data_utilities/msigdb/msigdb_mouse.2022.may.12.rds')
# 
# pathways <- pathways[pathways$gs_cat %in% c('H','C2', 'C3', 'C5', 'C6', 'C8'),]
# 
# h <- pathways[pathways$gs_cat=='H' , ]
# c2 <- pathways[pathways$gs_cat=='C2' & grepl('CP', pathways$gs_subcat), ]
# c3 <- pathways[pathways$gs_cat=='C3' & grepl('TFT', pathways$gs_subcat), ]
# c5 <- pathways[pathways$gs_cat=='C5' & pathways$gs_subcat %in% c('GO:BP', 'GO:MF'), ]
# c6 <- pathways[pathways$gs_cat=='C6' , ]
# c8 <- pathways[pathways$gs_cat=='C8' , ]
# 
# pathways <- rbind(h,c2,c3,c5,c6,c8)
# rm(h,c2,c3,c5,c6,c8)

term2gene <- pathways[,c('gs_name', 'gene_symbol')]

rm(pathways)


## run the enrichment analysis ##

up_msigdb_ora <- enricher(genenames,
                          TERM2GENE = term2gene, 
                          # TERM2NAME = term2name,
                          # OrgDb = org.Mm.eg.db,
                          # keyType = 'ENSEMBL',
                          pvalueCutoff = 0.05
)


#run termsimilarity
up_msigdb_ora <- enrichplot::pairwise_termsim(up_msigdb_ora)

#try to plot color = num genes in overexp gene list?
up_msigdb_ora@result$Percent_of_DEGs <- (up_msigdb_ora@result$Count / length(genenames))*100


set.seed(2021)
options(ggrepel.max.overlaps = 2000)

up_emap_total <- emapplot(up_msigdb_ora, color='Percent_of_DEGs', repel=T, showCategory = 200,  
                          cex_label_category=0.4, cex_category = 0.3,cex_line = 0.3)+
  scale_fill_distiller(palette = 'Reds', direction = 1, name='Percent of DEGs')


up_emap_total










# try to recover clusters....
# we need the object and the plot
obj <- up_msigdb_ora
plot <- up_emap_total


#get termsim and min_dist
mat <- obj@termsim
min_dist = 0.2 #default from package

sig <- res # set significant DE genes, with LFC and pvalues...



#keep only termsim matrix elements that are in final plot
# the top 30... by pvalue?
dat <- plot$data
# mat <- mat[rownames(mat) %in% dat$name ,  colnames(mat) %in% dat$name ]


#try to recapture clusters...

conectlist <- list()
for(i in c(rownames(mat)) ){
  
  #i don't know why the cols and rows are not the same
  # just check both
  x = na.omit( mat[i,] )
  xnames <- names(x[x>=min_dist])
  
  y = na.omit( mat[,i] )
  ynames <- names(y[y>=min_dist])
  
  conectlist[[i]] <- unique(c(i,xnames, ynames))
  
}



#select single=node clusters, remove from list
# singlenodes <- names( conectlist[which(sapply(conectlist, length)==0)] )
# 
# conectlist <- conectlist[!(names(conectlist) %in% singlenodes)]


#biggest first
# sort by number of pathways in each cluster
conectlist <- conectlist[order(sapply(conectlist, length), decreasing = T)]


#pairwise overlap?
# start with first cluster
reslist <- conectlist[1]

for(pway in names(conectlist)[-1] ){
  
  
  message(pway)
  cpws <- conectlist[[pway]]
  
  #thru cpws and reslist
  # if it's in, add it
  # if not, add new to reslist
  for(residx in c(1:length(reslist)) ){
    
    
    respways <- reslist[[residx]]
    
    # check if any of the cpws is in res for this residx and break from this subloop
    if( any(cpws %in% respways) ){
      
      message(' - matches with res ', residx)
      reslist[[residx]] <- unique(c(respways, cpws))
      
      break()
      
    } 
    
    
    #if not, check if there is another reslist index...
    # if not, add
    if(residx == length(reslist)){
      
      message(' - - - no match with res ', residx, ', adding another res index')
      reslist[[pway]] <- cpws
      
    }
    
    
  }
  
  
}




#it seems this can sometimes fail... not sure why...


#try to collapse lists wiht overlap
# maybe we could have done it this way above too...?

backup <- reslist
reslist2 <- list()
skiplist <- c()

for(residx in c(1:length(reslist)) ){
  
  # first, we check if we have previously already added the list...
  if( names(reslist)[residx] %in% skiplist ){next()}
  
  cpws <- reslist[[residx]]
  
  #check the lenght of intersects between this list nd all others
  # if there are, check the intersects of that list and all others and collapse recursively...
  # if not, just add to list
  
  intersectlens <- sapply(reslist[-residx], function(x){length( intersect(x, cpws) )})
  
  if( any(intersectlens > 0) ){
    
    
    message('Intersect found for reslist index : ', residx, '; intersecting indices:')
    
    which_have_intersects <- names(intersectlens[intersectlens>0])
    
    for(j in which_have_intersects){
      
      intersectingidx <- which(names(reslist)==j)
      message(' - reslist index: ',  intersectingidx)
      
      #add to skiplist
      skiplist <- c(skiplist, j)
      
      #update the reslist to include the intersecting pathway...
      cpws <- unique( c(cpws, reslist[[j]]) )
    }
    
    
  }
  
  # add to reslist2
  # only the intersects which were added get ignored...
  reslist2[[residx]] <- cpws
  
  
}


reslist <- reslist2
rm(reslist2)

#finally, we have the result!
# totally agnostic to cluster number


#sort by number of pathways in each cluster
reslist <- reslist[order(sapply(reslist, length), decreasing = T)]

#rename them
names(reslist) <- paste0('Cluster_', seq(1:length(reslist)))





#finally, do some actual analysis...
# get total genes in cluster
# get percent of DEGs from cluster

clustpercs <- list()
cum_perc <- 0
clustgeneslist <- list()

for(clustdex in 1:length(reslist) ){
  clust <- reslist[[clustdex]]
  
  #get the genes for each of these...
  clustergenes <- c()
  for(pw in clust){
    clustergenes <- c(clustergenes, 
                      intersect(pull(term2gene[term2gene$gs_name ==pw,2] ) ,
                                genenames)
    )
    
  }
  
  # get unique genes for this cluster
  # add them to list...
  clustergenes <- unique(clustergenes)
  clustgeneslist[[names(reslist)[clustdex]]] <- clustergenes
  
  #calculaate proportion:
  # length of unique genes in cluster / total significant DEGs
  perc <- (length(clustergenes) / nrow(sig) * 100)
  
  #add to running perc
  cum_perc <- cum_perc + perc
  
  clustpercs[[clustdex]] <- data.frame(clust = names(reslist)[clustdex],
                                       perc = perc,
                                       cum_perc = cum_perc)
  
  
  
}

clustpercs <- bind_rows(clustpercs)



#get unique genes
# also, get cumulative perc using this genes list...
unique_clustgeneslist <- list()
cum_perc_unique_list <- list()
cum_perc = 0

for(clustdex in 1:length(clustgeneslist)){
  
  #get unique genes
  # first get genes from this cluster
  this_clust_genes <- clustgeneslist[[clustdex]]
  
  #get all otger genes
  all_other_genes <- unique(unlist(clustgeneslist[-clustdex]))
  
  #select unique ones
  this_clust_genes <- this_clust_genes[!(this_clust_genes) %in% all_other_genes]
  
  #naame it
  name <- names(clustgeneslist)[clustdex]
  unique_clustgeneslist[[name]] <- this_clust_genes
  
  #calculaate proportion:
  # length of unique genes in cluster / total significant DEGs
  perc <- (length(this_clust_genes) / nrow(sig) * 100)
  
  #add to running perc
  cum_perc <- cum_perc + perc
  cum_perc_unique_list[[name]] <- data.frame(perc = perc, cum_perc=cum_perc)
  
  
}


bind_rows(cum_perc_unique_list)


#### this would be the place to rename clusters, if we want...
clustpercs
clustpercs[clustpercs$perc>=10,]



reslist[names(reslist) %in% clustpercs[clustpercs$perc>=3,]$clust]

clustpercs$relabel <- clustpercs$clust
clustpercs[1:3,'relabel'] <- c('Cluster_1_Immune_system', 'Cluster_2_Actin_and_muscle', 'Cluster_3_Endothelial')
clustpercs$relabel <- paste0(clustpercs$relabel, '\n', round(clustpercs$perc, 1), ' % of DEGs')



#just check any above 10%?
# get percents for the top clusters
sigclusts <- clustpercs[1:3,]
#sigclusts <- clustpercs[which(sapply(reslist, length)>1),]

# get actuaal genes of the top clusters
sigclustgenes <- clustgeneslist[names(clustgeneslist) %in% sigclusts$clust]
names(sigclustgenes) <- sigclusts$relabel

# get pathways in the top clusters
sigclustpways <- reslist[names(reslist) %in% sigclusts$clust]
names(sigclustpways) <- sigclusts$relabel

# get pathwyas in top clusters as a df
sigclustpwaysdf <- bind_rows( lapply(1:length(sigclustpways), function(x){
  i=sigclustpways[[x]]
  name=names(sigclustpways)[x]
  data.frame(pways=i, clust=name)
})
)

# venn of significant clusters
sig_ven <- ggVennDiagram::ggVennDiagram(sigclustgenes)+
  labs(title = 'Venn Diagram of genes from top clusters of enriched pathways')+
  theme(plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(expand = expansion(mult = .2))





### to make a table, get the object results, add in clustering info
#names(reslist) <- clustpercs$relabel
pwayclust <- unlist2(reslist)
pwayclust <- data.frame(Cluster=names(pwayclust),
                        ID = pwayclust)

#add in names and % 
pwayclust$PercDEGs <- clustpercs[match(pwayclust$Cluster, clustpercs$clust ),'perc']

#add clustering info to object res
pwayres <- obj@result
pwayres <- pwayres[pwayres$ID %in% pwayclust$ID,]
pwayres$Cluster <- pwayclust[match(pwayres$ID, pwayclust$ID), 'Cluster']
pwayres$Percent_of_DEGs_CLUSTER <- pwayclust[match(pwayclust$ID, pwayres$ID), "PercDEGs"]


#colors are defined above
levs <- table(pwayres$Cluster)
pal <- RColorBrewer::brewer.pal(Inf, palname)
pal <- pal[1:length(levs)]
pal <- c(pal, 'grey50')


### using coordinates, add labels to origninal emap
dat <- up_emap_total$data
dat$Cluster <- NA
dat[match(sigclustpwaysdf$pways, dat$name),'Cluster'] <- sigclustpwaysdf$clust


repelaggr <- aggregate(x ~ Cluster, dat, mean)
repelaggr$y <- aggregate(y ~ Cluster, dat, mean)[,2]
repelaggr$Color <- pal[1:nrow(repelaggr)]


up_emap_total_withclust <- up_emap_total +  ggrepel::geom_label_repel(inherit.aes = F,
                                                                      data = repelaggr, aes(x=x,y=y,label=Cluster, fill=Cluster),
                                                                      fill= repelaggr$Color, 
                                                                      box.padding = 5, max.overlaps = 200,
                                                                      # direction = 'x',
                                                                      min.segment.length = Inf)


up_emap_total_withclust



#add color to result emaapplot
obj@result$ClusterMain <- NA
obj@result[match(sigclustpwaysdf$pways, obj@result$ID),'ClusterMain'] <- sigclustpwaysdf$clust


#set.seed(2021)



emap_final <-  emapplot(obj, color='ClusterMain', repel=F, showCategory = 200, # layout='graphopt',
                        node_label='None', #cex_label_category=0.4, 
                        coords = dat[,1:2],
                        cex_category = 0.3,cex_line = 0.3 )+
  scale_fill_manual(values = pal)





### MAAKE SUB-PLOTS FOR THE TOP CLUSTERS...

namedfc <- res$log2FoldChange
names(namedfc) <- res$mgi_symbol

toppways <- list()
subclustlist <- list()
for(sigclustidx in c(1:nrow(sigclusts)) ){
  
  origlabel <- sigclusts[sigclustidx, 'clust']
  relabel <- sigclusts[sigclustidx, 'relabel']
  
  message(relabel)
  
  #get the pathways
  sigpways <- sigclustpways[[sigclustidx]]
  
  #recompute...
  
  subemap <- emapplot(obj, color='Percent_of_DEGs', repel=T,
                      #coords = subdat[,1:2],
                      showCategory = sigpways,
                      cex_label_category=0.4, 
                      cex_category = 0.3,cex_line = 0.3
  )+
    scale_fill_distiller(palette = 'Reds', direction = 1, name='Percent of DEGs')+
    ggtitle(relabel) +
    theme(plot.title = element_text(hjust = 0.5))
  
  
  
  #selct top pathways from cluster to show...
  pwayres_inthisclust <- pwayres[pwayres$ID %in% sigpways,]
  
  #select top 3 based n num DEGs
  num_to_pick <- ifelse(nrow(pwayres_inthisclust) >= 3, 3, nrow(pwayres_inthisclust))
  pwayres_inthisclust <- pwayres_inthisclust[order(pwayres_inthisclust$Percent_of_DEGs, decreasing = T),]
  pways_to_show <- pwayres_inthisclust$ID[1:num_to_pick]
  
  subcnet <- cnetplot(obj, showCategory = pways_to_show, foldChange = namedfc,
                      shadowtext='category',
                      cex_gene = 0.5, cex_label_gene = 0.5, cex_label_category = 0.7)+
    ggtitle(relabel) +
    theme(plot.title = element_text(hjust = 0.5))
  
  
  
  
  # subora@result$qvalue <- subora@result$Percent_of_DEGs
  subddotplot <-  dotplot(obj, showCategory=head(pwayres_inthisclust$ID,30),font.size=6)+
    ggtitle(relabel) +
    theme(plot.title = element_text(hjust = 0.5))
  
  
  
  
  # pdf(paste0(clustoutdir, '/subplots_', origlabel, '.pdf'), width = 8, height = 8)
  # print(subemap)
  # print(subcnet)
  # print(subddotplot)
  # 
  # dev.off()
  
  #for cross cnet: get top 3 pways?
  toppways[[sigclustidx]] <- pways_to_show
  
  subclustlist[[sigclustidx]] <- list(subemap,subcnet,subddotplot)
  
}


#get top pways from sig clusters
top_1_pway <- sapply(toppways, function(x){head(x,1)})
crosscnet <- cnetplot(obj, showCategory = top_1_pway, foldChange = namedfc,
                      shadowtext='category',
                      cex_gene = 0.5, cex_label_gene = 0.5, cex_label_category = 0.7)+
  ggtitle('top pathways from sig clusters') +
  theme(plot.title = element_text(hjust = 0.5))

crossdotplot <-  dotplot(obj, showCategory=unlist(toppways))+
  ggtitle('top pathways from sig clusters') +
  theme(plot.title = element_text(hjust = 0.5))+
  #scale_color_distiller(palette = 'Purples', direction = -1)+
  scale_color_gradient2(low="#54278F",high='gray', midpoint = 0.05)



#get gtop pways overall, default...
defcnet <- cnetplot(obj, showCategory = 3, foldChange = namedfc,
                    shadowtext='category',
                    cex_gene = 0.5, cex_label_gene = 0.5, cex_label_category = 0.7)+
  ggtitle('top pathways overall') +
  theme(plot.title = element_text(hjust = 0.5))

defdotplot <- dotplot(obj, showCategory=10)+
  ggtitle('top pathways overall') +
  theme(plot.title = element_text(hjust = 0.5))+
  #scale_color_distiller(palette = 'Purples', direction = -1)+
  scale_color_gradient2(low="#54278F",high='gray', midpoint = 0.05)




###save main plots

pdf(paste0(outdir, '/plots.pdf'), width = 8, height = 8)
print(up_emap_total)
print(up_emap_total_withclust)



print(emap_final)

crosscnet
crossdotplot

print(defcnet)
print(defdotplot)

print(sig_ven)

subclustlist

dev.off()




write.csv(file=paste0(outdir, '/pathway_resuts.csv'), pwayres, row.names = F)





#save dat and clusterprofilerobj

saveRDS(obj, paste0(outdir, '/clusterprofilerobj.rds'))
saveRDS(dat, paste0(outdir, '/emapcoords.rds'))



#update, write out the full pway res...
obj <- readRDS(paste0(outdir, '/clusterprofilerobj.rds'))
pwayres <- read.csv(paste0(outdir, '/pathway_resuts.csv'))

fullres <- obj@result
fullres <- fullres[fullres$p.adjust < 0.05,]
fullres$Cluster <- c( pwayres$Cluster, rep(NA, nrow(fullres)-nrow(pwayres)) )
fullres$Percent_of_DEGs_CLUSTER <- c( pwayres$Percent_of_DEGs_CLUSTER, rep(NA, nrow(fullres)-nrow(pwayres)) )
write.csv(file=paste0(outdir, '/pathway_results_fullsignificant.csv'), fullres, row.names = F)



beepr::beep()


### SPECIFIC PLOTS ###

pdf(paste0(outdir, '/finalplots_emap.pdf'), width = 5, height = 5)
print(emap_final)
dev.off()



pdf(paste0(outdir, '/finalplots_dotplot.pdf'), width = 5, height = 5)

# dotplot(obj, showCategory=10, label_format=25, font.size = 10)+
#   ggtitle('top pathways overall') +
#   theme(plot.title = element_text(hjust = 0.5))+
#   #scale_color_distiller(palette = 'Purples', direction = -1)+
#   scale_color_gradient2(low="#54278F",high='gray', midpoint = 0.05)

dotplot(obj, showCategory=10, label_format=20, font.size = 9)+
  #ggtitle('top pathways from sig clusters') +
  #theme(plot.title = element_text(hjust = 0.5))+
  #scale_color_distiller(palette = 'Purples', direction = -1)+
  scale_color_gradient2(low="firebrick",high='gray', midpoint = 0.05)

dev.off()





pdf(paste0(outdir, '/finalplots_cnet.pdf'), width = 8, height = 7)

pways <- c('DESCARTES_FETAL_CEREBRUM_MICROGLIA', 
           "DESCARTES_FETAL_CEREBRUM_VASCULAR_ENDOTHELIAL_CELLS",
           'GOMF_ACTIN_BINDING'
)
cnetplot(obj, showCategory = pways, foldChange = namedfc,
         shadowtext='category',
         cex_gene = 0.3, cex_label_gene = 0.4, cex_label_category = 0.8)+
  ggtitle('top pathways overall') +
  theme(plot.title = element_text(hjust = 0.5))


dev.off()


options(ggrepel.max.overlaps = Inf)
pdf(paste0(outdir, '/finalplots_cnet_UPDATED.pdf'), width = 7, height = 7)

pways <- c('GOBP_LEUKOCYTE_MIGRATION', 
           "GOBP_ENDOTHELIUM_DEVELOPMENT",
           'GOBP_MUSCLE_SYSTEM_PROCESS'
)


cnetplot(obj, showCategory = pways, foldChange = namedfc,
         shadowtext='category',
         node_label = 'all',
         cex_gene = 0.3, cex_label_gene = 0.5, cex_label_category = 0, 
         max.overlaps = Inf, layout = 'kk')



dev.off()




#### heatmaps stratify by pathway

pathways <- readRDS('~/Dropbox/data/general_data_utilities/msigdb/msigdb_mouse.2022.may.12.rds')

obj <- readRDS(paste0(outdir, '/clusterprofilerobj.rds'))


### get only significant pathways
# get significant enriched pways
pwayres <- obj@result
pwayres <- pwayres[pwayres$p.adjust < 0.05,]

#keep only pathway, category, subcat
pathways <- pathways[,c('gs_cat', 'gs_subcat', 'gs_name')]
pathways <- pathways[!duplicated(pathways$gs_name),]

#keep only pathways in signficant enriched res
pathways <- pathways[match(pwayres$ID, pathways$gs_name),]



## select subcats of interest
h <- pathways[pathways$gs_cat=='H' , ]
c2 <- pathways[pathways$gs_cat=='C2' & grepl('CP', pathways$gs_subcat), ]
c3 <- pathways[pathways$gs_cat=='C3' & grepl('TFT', pathways$gs_subcat), ]
c5 <- pathways[pathways$gs_cat=='C5' & pathways$gs_subcat %in% c('GO:BP', 'GO:MF', 'GO:CC'), ]
c6 <- pathways[pathways$gs_cat=='C6' , ]
c8 <- pathways[pathways$gs_cat=='C8' , ]

pathways_restricted <- rbind(h,c2,c3,c5,c6,c8)
pathways_restricted <- pathways_restricted[match(pwayres[pwayres$ID %in% pathways_restricted$gs_name,'ID'], pathways_restricted$gs_name),]

#make PDF
pdf(paste0(outdir, '/dotplots_bycategory.pdf'), width = 5, height = 5)


dotplot(obj, showCategory=20, label_format=25, font.size = 6)+
  ggtitle('top pathways overall') +
  theme(plot.title = element_text(hjust = 0.5))+
  #scale_color_distiller(palette = 'Purples', direction = -1)+
  scale_color_gradient2(low="firebrick",high='gray', midpoint = 0.05)



dotplot(obj, showCategory=pathways_restricted$gs_name[1:20], label_format=25, font.size = 6)+
  ggtitle('top pathways, restricted categories') + 
theme(plot.title = element_text(hjust = 0.5))+
  #scale_color_distiller(palette = 'Purples', direction = -1)+
  scale_color_gradient2(low="firebrick",high='gray', midpoint = 0.05)

### hallmarks
dotplot(obj, showCategory=h$gs_name[1:10], 
        label_format=20, font.size = 7)+
  ggtitle('Hallmarks') +
  #theme(plot.title = element_text(hjust = 0.5))+
  #scale_color_distiller(palette = 'Purples', direction = -1)+
  scale_color_gradient2(low="firebrick",high='gray', midpoint = 0.05)

### c2

dotplot(obj, showCategory=c2$gs_name[1:10], 
        label_format=20, font.size = 7)+
  ggtitle('C2 CP (curated gene sets, Canonical Pathways)') +
  #theme(plot.title = element_text(hjust = 0.5))+
  #scale_color_distiller(palette = 'Purples', direction = -1)+
  scale_color_gradient2(low="firebrick",high='gray', midpoint = 0.05)

c2subs <- unique(c2$gs_subcat)
for(sub in c2subs){
  subpwaydf <- c2[c2$gs_subcat %in% sub,]
  
  
  print( dotplot(obj, showCategory=subpwaydf$gs_name[1:10], 
          label_format=20, font.size = 7)+
    ggtitle(paste0('C2: ', sub)) +
    #theme(plot.title = element_text(hjust = 0.5))+
    #scale_color_distiller(palette = 'Purples', direction = -1)+
    scale_color_gradient2(low="firebrick",high='gray', midpoint = 0.05)
  )
  
}



### c3

dotplot(obj, showCategory=c3$gs_name[1:10], 
        label_format=20, font.size = 7)+
  ggtitle('C3 Regulatory target, TFs') +
  #theme(plot.title = element_text(hjust = 0.5))+
  #scale_color_distiller(palette = 'Purples', direction = -1)+
  scale_color_gradient2(low="firebrick",high='gray', midpoint = 0.05)




### c5 

dotplot(obj, showCategory=c5$gs_name[1:10], 
        label_format=20, font.size = 7)+
  ggtitle('C5, ontology') +
  #theme(plot.title = element_text(hjust = 0.5))+
  #scale_color_distiller(palette = 'Purples', direction = -1)+
  scale_color_gradient2(low="firebrick",high='gray', midpoint = 0.05)

c5subs <- unique(c5$gs_subcat)
for(sub in c5subs){
  subpwaydf <- c5[c5$gs_subcat %in% sub,]
  
  
  print( dotplot(obj, showCategory=subpwaydf$gs_name[1:10], 
                 label_format=20, font.size = 7)+
           ggtitle(paste0('c5: ', sub)) +
           #theme(plot.title = element_text(hjust = 0.5))+
           #scale_color_distiller(palette = 'Purples', direction = -1)+
           scale_color_gradient2(low="firebrick",high='gray', midpoint = 0.05)
  )
  
}



### c6


dotplot(obj, showCategory=c6$gs_name[1:10], 
        label_format=20, font.size = 7)+
  ggtitle('C6, oncogenic signatures') +
  #theme(plot.title = element_text(hjust = 0.5))+
  #scale_color_distiller(palette = 'Purples', direction = -1)+
  scale_color_gradient2(low="firebrick",high='gray', midpoint = 0.05)




### c8


dotplot(obj, showCategory=c8$gs_name[1:10], 
        label_format=20, font.size = 7)+
  ggtitle('C8 cell type signatures') +
  #theme(plot.title = element_text(hjust = 0.5))+
  #scale_color_distiller(palette = 'Purples', direction = -1)+
  scale_color_gradient2(low="firebrick",high='gray', midpoint = 0.05)



dev.off()






#### final heatmap: pick 2-3 from Hallmark, GOBP, C6


### hallmarks
dotplot(obj, showCategory=h$gs_name[1:10], 
        label_format=20, font.size = 7)+
  ggtitle('Hallmarks') +
  #theme(plot.title = element_text(hjust = 0.5))+
  #scale_color_distiller(palette = 'Purples', direction = -1)+
  scale_color_gradient2(low="firebrick",high='gray', midpoint = 0.05)


### gobp
sub = "GO:BP"
subpwaydf <- c5[c5$gs_subcat %in% sub,]


dotplot(obj, showCategory=subpwaydf$gs_name[1:20], 
               label_format=20, font.size = 7)+
         ggtitle(paste0('c5: ', sub)) +
         #theme(plot.title = element_text(hjust = 0.5))+
         #scale_color_distiller(palette = 'Purples', direction = -1)+
         scale_color_gradient2(low="firebrick",high='gray', midpoint = 0.05)


dotplot(obj, showCategory=c6$gs_name[1:30], 
        label_format=20, font.size = 7)+
  ggtitle('C6, oncogenic signatures') +
  #theme(plot.title = element_text(hjust = 0.5))+
  #scale_color_distiller(palette = 'Purples', direction = -1)+
  scale_color_gradient2(low="firebrick",high='gray', midpoint = 0.05)







#pick categories
mycats <- data.frame(category = rep(c( 'MSIGDB Hallmarks', 'GO: Biological Process', 'C6:Oncogenic'), each=3),
                     ID = c('HALLMARK_MYOGENESIS', 'HALLMARK_ALLOGRAFT_REJECTION', 'HALLMARK_KRAS_SIGNALING_UP',
                            'GOBP_MUSCLE_SYSTEM_PROCESS', 'GOBP_T_CELL_ACTIVATION', 'GOBP_LEUKOCYTE_MIGRATION',
                            "SNF5_DN.V1_UP", "MEK_UP.V1_UP" , "KRAS.LUNG_UP.V1_UP" )
                     )


### try adding in ifn
addcat <- data.frame(category = "MSIGDB Hallmarks",
                     ID = "HALLMARK_INTERFERON_GAMMA_RESPONSE")
mycats <- rbind(mycats, addcat)



#get pvalue, padj, and %DEGs
mycats <- cbind(mycats, pwayres[match(mycats$ID, pwayres$ID),c('Percent_of_DEGs', 'p.adjust', 'Count', 'GeneRatio')] )
mycats$GeneRatio <- DOSE::parse_ratio(mycats$GeneRatio)

#hallmark, then goBP, then c6
mycats$category <- factor(mycats$category, levels = c( 'MSIGDB Hallmarks', 'GO: Biological Process', 'C6:Oncogenic') )

#make sure order is right
mycats <- mycats[order(mycats$category),]

#fix long pway names
pnames <- mycats$ID 
pnewnames <- c()
for(p in pnames){
  
  
  if(str_length(p) > 20){
    
    #try to find and replace underscore closest to halfway...
    halfway <- round(str_length(p)/2)
    underscorepos <- str_locate_all(p,'_')[[1]][,1]
    
    distances <- abs(halfway-underscorepos)
    
    which_und_is_closest <- which.min(distances)
    
    split <- str_split_fixed(p,'_',Inf)
    
    newp <- paste0(paste(split[1,1:which_und_is_closest], collapse = '_'),
                   '\n',
                   paste(split[1,(which_und_is_closest+1):ncol(split)], collapse = '_')
    )
    
    
    
  } else{newp <- p}
  
  pnewnames <- c(pnewnames, newp)
  
}

mycats$ID <- pnewnames

#order pathways by pvalue
mycats$ID <- factor(mycats$ID, levels = mycats[order(mycats$GeneRatio, decreasing = F),'ID'])



  


pdf(paste0(outdir, '/finalplots_dotplot_STRAT-BY-CATEGORY_UPDATED_withIFN.pdf'), width = 5, height = 5)

ggplot(mycats, aes(GeneRatio, ID, size = Percent_of_DEGs, col = p.adjust))+
  geom_point() + 
  facet_wrap(~category, nrow = 3, scales='free_y')+
  scale_color_gradient2(low="firebrick",high='gray', midpoint = 0.05)+
  theme_linedraw()+
  theme(axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 10)) +
  ylab('')
  

dev.off()









#pick categories
mycats <- data.frame(category = c(rep( 'MSIGDB Hallmarks', 3),
                                  rep('C3: TF Targets',4),
                                  rep('C6: Oncogenic', 3)),
                     ID = c('HALLMARK_MYOGENESIS', 'HALLMARK_ALLOGRAFT_REJECTION', 'HALLMARK_KRAS_SIGNALING_UP',
                            'CTAWWWATA_RSRFC4_Q2', 'E12_Q6', 'COREBINDINGFACTOR_Q6', 'AP4_Q6_01',
                            "SNF5_DN.V1_UP", "MEK_UP.V1_UP" , "KRAS.LUNG_UP.V1_UP" )
)


#pick categories
mycats <- data.frame(category = c(rep( 'MSIGDB Hallmarks', 4),
                                  rep('C3: TF Targets',4),
                                  rep('C6: Oncogenic', 3)),
                     ID = c('HALLMARK_MYOGENESIS', 'HALLMARK_ALLOGRAFT_REJECTION', 'HALLMARK_KRAS_SIGNALING_UP',
                            'CTAWWWATA_RSRFC4_Q2', 'E12_Q6', 'COREBINDINGFACTOR_Q6', 'AP4_Q6_01',
                            "SNF5_DN.V1_UP", "MEK_UP.V1_UP" , "KRAS.LUNG_UP.V1_UP" )
)



#get pvalue, padj, and %DEGs
mycats <- cbind(mycats, pwayres[match(mycats$ID, pwayres$ID),c('Percent_of_DEGs', 'p.adjust', 'Count', 'GeneRatio')] )
mycats$GeneRatio <- DOSE::parse_ratio(mycats$GeneRatio)

#hallmark, then goBP, then c6
mycats$category <- factor(mycats$category, levels = c( 'MSIGDB Hallmarks', 'C3: TF Targets', 'C6: Oncogenic') )

#save ID for cnetplots
mycats$oldid <- mycats$ID

# #fix long pway names
# pnames <- mycats$ID 
# pnewnames <- c()
# for(p in pnames){
#   
#   
#   if(str_length(p) > 20){
#     
#     #try to find and replace underscore closest to halfway...
#     halfway <- round(str_length(p)/2)
#     underscorepos <- str_locate_all(p,'_')[[1]][,1]
#     
#     distances <- abs(halfway-underscorepos)
#     
#     which_und_is_closest <- which.min(distances)
#     
#     split <- str_split_fixed(p,'_',Inf)
#     
#     newp <- paste0(paste(split[1,1:which_und_is_closest], collapse = '_'),
#                    '\n',
#                    paste(split[1,(which_und_is_closest+1):ncol(split)], collapse = '_')
#     )
#     
#     
#     
#   } else{newp <- p}
#   
#   pnewnames <- c(pnewnames, newp)
#   
# }
# 
# mycats$ID <- pnewnames


### edit the tf target names

mycats[mycats$ID == 'CTAWWWATA_RSRFC4_Q2', 'ID'] <- 'MEF2A targets - "CTAWWWATA_RSRFC4_Q2"'
mycats[mycats$ID == 'E12_Q6', 'ID'] <- 'TCF3 targets - "E12_Q6"'
mycats[mycats$ID == 'COREBINDINGFACTOR_Q6', 'ID'] <- 'CBFA2T2/3 targets - "COREBINDINGFACTOR_Q6"'
mycats[mycats$ID == 'AP4_Q6_01', 'ID'] <- 'TFAP4 targets - "AP4_Q6_01"'

#order pathways by pvalue
mycats$ID <- factor(mycats$ID, levels = mycats[order(mycats$GeneRatio, decreasing = F),'ID'])






pdf(paste0(outdir, '/finalplots_dotplot_STRAT-BY-CATEGORY_UPDATED.pdf'), width = 9, height = 5)

ggplot(mycats, aes(GeneRatio, ID, size = Percent_of_DEGs, col = p.adjust))+
  geom_point() + 
  facet_wrap(~category, nrow = 3, scales='free_y')+
  scale_color_gradient2(low="firebrick",high='gray', midpoint = 0.05)+
  theme_linedraw()+
  theme(axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 15)) +
  ylab('')




dev.off()





cnetplot(obj, showCategory = mycats[mycats$category == 'C3: TF Targets', 'oldid'], 
         foldChange = namedfc,
         shadowtext='category',
         node_label = 'all',
         cex_gene = 0.3, cex_label_gene = 0.5, cex_label_category = 0.8, 
         max.overlaps = Inf, layout = 'kk')




cnetplot(obj, showCategory = mycats[mycats$category == 'C6: Oncogenic', 'oldid'], 
         foldChange = namedfc,
         shadowtext='category',
         node_label = 'all',
         cex_gene = 0.3, cex_label_gene = 0.5, cex_label_category = 0.8, 
         max.overlaps = Inf, layout = 'kk')
