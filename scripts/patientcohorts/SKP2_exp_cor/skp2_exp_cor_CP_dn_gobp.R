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
outdir <- "nonsurvival/skp2cor_tkoenrichment/cor_dn_updated"

dir.create(outdir)

### THINGS TO ADJUST ###

# INPUT GENE LIST

# COLOR PALETTE


# input gene list
#read in DEGs
res <- read.csv('nonsurvival/skp2cor_tkoenrichment/skp2_cor_res.csv')

#fillter
res <- res[res$fdr < 0.05,]
res <- res[res$cor < 0,]

#sort bby effect size
res <- res[order(abs(res$cor), decreasing = T),]

genenames <- res$gene



# COLOR PALETTE
# up = Accent, down = 'Set3'
palname = 'Set3'














#get pathways:
# just GO BP


# pathways <- msigdbr::msigdbr(species = 'Homo sapiens')
# pathways <- msigdbr::msigdbr(species = 'Mus musculus')

pathways <- readRDS('~/Dropbox/data/general_data_utilities/msigdb/msigdb_human.2022.may.12.rds')


pathways <- pathways[pathways$gs_cat=='C5',]
pathways <- pathways[grepl('GO', pathways$gs_subcat),]
pathways <- pathways[pathways$gs_subcat!='GO:CC',]

#just GO BP; can revist others later
pathways <- pathways[pathways$gs_subcat=='GO:BP',]



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


reslist[names(reslist) %in% clustpercs[clustpercs$perc>=10,]$clust]

clustpercs$relabel <- clustpercs$clust
# clustpercs[1,'relabel'] <- c('Cluster_1_Cardiac_Development', 'Cluster_2_MTOR_signalling', 'Cluster_3_Carbohydryate_Metabolism')
clustpercs$relabel <- paste0(clustpercs$relabel, '\n', round(clustpercs$perc, 1), ' % of DEGs')



#just check any above 10%?
# get percents for the top clusters
# sigclusts <- clustpercs[clustpercs$perc>=10,]
sigclusts <- clustpercs[1,]

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
# sig_ven <- ggVennDiagram::ggVennDiagram(sigclustgenes)+
#   labs(title = 'Venn Diagram of genes from top clusters of enriched pathways')+
#   theme(plot.title = element_text(hjust = 0.5))+
#   scale_x_continuous(expand = expansion(mult = .2))





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
obj@result$ClusterMain <- 'Other'
obj@result[match(sigclustpwaysdf$pways, obj@result$ID),'ClusterMain'] <- sigclustpwaysdf$clust


#set.seed(2021)



emap_final <-  emapplot(obj, color='ClusterMain', repel=F, showCategory = 200, # layout='graphopt',
                        node_label='None', #cex_label_category=0.4, 
                        coords = dat[,1:2],
                        cex_category = 0.3,cex_line = 0.3 )+
  scale_fill_manual(values = pal, name = 'Cluster')





### MAAKE SUB-PLOTS FOR THE TOP CLUSTERS...

namedfc <- res$cor
names(namedfc) <- res$gene

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
  
  
  subclustlist[[sigclustidx]] <- list(subemap,subcnet,subddotplot)
  
}





###save main plots

pdf(paste0(outdir, '/plots.pdf'), width = 8, height = 8)
print(up_emap_total)
print(up_emap_total_withclust)

print(emap_final)
#print(sig_ven)

subclustlist

dev.off()


write.csv(file=paste0(outdir, '/pathway_resuts.csv'), pwayres, row.names = F)



#save dat and clusterprofilerobj

saveRDS(obj, paste0(outdir, '/clusterprofilerobj.rds'))
saveRDS(dat, paste0(outdir, '/emapcoords.rds'))



beepr::beep()


### downsteam###


obj <- readRDS(paste0(outdir, '/clusterprofilerobj.rds'))


pdf(paste0(outdir, '/alldotplot.pdf'), width = 5, height = 5)
dotplot(obj, showCategory=20, font.size=5)
dev.off()