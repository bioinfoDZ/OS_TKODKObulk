library(tidyverse)
library(clusterProfiler)
library(msigdbr)
library(enrichplot)
library(org.Mm.eg.db)

set.seed(2021)


#make out dir
outdir <- 'results/comparative-de/TKO-vs-DKO/4.downstream/'

dir.create(outdir)

#read in DE res

res <- read.csv('results/comparative-de/TKO-vs-DKO/2.deresults/DEresults-normalized_count_matrix.csv')

#fix ensembl...
res$ensembl_gene_id <- str_split_fixed(res$ensembl_gene_id, '\\.', 2)[,1]


# try to set ggrepel
# not sure this actually fixes things though...
options(ggrepel.max.overlaps = 20)


#wholemsigdb

pathways <- msigdbr::msigdbr(species = 'Mus musculus')



# just keep all pathways
# it will match the msigdb analysis we already did...

#exclude
# # we exclude any 'legacy' genesets, maybe old and superceded...
# # we also exclude CGP from C2, often pretty confusing...
# # we also exclude the immunologic ones, also often veyr confusing and difficult to interpret
# subcatexclude <- c('MIR:MIR_Legacy', 'TFT:TFT_Legacy', 'CGP')
# catexclude <- c('C7')
# 
# 
# pathways <- pathways[!(pathways$gs_cat %in% catexclude),]
# pathways <- pathways[!(pathways$gs_subcat %in% subcatexclude),]


#prep for gsea...
#term2gene: pathways and genes




term2gene <- pathways[,c('gs_name', 'ensembl_gene')]

rm(pathways)


#term2name: pathways and categories?
# term2name <- pathways[!duplicated(pathways$gs_name), c('gs_name', 'gs_subcat')]
# not sure what this does
# it does not seem to add categories to emapplot






### first, do overexp genes

suboutdir <- paste0(outdir, '/overexp')
dir.create(suboutdir)



#select significant genes...
sig <- res[res$padj < 0.05,]
sig <- sig[sig$log2FoldChange > 1,]


#rank by FC
genes <- sig$log2FoldChange
names(genes) <- sig$ensembl_gene_id

#sort them by decreasing magnitude
genes <- genes[order(abs(genes), decreasing = T)]

# just in case, remove any NAs...
genes <- na.omit(genes)


#just give the actual gene names as input...
genenames <- names(genes)






## run the enrichment analysis ##

up_msigdb_ora <- enricher(genenames,
                          TERM2GENE = term2gene, 
                          # TERM2NAME = term2name,
                          # OrgDb = org.Mm.eg.db,
                          # keyType = 'ENSEMBL',
                          pvalueCutoff = 0.001
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


### for the real analysis:
# network with all pathways in top 200 by pval
# the top 200 is the default in the termsim...

set.seed(2021)
up_emap_total <- emapplot(up_msigdb_ora, color='Percent_of_DEGs', repel=T, showCategory = 200,  layout='graphopt',
                          cex_label_category=0.4, cex_category = 0.3,cex_line = 0.3,seed=2021)+
  scale_fill_distiller(palette = 'Reds', direction = 1, name='Percent of DEGs')


up_emap_total


# try to recover clusters....
# we need the object and the plot
obj <- up_msigdb_ora
plot <- up_emap_total


#get termsim and min_dist
mat <- obj@termsim
min_dist = 0.2 #default from package



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
                      intersect(term2gene[term2gene$gs_name ==pw,"ensembl_gene"]$ensembl_gene, 
                                sig$ensembl_gene_id)
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

clustpercs[clustpercs$perc>=10,]



reslist[names(reslist) %in% clustpercs[clustpercs$perc>=10,]$clust]

#get the genes for unified pathway analysis...
as.data.frame( res[match(clustgeneslist[[3]], res$ensembl_gene_id), "mgi_symbol"] )


clustpercs$relabel <- clustpercs$clust
clustpercs[clustpercs$perc>=10,'relabel'] <- c('Cluster_1_Immune', 'Cluster_2_Muscle', 'Cluster_3_Endothelial')
clustpercs$relabel <- paste0(clustpercs$relabel, '\n', round(clustpercs$perc, 1), ' % of DEGs')



#the cumulative percentage is not relly helpful, 
# gene lists can have some minimal overlap (defined by min_dist)

#just check any above 10%?
# get percents for the top clusters
sigclusts <- clustpercs[clustpercs$perc>=10,]

# get actuaal genes of the top clusters
sigclustgenes <- clustgeneslist[names(clustgeneslist) %in% sigclusts$clust]
names(sigclustgenes) <- sigclusts$relabel

# get pathways in the top clusters
sigclustpways <- reslist[names(reslist) %in% sigclusts$clust]
names(sigclustpways) <- sigclusts$relabel



#ggvenndiagram of the top ones....
up_venn <- ggVennDiagram::ggVennDiagram(sigclustgenes)+
  labs(title = 'Venn Diagram of genes from top clusters of enriched pathways')+
  theme(plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(expand = expansion(mult = .2))


#try to plot...
# get sigclustpwasys as a df
sigclustpwaysdf <- bind_rows( lapply(1:length(sigclustpways), function(x){
  i=sigclustpways[[x]]
  name=names(sigclustpways)[x]
  data.frame(pways=i, clust=name)
})
)


#add color to result emaapplot
up_msigdb_ora@result$Cluster <- NA
up_msigdb_ora@result[match(sigclustpwaysdf$pways, up_msigdb_ora@result$ID),'Cluster'] <- sigclustpwaysdf$clust


set.seed(2021)
up_emap_total_redone <-  emapplot(up_msigdb_ora, color='Cluster', repel=T, showCategory = 200, layout='graphopt',
                                  node_label='None', #cex_label_category=0.4, 
                                  cex_category = 0.3,cex_line = 0.3 )+
  scale_fill_brewer(palette = 'Set2', direction = 1, name='Cluster')

#add ggrepel labels to emapplot
dat <- up_emap_total_redone$data
dat$Cluster <- NA
dat[match(sigclustpwaysdf$pways, dat$name),'Cluster'] <- sigclustpwaysdf$clust


repelaggr <- aggregate(x ~ Cluster, dat, mean)
repelaggr$y <- aggregate(y ~ Cluster, dat, mean)[,2]

#add percentage to the label name...

up_emap_total_redone_withrepel <- up_emap_total_redone +  ggrepel::geom_label_repel(inherit.aes = F, 
                                                                                    data = repelaggr, aes(x=x,y=y,fill=Cluster,label=Cluster), 
                                                                                    min.segment.length = Inf)




#plot just the sig clusters...
subdat <- dat[!is.na(dat$Cluster),]

# to be able to "fill" by cluster in the ggrepel labels,
# we explicitly set fill outside of aes with colors defined before plot call
repelaggr$Color <- RColorBrewer::brewer.pal(nrow(repelaggr), 'Set2')


up_emap_total_justclusters <- emapplot(up_msigdb_ora, color='Percent_of_DEGs', repel=T, showCategory = subdat$name, #layout='graphopt',
                                       cex_label_category=0.4, 
                                       cex_category = 0.3,cex_line = 0.3,
                                       coords = subdat[,1:2])+
  scale_fill_distiller(palette = 'Reds', direction = 1, name='Percent_of_DEGs')




up_emap_total_justclusters_withrepel <- up_emap_total_justclusters +
  ggrepel::geom_label_repel(inherit.aes = F,
                            data = repelaggr, aes(x=x,y=y,label=Cluster),
                            fill= repelaggr$Color, 
                            #box.padding = 10, max.overlaps = 200,
                            # direction = 'x',
                            min.segment.length = Inf)





up_emap_total_justclusters_colorbyclust <- emapplot(up_msigdb_ora, color='Cluster', repel=T, 
                                                    showCategory = subdat$name, #layout='graphopt',
                                                    cex_label_category=0.4, 
                                                    cex_category = 0.3,cex_line = 0.3,
                                                    # coords = subdat[,1:2]
)+
  scale_fill_brewer(palette = 'Set2', direction = 1, name='Cluster')



# up_emap_total_redone
# up_emap_total_redone_withrepel
# 
# up_emap_total / up_emap_total_redone

up_emap_total + up_emap_total_redone_withrepel + 
  up_emap_total_justclusters + up_emap_total_justclusters_withrepel



#plot the basic plots....
ggsave(paste0(suboutdir, '/emap_pval.jpg'), up_emap_pval, width= 10, height= 10)
ggsave(paste0(suboutdir, '/emap_perc.jpg'), up_emap_perc, width= 10, height= 10)

#olot the total plot...
ggsave(paste0(suboutdir, '/emap_total.jpg'), up_emap_total, width= 10, height= 10)


#plot the cluster plots
ggsave(paste0(suboutdir, '/emap_total_cluster.jpg'), up_emap_total_redone, width= 10, height= 10)
ggsave(paste0(suboutdir, '/emap_total_cluster_withrepel.jpg'), up_emap_total_redone_withrepel, width= 10, height= 10)


# plot the total plots with just the top cluster nodes...
ggsave(paste0(suboutdir, '/emap_total_justsig.jpg'), up_emap_total_justclusters, width= 10, height= 10)
ggsave(paste0(suboutdir, '/emap_total_justsig_withrepel.jpg'), up_emap_total_justclusters_withrepel, width= 10, height= 10)

ggsave(paste0(suboutdir, '/emap_total_justsig_colorbyclust.jpg'), up_emap_total_justclusters_colorbyclust, width= 10, height= 10)

#save the venn
ggsave(paste0(suboutdir, '/venn.jpg'), up_venn, width= 10, height= 10, bg='white')


#plot the nice combo
combo <- up_emap_total + up_emap_total_redone_withrepel + 
  up_emap_total_justclusters + up_emap_total_justclusters_withrepel

ggsave(paste0(suboutdir, '/emap_combo.jpg'), combo, width = 14, height = 10)






#try to make one big one....
options(ggrepel.max.overlaps = 4)

#add color to result emaapplot
up_msigdb_ora@result$Cluster <- 'Other'
up_msigdb_ora@result[match(sigclustpwaysdf$pways, up_msigdb_ora@result$ID),'Cluster'] <- sigclustpwaysdf$clust


#set.seed(2021)
levs <- table(up_msigdb_ora@result$Cluster)
pal <- RColorBrewer::brewer.pal(name = 'Set2', n = (length(levs) - 1))
pal <- c(pal, 'grey50')

emap_final <-  emapplot(up_msigdb_ora, color='Cluster', repel=F, showCategory = 200, # layout='graphopt',
                                  node_label='None', #cex_label_category=0.4, 
                                  cex_category = 0.5,cex_line = 0.5 )+
  scale_fill_manual(values = pal, name = 'Cluster')


ggsave(paste0(suboutdir, '/emap_final.jpg'), emap_final, width= 7, height= 7)


options(ggrepel.max.overlaps = 20)





### MAAKE SUB-PLOTS FOR THE TOP CLUSTERS...
clustoutdir <- paste0(suboutdir, '/pathway_clusters/')
dir.create(clustoutdir)


for(sigclustidx in c(1:nrow(sigclusts))){
  
  origlabel <- sigclusts[sigclustidx, 'clust']
  relabel <- sigclusts[sigclustidx, 'relabel']
  
  message(relabel)
  
  #get the pathways
  sigpways <- sigclustpways[[sigclustidx]]
  
  #recompute...
  
  subterm2gene <- term2gene[term2gene$gs_name %in% sigpways,]
  subora <- enricher(genenames,
                     TERM2GENE = subterm2gene
                     # TERM2NAME = term2name,
                     # OrgDb = org.Mm.eg.db,
                     # keyType = 'ENSEMBL',
                     # pvalueCutoff = 0.001,
  )
  
  
  subora <- enrichplot::pairwise_termsim(subora)
  
  
  #try to plot color = num genes in overexp gene list?
  subora@result$Percent_of_DEGs <- (subora@result$Count / length(genenames))*100
  
  #get dat with coords for this cluster
  subdat <- dat[dat$name %in% sigpways,]
  
  subemap <- emapplot(subora, color='Percent_of_DEGs', repel=T,
                      #coords = subdat[,1:2],
                      showCategory = length(sigpways),
                      cex_label_category=0.4, 
                      cex_category = 0.3,cex_line = 0.3
  )+
    scale_fill_distiller(palette = 'Reds', direction = 1, name='Percent of DEGs')+
    ggtitle(relabel) +
    theme(plot.title = element_text(hjust = 0.5))
  
  
  subora <- setReadable(subora, 'org.Mm.eg.db', keyType = 'ENSEMBL')
  
  subcnet <- cnetplot(subora, showCategory = 5, foldChange = genes,
                      shadowtext='category',
                      cex_gene = 0.5, cex_label_gene = 0.5)
  
  
  
  # subora@result$qvalue <- subora@result$Percent_of_DEGs
  subddotplot <-  dotplot(subora, showCategory=30,font.size=10)+
    scale_color_gradient2(low = 'Firebrick', high = 'Steelblue', midpoint = 0.05)
  
  
  
  clusteroutdir_sub <- paste0(clustoutdir, '/', origlabel)
  dir.create(clusteroutdir_sub)
  
  ggsave(paste0(clusteroutdir_sub, '/emap.jpg'), subemap, width = 10, height = 10)
  ggsave(paste0(clusteroutdir_sub, '/cnet.jpg'), subcnet, width = 10, height = 10,bg = 'white')
  ggsave(paste0(clusteroutdir_sub, '/dotplot.jpg'), subddotplot, width = 10, height = 10)
  
  
}

#if needed, do any downstream stuff...
# ie, compare sub-clusters, special cnetpots...
sigclustidx <- 1

origlabel <- sigclusts[sigclustidx, 'clust']
relabel <- sigclusts[sigclustidx, 'relabel']

message(relabel)

#get the pathways
sigpways <- sigclustpways[[sigclustidx]]

#recompute...

subterm2gene <- term2gene[term2gene$gs_name %in% sigpways,]
subora <- enricher(genenames,
                   TERM2GENE = subterm2gene
                   # TERM2NAME = term2name,
                   # OrgDb = org.Mm.eg.db,
                   # keyType = 'ENSEMBL',
                   # pvalueCutoff = 0.001,
)


subora <- enrichplot::pairwise_termsim(subora)


#try to plot color = num genes in overexp gene list?
subora@result$Percent_of_DEGs <- (subora@result$Count / length(genenames))*100

#get dat with coords for this cluster
subdat <- dat[dat$name %in% sigpways,]

subemap <- emapplot(subora, color='Percent_of_DEGs', repel=T,
                    #coords = subdat[,1:2],
                    showCategory = length(sigpways),
                    cex_label_category=0.4, 
                    cex_category = 0.3,cex_line = 0.3
)+
  scale_fill_distiller(palette = 'Reds', direction = 1, name='Percent of DEGs')+
  ggtitle(relabel) +
  theme(plot.title = element_text(hjust = 0.5))


subora <- setReadable(subora, 'org.Mm.eg.db', keyType = 'ENSEMBL')

subcnet <- cnetplot(subora, showCategory = 5, foldChange = genes,
                    shadowtext='category',
                    cex_gene = 0.5, cex_label_gene = 0.5)

#for specific cnet...
pways_for_cnet <- c('DESCARTES_FETAL_INTESTINE_MYELOID_CELLS',
                    'DESCARTES_FETAL_ADRENAL_LYMPHOID_CELLS',
                    'GOBP_DIVALENT_INORGANIC_CATION_HOMEOSTASIS')


spec_cnet <- cnetplot(subora, showCategory = pways_for_cnet, foldChange = genes,
         shadowtext='category',
         cex_gene = 0.5, cex_label_gene = 0.5, cex_label_category=0.8)


clusteroutdir_sub <- paste0(clustoutdir, '/', origlabel)
# dir.create(clusteroutdir_sub)

ggsave(paste0(clusteroutdir_sub, '/speccnet.jpg'), spec_cnet, width = 7, height = 7,bg = 'white')









dev.off()






























### second, do dowreg genes

suboutdir <- paste0(outdir, '/underexp')
dir.create(suboutdir)



#select significant genes...
sig <- res[res$padj < 0.05,]
sig <- sig[sig$log2FoldChange < -1,]

#rank by FC
genes <- sig$log2FoldChange
names(genes) <- sig$ensembl_gene_id

#sort them by decreasing magnitude
genes <- genes[order(abs(genes), decreasing = T)]

# just in case, remove any NAs...
genes <- na.omit(genes)


#just give the actual gene names as input...
genenames <- names(genes)






## run the enrichment analysis ##

down_msigdb_ora <- enricher(genenames,
                          TERM2GENE = term2gene, 
                          # TERM2NAME = term2name,
                          # OrgDb = org.Mm.eg.db,
                          # keyType = 'ENSEMBL',
                          pvalueCutoff = 0.001
)


#run termsimilarity
down_msigdb_ora <- enrichplot::pairwise_termsim(down_msigdb_ora)

#emaapplot, with just the top 30...
down_emap_pval <- emapplot(down_msigdb_ora, layout='graphopt', repel=T)


#try to plot color = num genes in overexp gene list?
down_msigdb_ora@result$Percent_of_DEGs <- (down_msigdb_ora@result$Count / length(genenames))*100



down_emap_perc <- emapplot(down_msigdb_ora, color='Percent_of_DEGs', layout='graphopt', repel=T,
                         coords = down_emap_pval$data[,1:2])+
  scale_fill_distiller(palette = 'Reds', direction = 1, name='Percent of DEGs')


### for the real analysis:
# network with all pathways in top 200 by pval
# the top 200 is the default in the termsim...

set.seed(2021)
down_emap_total <- emapplot(down_msigdb_ora, color='Percent_of_DEGs', repel=T, showCategory = 200,  layout='graphopt',
                          cex_label_category=0.4, cex_category = 0.3,cex_line = 0.3,seed=2021)+
  scale_fill_distiller(palette = 'Reds', direction = 1, name='Percent of DEGs')


down_emap_total


# try to recover clusters....
# we need the object and the plot
obj <- down_msigdb_ora
plot <- down_emap_total


#get termsim and min_dist
mat <- obj@termsim
min_dist = 0.2 #default from package



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
                      intersect(term2gene[term2gene$gs_name ==pw,"ensembl_gene"]$ensembl_gene, 
                                sig$ensembl_gene_id)
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

clustpercs[clustpercs$perc>=10,]



reslist[names(reslist) %in% clustpercs[clustpercs$perc>=10,]$clust]

#get the genes for unified pathway analysis...
as.data.frame( res[match(clustgeneslist[[3]], res$ensembl_gene_id), "mgi_symbol"] )


clustpercs$relabel <- clustpercs$clust
clustpercs[clustpercs$perc>=10,'relabel'] <- c('Cluster_1_Hypoxia-and-Inflammation', 'Cluster_2_Cholesterol-related', 'Cluster_3_Mesenchymal')
clustpercs$relabel <- paste0(clustpercs$relabel, '\n', round(clustpercs$perc, 1), ' % of DEGs')



#the cumulative percentage is not relly helpful, 
# gene lists can have some minimal overlap (defined by min_dist)

#just check any above 10%?
# get percents for the top clusters
sigclusts <- clustpercs[clustpercs$perc>=10,]

# get actuaal genes of the top clusters
sigclustgenes <- clustgeneslist[names(clustgeneslist) %in% sigclusts$clust]
names(sigclustgenes) <- sigclusts$relabel

# get pathways in the top clusters
sigclustpways <- reslist[names(reslist) %in% sigclusts$clust]
names(sigclustpways) <- sigclusts$relabel



#ggvenndiagram of the top ones....
down_venn <- ggVennDiagram::ggVennDiagram(sigclustgenes)+
  labs(title = 'Venn Diagram of genes from top clusters of enriched pathways')+
  theme(plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(expand = expansion(mult = .2))


#try to plot...
# get sigclustpwasys as a df
sigclustpwaysdf <- bind_rows( lapply(1:length(sigclustpways), function(x){
  i=sigclustpways[[x]]
  name=names(sigclustpways)[x]
  data.frame(pways=i, clust=name)
})
)


#add color to result emaapplot
down_msigdb_ora@result$Cluster <- NA
down_msigdb_ora@result[match(sigclustpwaysdf$pways, down_msigdb_ora@result$ID),'Cluster'] <- sigclustpwaysdf$clust


set.seed(2021)
down_emap_total_redone <-  emapplot(down_msigdb_ora, color='Cluster', repel=T, showCategory = 200, layout='graphopt',
                                  node_label='None', #cex_label_category=0.4, 
                                  cex_category = 0.3,cex_line = 0.3 )+
  scale_fill_brewer(palette = 'Set2', direction = 1, name='Cluster')

#add ggrepel labels to emapplot
dat <- down_emap_total_redone$data
dat$Cluster <- NA
dat[match(sigclustpwaysdf$pways, dat$name),'Cluster'] <- sigclustpwaysdf$clust


repelaggr <- aggregate(x ~ Cluster, dat, mean)
repelaggr$y <- aggregate(y ~ Cluster, dat, mean)[,2]

#add percentage to the label name...

down_emap_total_redone_withrepel <- down_emap_total_redone +  ggrepel::geom_label_repel(inherit.aes = F, 
                                                                                    data = repelaggr, aes(x=x,y=y,fill=Cluster,label=Cluster), 
                                                                                    min.segment.length = Inf)




#plot just the sig clusters...
subdat <- dat[!is.na(dat$Cluster),]

# to be able to "fill" by cluster in the ggrepel labels,
# we explicitly set fill outside of aes with colors defined before plot call
repelaggr$Color <- RColorBrewer::brewer.pal(nrow(repelaggr), 'Set2')


down_emap_total_justclusters <- emapplot(down_msigdb_ora, color='Percent_of_DEGs', repel=T, showCategory = subdat$name, #layout='graphopt',
                                       cex_label_category=0.4, 
                                       cex_category = 0.3,cex_line = 0.3,
                                       coords = subdat[,1:2])+
  scale_fill_distiller(palette = 'Reds', direction = 1, name='Percent_of_DEGs')




down_emap_total_justclusters_withrepel <- down_emap_total_justclusters +
  ggrepel::geom_label_repel(inherit.aes = F,
                            data = repelaggr, aes(x=x,y=y,label=Cluster),
                            fill= repelaggr$Color, 
                            #box.padding = 10, max.overlaps = 200,
                            # direction = 'x',
                            min.segment.length = Inf)





down_emap_total_justclusters_colorbyclust <- emapplot(down_msigdb_ora, color='Cluster', repel=T, 
                                                    showCategory = subdat$name, #layout='graphopt',
                                                    cex_label_category=0.4, 
                                                    cex_category = 0.3,cex_line = 0.3,
                                                    # coords = subdat[,1:2]
)+
  scale_fill_brewer(palette = 'Set2', direction = 1, name='Cluster')



# down_emap_total_redone
# down_emap_total_redone_withrepel
# 
# down_emap_total / down_emap_total_redone

down_emap_total + down_emap_total_redone_withrepel + 
  down_emap_total_justclusters + down_emap_total_justclusters_withrepel



#plot the basic plots....
ggsave(paste0(suboutdir, '/emap_pval.jpg'), down_emap_pval, width= 10, height= 10)
ggsave(paste0(suboutdir, '/emap_perc.jpg'), down_emap_perc, width= 10, height= 10)

#olot the total plot...
ggsave(paste0(suboutdir, '/emap_total.jpg'), down_emap_total, width= 10, height= 10)


#plot the cluster plots
ggsave(paste0(suboutdir, '/emap_total_cluster.jpg'), down_emap_total_redone, width= 10, height= 10)
ggsave(paste0(suboutdir, '/emap_total_cluster_withrepel.jpg'), down_emap_total_redone_withrepel, width= 10, height= 10)


# plot the total plots with just the top cluster nodes...
ggsave(paste0(suboutdir, '/emap_total_justsig.jpg'), down_emap_total_justclusters, width= 10, height= 10)
ggsave(paste0(suboutdir, '/emap_total_justsig_withrepel.jpg'), down_emap_total_justclusters_withrepel, width= 10, height= 10)

ggsave(paste0(suboutdir, '/emap_total_justsig_colorbyclust.jpg'), down_emap_total_justclusters_colorbyclust, width= 10, height= 10)

#save the venn
ggsave(paste0(suboutdir, '/venn.jpg'), down_venn, width= 10, height= 10, bg='white')


#plot the nice combo
combo <- down_emap_total + down_emap_total_redone_withrepel + 
  down_emap_total_justclusters + down_emap_total_justclusters_withrepel

ggsave(paste0(suboutdir, '/emap_combo.jpg'), combo, width = 14, height = 10)






#try to make one big one....
options(ggrepel.max.overlaps = 4)

#add color to result emaapplot
down_msigdb_ora@result$Cluster <- 'Other'
down_msigdb_ora@result[match(sigclustpwaysdf$pways, down_msigdb_ora@result$ID),'Cluster'] <- sigclustpwaysdf$clust


#set.seed(2021)
levs <- table(down_msigdb_ora@result$Cluster)
pal <- RColorBrewer::brewer.pal(name = 'Accent', n = (length(levs) - 1))
pal <- c(pal, 'grey50')

emap_final <-  emapplot(down_msigdb_ora, color='Cluster', repel=F, showCategory = 200, # layout='graphopt',
                        node_label='None', #cex_label_category=0.4, 
                        cex_category = 0.5,cex_line = 0.5 )+
  scale_fill_manual(values = pal, name = 'Cluster')


ggsave(paste0(suboutdir, '/emap_final.jpg'), emap_final, width= 7, height= 7)


options(ggrepel.max.overlaps = 20)





### MAAKE SUB-PLOTS FOR THE TOP CLUSTERS...
clustoutdir <- paste0(suboutdir, '/pathway_clusters/')
dir.create(clustoutdir)


for(sigclustidx in c(1:nrow(sigclusts))){
  
  origlabel <- sigclusts[sigclustidx, 'clust']
  relabel <- sigclusts[sigclustidx, 'relabel']
  
  message(relabel)
  
  #get the pathways
  sigpways <- sigclustpways[[sigclustidx]]
  
  #recompute...
  
  subterm2gene <- term2gene[term2gene$gs_name %in% sigpways,]
  subora <- enricher(genenames,
                     TERM2GENE = subterm2gene
                     # TERM2NAME = term2name,
                     # OrgDb = org.Mm.eg.db,
                     # keyType = 'ENSEMBL',
                     # pvalueCutoff = 0.001,
  )
  
  
  subora <- enrichplot::pairwise_termsim(subora)
  
  
  #try to plot color = num genes in overexp gene list?
  subora@result$Percent_of_DEGs <- (subora@result$Count / length(genenames))*100
  
  #get dat with coords for this cluster
  subdat <- dat[dat$name %in% sigpways,]
  
  subemap <- emapplot(subora, color='Percent_of_DEGs', repel=T,
                      #coords = subdat[,1:2],
                      showCategory = length(sigpways),
                      cex_label_category=0.4, 
                      cex_category = 0.3,cex_line = 0.3
  )+
    scale_fill_distiller(palette = 'Reds', direction = 1, name='Percent of DEGs')+
    ggtitle(relabel) +
    theme(plot.title = element_text(hjust = 0.5))
  
  
  subora <- setReadable(subora, 'org.Mm.eg.db', keyType = 'ENSEMBL')
  
  subcnet <- cnetplot(subora, showCategory = 5, foldChange = genes,
                      shadowtext='category',
                      cex_gene = 0.5, cex_label_gene = 0.5)
  
  
  
  # subora@result$qvalue <- subora@result$Percent_of_DEGs
  subddotplot <-  dotplot(subora, showCategory=30,font.size=10)+
    scale_color_gradient2(low = 'Firebrick', high = 'Steelblue', midpoint = 0.05)
  
  
  
  clusteroutdir_sub <- paste0(clustoutdir, '/', origlabel)
  dir.create(clusteroutdir_sub)
  
  ggsave(paste0(clusteroutdir_sub, '/emap.jpg'), subemap, width = 10, height = 10)
  ggsave(paste0(clusteroutdir_sub, '/cnet.jpg'), subcnet, width = 10, height = 10,bg = 'white')
  ggsave(paste0(clusteroutdir_sub, '/dotplot.jpg'), subddotplot, width = 10, height = 10)
  
  
}

#if needed, do any downstream stuff...
# ie, compare sub-clusters, special cnetpots...
sigclustidx <- 1

origlabel <- sigclusts[sigclustidx, 'clust']
relabel <- sigclusts[sigclustidx, 'relabel']

message(relabel)

#get the pathways
sigpways <- sigclustpways[[sigclustidx]]

#recompute...

subterm2gene <- term2gene[term2gene$gs_name %in% sigpways,]
subora <- enricher(genenames,
                   TERM2GENE = subterm2gene
                   # TERM2NAME = term2name,
                   # OrgDb = org.Mm.eg.db,
                   # keyType = 'ENSEMBL',
                   # pvalueCutoff = 0.001,
)


subora <- enrichplot::pairwise_termsim(subora)


#try to plot color = num genes in overexp gene list?
subora@result$Percent_of_DEGs <- (subora@result$Count / length(genenames))*100

#get dat with coords for this cluster
subdat <- dat[dat$name %in% sigpways,]

subemap <- emapplot(subora, color='Percent_of_DEGs', repel=T,
                    #coords = subdat[,1:2],
                    showCategory = length(sigpways),
                    cex_label_category=0.4, 
                    cex_category = 0.3,cex_line = 0.3
)+
  scale_fill_distiller(palette = 'Reds', direction = 1, name='Percent of DEGs')+
  ggtitle(relabel) +
  theme(plot.title = element_text(hjust = 0.5))


subora <- setReadable(subora, 'org.Mm.eg.db', keyType = 'ENSEMBL')

subcnet <- cnetplot(subora, showCategory = 5, foldChange = genes,
                    shadowtext='category',
                    cex_gene = 0.5, cex_label_gene = 0.5)

#for specific cnet...
pways_for_cnet <- c('HALLMARK_HYPOXIA',
                    'SEKI_INFLAMMATORY_RESPONSE_LPS_UP',
                    'MTOR_UP.N4.V1_UP')


spec_cnet <- cnetplot(subora, showCategory = pways_for_cnet, foldChange = genes,
                      shadowtext='category',
                      cex_gene = 0.5, cex_label_gene = 0.5, cex_label_category=0.8)


clusteroutdir_sub <- paste0(clustoutdir, '/', origlabel)
# dir.create(clusteroutdir_sub)

ggsave(paste0(clusteroutdir_sub, '/speccnet.jpg'), spec_cnet, width = 7, height = 7,bg = 'white')









dev.off()




beepr::beep()