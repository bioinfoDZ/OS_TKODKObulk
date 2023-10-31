library(tidyverse)
library(preprocessCore)
library(survival)
library(survminer)
library(DESeq2)
library(patchwork)

library(biomaRt) #mouse-human homologs

library(MCPcounter) #for calculating microenv scores in human

set.seed(2020)


#### prep run #####

#biomart literally failed completely, wonderful

### as of 2022/05/05 we use biomart archive from feb2022, when we did the analysis to begin with...

hom <- read.csv('data/biomart_nodups_may05-2022_feb2021archive.csv')



# set up module score function

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



### fix long names function
fixlongnames <- function(pnames){
  
  pnewnames <- c()
  for(p in pnames){
    
    
    if(str_length(p) > 50){
      
      #try to find and replace underscore closest to halfway...
      halfway <- round(str_length(p)/2)
      underscorepos <- str_locate_all(p,'_')[[1]][,1]
      
      distances <- abs(halfway-underscorepos)
      
      which_und_is_closest <- which.min(distances)
      
      split <- str_split_fixed(p,'_',Inf)
      
      newp <- paste0(paste(split[1,1:which_und_is_closest], collapse = '_'),
                     '        ', '\n',
                     paste(split[1,(which_und_is_closest+1):ncol(split)], collapse = '_')
      )
      
      
      
    } else{newp <- p}
    
    pnewnames <- c(pnewnames, newp)
    
  }
  
  pnewnames
}














#tpm or counts
tpm = T


if(tpm == T){
  l <- readRDS('data/parsed/parsed_tpm.rds')
} else{
  l <- readRDS('data/parsed/parsed.rds')
}

#read the list:
# names(l)
# 1           gem
# 2 samp_case_ids
# 3 clinical_data
# 4      metadata
# 5  ensembl_hugo


gem  <- l[[1]]
sdf  <- l[[2]]
clin <- l[[3]]
sampmat <- l[[4]]
genes <- l[[5]]

rm(l, sampmat)



#fix clin

# binarize First Event status to censored / event (ie no / yes or 0/1)
firstevent <- clin$`First Event`
firstevent[is.na(firstevent)] <- 'NA'
firstevent[ firstevent == 'Censored' | firstevent == 'None' ] <- 'Censored'
firstevent[ firstevent!= 'NA' & firstevent != 'Censored' ] <- 'Event'
firstevent[firstevent == 'NA'] <- NA

clin$event_recoded <- factor(firstevent, levels = c('Event', 'Censored') )
rm(firstevent)

# recode time to first event
firsteventtime <- clin$`Event Free Survival Time in Days`
firsteventtime <- firsteventtime / 365.25

clin$event_recoded_time <- firsteventtime 
rm(firsteventtime)

# recode vital status to death 1 = yes, 0 = no
vitalstat <- clin$`Vital Status`
vitalstat[is.na(vitalstat)] <- 'NA'
vitalstat[ vitalstat == 'Alive'] <- 'Censored'
vitalstat[ vitalstat!= 'NA' & vitalstat != 'Censored' ] <- 'Dead'
vitalstat[vitalstat == 'NA'] <- NA

clin$vitalstatus_recoded <- factor(vitalstat, levels = c('Dead', 'Censored') )
rm(vitalstat)

# recode survival time
ostime <- clin$`Overall Survival Time in Days`
ostime <- ostime / 365.25

clin$vitalstatus_recoded_time <- ostime 
rm(ostime)


# recode metastasis, remove parenthesis etc
clin$MetastasisAtDiagnosis <- "Metastatic"
clin[grepl('Non', clin$`Disease at diagnosis`),"MetastasisAtDiagnosis"] <- "Non-metastatic"

clin$MetastasisAtDiagnosis <- factor(clin$MetastasisAtDiagnosis, levels = c('Non-metastatic', 'Metastatic'))

#recode age to years
clin$age_years <- clin$`Age at Diagnosis in Days` / 365.25



#set up cleaner clinical variable dataframe

clindat <- data.frame(Sample = clin$`TARGET USI`,
                      Gender = clin$Gender,
                      Age_at_diagnosis = clin$age_years,
                      
                      Mestastasis_at_diagnosis = clin$MetastasisAtDiagnosis,
                      
                      First_event = clin$event_recoded,
                      First_event_time = clin$event_recoded_time,
                      
                      Vital_status = clin$vitalstatus_recoded,
                      Vital_status_time = clin$vitalstatus_recoded_time
)








rm(clin)


#check NAs
clindat[!complete.cases(clindat),]

#one sample can be fixed: TARGET-40-0A4I9K
# marked as dead, but event is missing
# even tho event is missing, first event time is there
# vital status time and first event time are both >= 5
# we will set first event as true for this sample.
clindat[clindat$Sample == 'TARGET-40-0A4I9K',"First_event"] <- "Event"

#other samples are missing any survival, just remove them.

clindat <- clindat[complete.cases(clindat),]



#remove samples missing survival

sdf <- sdf[match(clindat$Sample, sdf$cases),]
gem <- gem[,match(sdf$samps, colnames(gem))]

# also rename gem colnames... from "sample" to "case"
colnames(gem) <- clindat$Sample
rm(sdf)




#remove extreme low TPM samples
badsamps <- c('TARGET-40-PALFYN',
              'TARGET-40-PASKZZ')

clindat <- clindat[!( clindat$Sample %in% badsamps),]
gem <- gem[,colnames(gem) %in% clindat$Sample]



#norm: quant norm for TPM, DESeq2 size factors if counts

if(tpm == T){
  
  #quantile norm; 
  gemt <- as.data.frame(preprocessCore::normalize.quantiles(as.matrix(gem)))
  rownames(gemt) <- rownames(gem)
  colnames(gemt) <- colnames(gem)
  
  gem <- gemt
  rm(gemt)
  
} else{
  
  gem <- t( t(gem) / DESeq2::estimateSizeFactorsForMatrix(gem) )
  
}




#### calculate scores ####

# do this before PCA, we want to see how they loook in PCA



#calculate the module scores: TKO

#read in DE genes
res <- read.csv('~/Dropbox/data/bangdata/2021october-TKOvsDKO/results/comparative-de/TKO-vs-DKO/2.deresults/DEresults-normalized_count_matrix.csv')

#get only sig res
res <- res[abs(res$log2FoldChange) > 1 & res$padj < 0.05,]


#get the human IDs...

#remove dup genes (if any)
dups <- names( table(res$mgi_symbol)[table(res$mgi_symbol)>1] )
res <- res[!(res$mgi_symbol %in% dups),]


#try to match up...
res <- res[res$mgi_symbol %in% hom$MGI.symbol,]
hom_in_res <- hom[match(res$mgi_symbol, hom$MGI.symbol),]
res <- cbind(hom_in_res$HGNC.symbol, res)
colnames(res)[1] <- 'HGNC_ortholog'

#split into up/down lists
up <- res[res$log2FoldChange>0,"HGNC_ortholog"]
down <- res[res$log2FoldChange<0,"HGNC_ortholog"]

mm <- data.frame(TKO_overexpressed = modulescore(gem, up),
                 TKO_underexpressed = modulescore(gem, down))

## save it for later myeloid in/leave out
tko_up <- up


rm(res, hom_in_res,up,down,dups)





### do it again with pevon, c1 treated...


## pevon
#read in DE genes
res <- read.csv('~/Dropbox/data/bangdata/may2021-DKOvsTKOvsSkp2inhib-bulk/results/celllines/comparative-de/DKO-culture-Pevon-vs-DKO-culture/2.deresults/DEresults-normalized_count_matrix.csv')

#get only sig res
res <- res[abs(res$log2FoldChange) > 1 & res$padj < 0.05,]


#get the human IDs...

#remove dup genes (if any)
dups <- names( table(res$mgi_symbol)[table(res$mgi_symbol)>1] )
res <- res[!(res$mgi_symbol %in% dups),]


#try to match up...
res <- res[res$mgi_symbol %in% hom$MGI.symbol,]
hom_in_res <- hom[match(res$mgi_symbol, hom$MGI.symbol),]
res <- cbind(hom_in_res$HGNC.symbol, res)
colnames(res)[1] <- 'HGNC_ortholog'

#split into up/down lists
up <- res[res$log2FoldChange>0,"HGNC_ortholog"]
down <- res[res$log2FoldChange<0,"HGNC_ortholog"]


mm_add <- data.frame(Pevon_treated_up = modulescore(gem, up),
                     Pevon_treated_down = modulescore(gem, down))

mm <- cbind(mm, mm_add)

rm(res, hom_in_res,up,down,dups, mm_add)







## C1
#read in DE genes
res <- read.csv('~/Dropbox/data/bangdata/may2021-DKOvsTKOvsSkp2inhib-bulk/results/celllines/comparative-de/DKO-culture-C1-vs-DKO-culture/2.deresults/DEresults-normalized_count_matrix.csv')

#get only sig res
res <- res[abs(res$log2FoldChange) > 1 & res$padj < 0.1,]


#get the human IDs...

#remove dup genes (if any)
dups <- names( table(res$mgi_symbol)[table(res$mgi_symbol)>1] )
res <- res[!(res$mgi_symbol %in% dups),]


#try to match up...
res <- res[res$mgi_symbol %in% hom$MGI.symbol,]
hom_in_res <- hom[match(res$mgi_symbol, hom$MGI.symbol),]
res <- cbind(hom_in_res$HGNC.symbol, res)
colnames(res)[1] <- 'HGNC_ortholog'

#split into up/down lists
up <- res[res$log2FoldChange>0,"HGNC_ortholog"]
down <- res[res$log2FoldChange<0,"HGNC_ortholog"]


mm_add <- data.frame(C1_treated_up = modulescore(gem, up),
                     C1_treated_down = modulescore(gem, down))

mm <- cbind(mm, mm_add)

rm(res, hom_in_res,up,down,dups, mm_add)









#MCP counter scores

# links to a live github: 2022.01.05
# mcp <- MCPcounter.estimate(gem, featuresType = 'HUGO_symbols')
# saveRDS(mcp, 'data/mcpresults_2022.rds')
mcp <- readRDS("data/mcpresults_2022.rds")

#check mcp gene markers...
# links to a live github: 2022.01.05
#genes=read.table(curl("http://raw.githubusercontent.com/ebecht/MCPcounter/master/Signatures/genes.txt"),sep="\t",stringsAsFactors=FALSE,header=TRUE,colClasses="character",check.names=FALSE)



#add to mm
rownames(mcp) <- gsub(' ', '_', rownames(mcp))
mm <- cbind(mm,t(mcp))



rm(mcp)










############ prep the special TKO gene list, minus myeloid #############

up <- tko_up

c8 <- msigdbr::msigdbr(category = 'C8')


#remove all myeloid pathways
myeloid <- c8[grepl("myeloid", ignore.case = T, c8$gs_name),]

up_nomyeloid <- up[!(up %in% myeloid$human_gene_symbol)]



# try keeping ONLY the myeloid....
#remove all myeloid pathways
myeloid <- c8[grepl("myeloid", ignore.case = T, c8$gs_name),]

up_onlymyeloid <- up[up %in% myeloid$human_gene_symbol]






#put in df

editedmodules <- data.frame(TKO_overexpressed_nomyeloid = modulescore(gem, up_nomyeloid),
                            TKO_overexpressed_onlymyeloid = modulescore(gem, up_onlymyeloid) )

#bind
mm <- cbind(mm,editedmodules)

rm(editedmodules)







### get SKP2 low patienrts

skp2 <- t(gem['SKP2',])


mm$SKP2 <- skp2


x = cor.test(mm$SKP2, mm$TKO_overexpressed)
ggplot(mm, aes(SKP2, TKO_overexpressed))+
  geom_point()+
  labs(subtitle = paste0('Cor = ', formatC(x$estimate, 3, format = "f"), '\nP = ', formatC(x$p.value, 3, format = "f")))

x = cor.test(mm$SKP2, mm$TKO_underexpressed)
ggplot(mm, aes(SKP2, TKO_underexpressed))+
  geom_point()+
  labs(subtitle = paste0('Cor = ', formatC(x$estimate, 3, format = "f"), '\nP = ', formatC(x$p.value, 3, format = "f")))



pdf('finalized_results/skp2_correlations/skp2_monocytelineaage_cor.pdf', 
    height = 4, width = 4)

x = cor.test(mm$SKP2, mm$Monocytic_lineage)
ggplot(mm, aes(SKP2, Monocytic_lineage))+
  geom_point()+
  labs(subtitle = paste0('Cor = ', formatC(x$estimate, 3, format = "f"), '\nP = ', formatC(x$p.value, 3, format = "f")))+
  ylab('MCP-Counter Monocyte Lineage Score') + xlab('SKP2 (TPM)')+
  theme_classic2()

dev.off()



x = cor.test(mm$SKP2, mm$CD8_T_cells)
ggplot(mm, aes(SKP2, CD8_T_cells))+
  geom_point()+
  labs(subtitle = paste0('Cor = ', formatC(x$estimate, 3, format = "f"), '\nP = ', formatC(x$p.value, 3, format = "f")))+
  ylab('MCP-Counter Monocyte Lineage Score') + xlab('SKP2 (TPM)')+
  theme_classic2()



ComplexHeatmap::Heatmap(cor(mm))






### find all skp2 cor genes


genes <- rownames(gem)
genes <- genes[!(genes %in% c('SKP2'))]

reslist <- list()
for(gene in genes){
  genevec <- t(gem[gene,])
  y <- cor.test(genevec[,1], skp2[,1])
  reslist[[gene]] <- data.frame(gene = gene, cor = y$estimate, p = y$p.value)
}

res <- dplyr::bind_rows(reslist)

res$fdr <- p.adjust(res$p, method = 'fdr')

res <- res[order(res$cor, decreasing = T),]



sigres <- res[res$fdr<0.05,]



outdir <- paste0('nonsurvival/skp2cor_tkoenrichment/')
write.csv(paste0(outdir, '/skp2_cor_res.csv'), x = res, row.names = F, quote = F)






### check dist of TKO up genes...

#calculate the module scores: TKO

#read in DE genes
res <- read.csv('~/Dropbox/data/bangdata/2021october-TKOvsDKO/results/comparative-de/TKO-vs-DKO/2.deresults/DEresults-normalized_count_matrix.csv')

#get only sig res
res <- res[abs(res$log2FoldChange) > 1 & res$padj < 0.05,]


#get the human IDs...

#remove dup genes (if any)
dups <- names( table(res$mgi_symbol)[table(res$mgi_symbol)>1] )
res <- res[!(res$mgi_symbol %in% dups),]


#try to match up...
res <- res[res$mgi_symbol %in% hom$MGI.symbol,]
hom_in_res <- hom[match(res$mgi_symbol, hom$MGI.symbol),]
res <- cbind(hom_in_res$HGNC.symbol, res)
colnames(res)[1] <- 'HGNC_ortholog'

#split into up/down lists
up <- res[res$log2FoldChange>0,"HGNC_ortholog"]
down <- res[res$log2FoldChange<0,"HGNC_ortholog"]

rm(res,hom_in_res, dups)


##read in cor genes
res <- read.csv('nonsurvival/skp2cor_tkoenrichment//skp2_cor_res.csv')


#up genes
res_tkoup <- res[res$gene %in% up,]

#dn genes
res_tkodown <- res[res$gene %in% down,]

up_hist <- ggplot(res_tkoup, aes(cor))+
  geom_histogram(fill='steelblue', col='black')+ggtitle('TKO Up genes')


dn_hist <- ggplot(res_tkodown, aes(cor))+
  geom_histogram(fill='steelblue', col='black')+ ggtitle('TKO Dn genes')



up_hist_sig <- ggplot(res_tkoup[res_tkoup$p < 0.05,], aes(cor))+
  geom_histogram(fill='steelblue', col='black')+ggtitle('TKO Up genes, sig cor only') 


dn_hist_sig <- ggplot(res_tkodown[res_tkodown$p < 0.05,], aes(cor))+
  geom_histogram(fill='steelblue', col='black')+ggtitle('TKO Dn genes, sig cor only')



up_hist + up_hist_sig + dn_hist + dn_hist_sig








#### enrichment of skp2cor genes...

sigres <- res[res$fdr < 0.05,]


##### first, do it for os cor genes...


#set up term2gene
term2gene <- data.frame(cluster = 'tko_up',
                        gene = up)

term2gene <- rbind(term2gene,
                   data.frame(cluster='tko_dn', gene = down))



#get universe genes
univ <- res$gene

#keep only pway gnees in the universe?
term2gene <- term2gene[term2gene$gene %in% univ,]
pwaylens <- table(term2gene$cluster)

#get my degs
mydegs <- sigres[sigres$cor>0, "gene"]


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



ggplot(fres, aes(ID, OR, col = -log10(P), size = overlapsize))+
  geom_point() + scale_color_gradient(low = 'grey', high = 'purple')


#keep the res for up cor genes
fres_upcorgenes <- fres









### repeat for dn cor genes



#get my degs
mydegs <- sigres[sigres$cor<0, "gene"]


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


ggplot(fres, aes(ID, OR, col = -log10(P), size = overlapsize))+
  geom_point() + scale_color_gradient(low = 'grey', high = 'purple')


#keep the res for up cor genes
fres_dncorgenes <- fres



### plot togehter...
fres_upcorgenes$CorList <- 'SKP2 positive\ncoexpressed genes'
fres_dncorgenes$CorList <- 'SKP2 negative\ncoexpressed genes'

fres_comb <- rbind(fres_upcorgenes, fres_dncorgenes)


ggplot(fres_comb, aes(ID, OR, col = -log10(P), size = overlapsize))+
  geom_point() + scale_color_gradient(low = 'grey', high = 'purple')+
  facet_wrap(~CorList)+
  theme_light()




### show enrichment of TKO up genes in negatively cor gene list

pdf('finalized_results/skp2_correlations/enrichmet_skp2negcor_tko.pdf',
    height = 4, width = 4)


fres$ID <- c('TKO down', 'TKO up')
ggplot(fres, aes(ID, OR, col = -log10(P), size = overlapsize))+
  geom_point() + scale_color_gradient(low = 'grey', high = 'purple')+
  theme_light()


dev.off()



## finalize, for table, match ensembl ID with gene name



#tpm or counts
tpm = T


if(tpm == T){
  l <- readRDS('data/parsed/parsed_tpm.rds')
} else{
  l <- readRDS('data/parsed/parsed.rds')
}

#read the list:
# names(l)
# 1           gem
# 2 samp_case_ids
# 3 clinical_data
# 4      metadata
# 5  ensembl_hugo


gem  <- l[[1]]
sdf  <- l[[2]]
clin <- l[[3]]
sampmat <- l[[4]]
genes <- l[[5]]

rm(l, sampmat)

rm(gem,sdf,clin,sampmat)



res <- read.csv('nonsurvival/skp2cor_tkoenrichment//skp2_cor_res.csv')



#match up genes and res
genes <- genes[genes$external_gene_name %in% res$gene,]

genes <- genes[match(res$gene, genes$external_gene_name),]

res <- cbind(genes, res)

res <- res[,c('ensembl_gene_id', 'external_gene_name', 'cor', 'p', 'fdr')]
colnames(res)[2] <- 'gene_symbol'

write.csv(res, '~/Dropbox/data/bangdata/2021october-TKOvsDKO/manuscript/manuscript/supplement_tables/supptab3_TARGET_SKP2_cor.csv',
          quote = F, row.names = F)

