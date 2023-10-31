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







#readin gem, make first col row names
gem <- read.csv('data/kuijjer/parsed-gem.csv')
rownames(gem) <- gem[,1] ; gem <- gem[,-1]



#read in kuijjer md
md <- read.csv('data/kuijjer/parsed-md.csv')



#read in genes
genes <- read.csv('data/kuijjer/parsed-genesymbols.csv')


### prep them ###


#deal with GEM: average the multi-probe genes
#first, get gem with non-dup genes and hugo rownames

dups <- names( table(genes$hugo)[table(genes$hugo)>1] )
nondupgenes <- genes[!(genes$hugo %in% dups),]

savegem <- gem[match(nondupgenes$probeset, rownames(gem)), ]
rownames(savegem) <- nondupgenes$hugo

#next, get the gem with avg only genes
dupgenes <- genes[genes$hugo %in% dups,]

dupgem <- gem[match(dupgenes$probeset, rownames(gem)),]

#for each gene, get avg exp of each multi-probe
# ie, for each gene, get the rows, and take the colMean
total = length(dups)
pb <- txtProgressBar(min = 0, max = total, style = 3)


genemeans <- list()
for(geneidx in 1:length(dups)) {
  gene <- dups[geneidx]
  sdg <- dupgenes[dupgenes$hugo==gene,]
  sg <- dupgem[rownames(dupgem) %in% sdg$probeset,]
  cm <- data.frame(colMeans(sg))
  colnames(cm) <- gene
  genemeans[[gene]] <- cm
  rm(cm, sdg, sg)
  
  setTxtProgressBar(pb, geneidx)
}

genemeansdf <- dplyr::bind_cols(genemeans)
genemeansdf <- t(genemeansdf)

#rbind the averaged multi-probe gene matrix to the unique probe gene matrix
savegem <- rbind(savegem, genemeansdf)



#get md with only biopsies and remove the one mising TOF
# also, subset gem...

#one person not deceased a nd no met is missing tim, just mark is as followup=0...
md[is.na(md$tof),'tof'] <- 0

savemd <- md[md$type=='biopsy',]
savemd <- savemd[!is.na(savemd$tof),]  


savegem <- savegem[,match(savemd$sample, colnames(savegem))]


gem <- savegem
md <- savemd



#clean env
rm(dupgem, dupgenes, genemeans, genemeansdf, genes, nondupgenes,savegem, savemd, dups)
rm(gene,geneidx,total,pb)



### fix up clinical variables ###

##OVERALL SURVIVAL:
# recode vital status to death 1 = yes, 0 = no
vitalstat <- md$deceased
vitalstat[ vitalstat == 'false'] <- 'Censored'
vitalstat[ vitalstat != 'Censored' ] <- 'Dead'

md$vitalstatus_recoded <- factor(vitalstat, levels = c('Dead', 'Censored') )
rm(vitalstat)

# recode survival time
ostime <- md$tof

#make overall surv months to years
ostime <- ostime / 12

md$vitalstatus_recoded_time <- ostime 
rm(ostime)

###MET-FREE SURVIVAL:
# for met=true, use the tom1
# for met = false, use tof

#set time to met = tof, this is for met=false
md$time_to_met <- md$tof

#for met=true, set time_to_met=tom1
md[md$met=='true',"time_to_met"] <- md[md$met=='true',"tom1"]

#make time to met months to years
md$time_to_met <- md$time_to_met / 12

#set up actualy status var
md$met_recode <- 'Censored'
md[md$metastasis=='true', "met_recode"] <- 'Event'
md$met_recode <- factor(md$met_recode, levels = c('Event', 'Censored'))



#most important clin var is huvos score
# dichotomize, 1/2 vs 3/4
huvos <- md$huvos
huvos[!is.na(huvos) & huvos<=2] <- 'Huvos 1/2'
huvos[!is.na(huvos) & huvos!='Huvos 1/2'] <- 'Huvos 3/4'

huvos <- factor(huvos, levels = c('Huvos 1/2', 'Huvos 3/4'))
md$huvos_recode <- huvos; rm(huvos)


#recode age to years
md$age_years <- md$age / 12


#set up cleaner clinical variable dataframe

clindat <- data.frame(Sample = md$sample,
                      
                      Huvos = md$huvos_recode,
                      
                      Met_recode = md$met_recode,
                      Time_to_met = md$time_to_met,
                      
                      
                      Vital_status = md$vitalstatus_recoded,
                      Vital_status_time = md$vitalstatus_recoded_time
)





#check NAs --> they are missing huvos...
clindat[!complete.cases(clindat),]


#norm: quant norm for TPM, DESeq2 size factors if counts


#quantile norm; 
gemt <- as.data.frame(preprocessCore::normalize.quantiles(as.matrix(gem)))
rownames(gemt) <- rownames(gem)
colnames(gemt) <- colnames(gem)

gem <- gemt
rm(gemt)





#perform batch correction, no log though
gem <- as.matrix(gem)

gem_limma <- limma::removeBatchEffect(log2(gem+1), batch = md$lab)

gem_limma <- as.matrix(gem_limma)

gem <- gem_limma
rm(gem_limma)








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
mcp <- MCPcounter.estimate(gem, featuresType = 'HUGO_symbols')
# saveRDS(mcp, 'data/mcpresults_2022.rds')
# mcp <- readRDS("data/mcpresults_2022.rds")

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

skp2 <- gem['SKP2',]


mm$SKP2 <- skp2


x = cor.test(mm$SKP2, mm$TKO_overexpressed)
ggplot(mm, aes(SKP2, TKO_overexpressed))+
  geom_point()+
  labs(subtitle = paste0('Cor = ', formatC(x$estimate, 3, format = "f"), '\nP = ', formatC(x$p.value, 3, format = "f")))

x = cor.test(mm$SKP2, mm$TKO_underexpressed)
ggplot(mm, aes(SKP2, TKO_underexpressed))+
  geom_point()+
  labs(subtitle = paste0('Cor = ', formatC(x$estimate, 3, format = "f"), '\nP = ', formatC(x$p.value, 3, format = "f")))



pdf('finalized_results/skp2_correlations/skp2_monocytelineaage_cor_KUIJJER.pdf', 
    height = 4, width = 4)

x = cor.test(mm$SKP2, mm$Monocytic_lineage)
ggplot(mm, aes(SKP2, Monocytic_lineage))+
  geom_point()+
  labs(subtitle = paste0('Cor = ', formatC(x$estimate, 3, format = "f"), '\nP = ', formatC(x$p.value, 2, format = "e")))+
  ylab('MCP-Counter Monocyte Lineage Score') + xlab('SKP2 (TPM)')+
  theme_classic2()

dev.off()






### find all skp2 cor genes


genes <- rownames(gem)
genes <- genes[!(genes %in% c('SKP2'))]

reslist <- list()
for(gene in genes){
  genevec <- gem[gene,]
  y <- cor.test(genevec, skp2)
  reslist[[gene]] <- data.frame(gene = gene, cor = y$estimate, p = y$p.value)
}

res <- dplyr::bind_rows(reslist)

res$fdr <- p.adjust(res$p, method = 'fdr')

res <- res[order(res$cor, decreasing = T),]



sigres <- res[res$fdr<0.05,]



outdir <- paste0('nonsurvival/skp2cor_tkoenrichment/')
write.csv(paste0(outdir, '/skp2_cor_res_KUIJJER.csv'), x = res, row.names = F, quote = F)


