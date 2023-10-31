library(tidyverse)
library(preprocessCore)
library(survival)
library(survminer)
library(DESeq2)
library(patchwork)

library(biomaRt) #mouse-human homologs

library(MCPcounter) #for calculating microenv scores in human

set.seed(2020)



### final followup from tko vs dko bulk rnaseq paper

#Multivar survival model 
# with both monocyte MCP score AND 
# TKO overexp gene list


#Subtract myeloid genes from TKO: 
# remove DE genes present in msigdb myeloid pathways, 
# check survival








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




#cap at 5 years
# set time above 5 to just 5
# also we need to censor after 5 years before doing that

clindat[clindat$First_event_time >5,'First_event'] <- 'Censored'
clindat[clindat$First_event_time >5,'First_event_time'] <- 5


clindat[clindat$Vital_status_time >5,'Vital_status'] <- 'Censored'
clindat[clindat$Vital_status_time >5,'Vital_status_time'] <- 5




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




#clean up
rm(c8,hom,myeloid,tko_up,up,up_nomyeloid,up_onlymyeloid, genes)





#lets use vital status
status <- 'Vital_status'
time <- 'Vital_status_time'



#set outdir
outdir <- 'survivalanalysis/FINAL_TKO_VS_DKO/final_target_tkodkobulkpaper_5year'
dir.create(outdir)







### editied modules

suboutdir <- paste0(outdir, '/modulescores_vitalstatus')
dir.create(suboutdir)

f#get module score names
scorenames <- names(mm)


total = ncol(mm)
pb <- txtProgressBar(min = 0, max = total, style = 3)


plotlist <- list()
modellist <- list()
datalist <- list() 

for(scoreidx in 1:length(scorenames)) {
  
  #get the score name
  scorename <- colnames(mm)[scoreidx]
  
  #get the score valeues, numeric continuopus
  scorevec <- mm[,scorename]
  
  ## dichotomize the score ##
  highname <- paste0('High\nN=', table(na.omit(scorevec) >= median(na.omit(scorevec)) )[2] )
  lowname <- paste0('Low\nN=', table(na.omit(scorevec) >= median(na.omit(scorevec)) )[1] )
  
  dichscore <- ifelse(scorevec >= median(na.omit(scorevec)), yes = highname, no = lowname)
  dichscore <- factor(dichscore, levels = c(lowname, highname))
  
  ### split score to tertiles ###
  terciles <- factor(ntile(scorevec, 3))
  terciles <- plyr::mapvalues(x = terciles, from = levels(terciles),
                              to = c( paste0('Low\nN=', table(terciles)[1]),
                                      paste0('Med\nN=', table(terciles)[2]),
                                      paste0('High\nN=', table(terciles)[2]) )
  )
  
  
  
  ### set up survival df ###
  df <- data.frame(samps = clindat$Sample,
                   continuous = scorevec,
                   dichotomized = dichscore,
                   terciles = terciles,
                   surv = clindat[,time],
                   status = clindat[,status])
  
  
  
  #check dist
  #check dist over dead/alive
  hist <- ggplot(df, aes(continuous, fill=status))+
    geom_histogram(col='black',position='identity', alpha=0.5) +
    facet_wrap(~status, nrow=2)
  
  #check dist cor with survival?
  
  corCensored <- cor.test(df[df$status=='Censored', "surv"], df[df$status=='Censored', "continuous"])
  corNonCensored <- cor.test(df[df$status!='Censored', "surv"], df[df$status!='Censored', "continuous"])
  corlab <- paste0('Non-Censored cor: ', round(corNonCensored$estimate, 3), ', Non-Censored P:', round(corNonCensored$p.value, 3),
                   '\nCensored cor: ', round(corCensored$estimate, 3), ', Censored P: ', round(corCensored$p.value, 3))
  corp <- ggplot(df, aes(continuous, surv, col = status))+
    geom_point()+
    geom_smooth()+
    facet_wrap(~status, nrow=2)+
    labs(caption = corlab)
  
  distplots <- hist+corp
  
  # code status as numeric
  df$status_code <- ifelse(df$status == 'Censored', yes=0, no = 1)
  
  
  
  ### tertile analysis ###
  fit <- survfit(Surv(surv, status_code) ~ terciles, data = df)
  
  gterc <- ggsurvplot(fit,
                      conf.int = F,
                      pval = TRUE, 
                      pval.method = T, 
                      linetype = "strata",  
                      test.for.trend = F,
                      title=scorename,
                      ggtheme = theme_classic2(), 
                      palette =  c( "#8DA0CB", "#66C2A5", "#FC8D62" ) )+ 
    xlab('Time (Years)')
  
  
  ### dichotomized analysis ###
  fit <- survfit(Surv(surv, status_code) ~ dichotomized, data = df)
  
  gdich <- ggsurvplot(fit,
                      conf.int = F,
                      pval = TRUE, #pval.coord = c(0, 0.15),
                      pval.method = T, #pval.method.coord = c(0, 0.25),
                      linetype = "strata",  
                      test.for.trend = F,
                      title=scorename,
                      ggtheme = theme_classic2(), 
                      palette = c("#377EB8" , "#E41A1C") )+ 
    xlab('Time (Years)')
  
  
  gterc <- gterc$plot
  gdich <- gdich$plot
  
  ### cox modelling ###
  #dich model, will match plot
  dich <- coxph(Surv(surv, status_code) ~ dichotomized, data = df)
  
  #terc
  terc <- coxph(Surv(surv, status_code) ~ terciles, data = df)
  
  #continuous univariate cox with cibersort score
  continuous <- coxph(Surv(surv, status_code) ~ continuous, data = df)
  
  #prep for multivar
  df <- cbind(df, clindat)
  
  #multivar cox test 
  multivar <- coxph(Surv(surv, status_code) ~ continuous + Mestastasis_at_diagnosis , data = df)
  
  #bind all outputs
  models <- list(dich = dich,
                 terc = terc,
                 continuous = continuous,
                 multivar = multivar
  )
  
  plots <- list(gdich, gterc, distplots)
  
  #temporarily save outs to plot list...
  plotlist[[scorename]] <- plots
  modellist[[scorename]] <- models
  datalist[[scorename]] <- df
  
  
  setTxtProgressBar(pb, scoreidx)
  
}




### check the hazard ratios ###


#univar first
univarlist <- list()
for(score in names(plotlist) ){
  summx <- summary(modellist[[score]]$continuous)
  univar <- as.data.frame(summx$conf.int)
  univar$p <- summx$coefficients[,5]
  colnames(univar) <- c('HR', 'hrneg', 'lower95', 'upper95', 'p')
  univar<- data.frame(score=score,univar, row.names = NULL)
  univarlist[[score]]  <- univar
}
univardf <- dplyr::bind_rows(univarlist)

univardf$score_rename <- fixlongnames(univardf$score)


#if more than 30, keep only 30 most significant on each side
# no filtering for genes...
univardf <- univardf[order(univardf$p),]
# up <- univardf[univardf$HR>1,] ; down <- univardf[univardf$HR<1,]
# up <- up[up$p < 0.05,] ; down <- down[down$p < 0.05,]
# 
# univardf <- rbind(up,down)



#sort by HR
univardf <- univardf[order(univardf$HR),]

#if more than 30, keep only 30 strongest...
# if(nrow(univardf) > 30){
#   univardf <- rbind( head(univardf, 15), tail(univardf, 15))
# }

#significance...
univardf$test <- ifelse(univardf$p < 0.05, 'significant', 'non-sig')

#factorize..
univardf$score_rename <- factor(univardf$score_rename, levels=univardf$score_rename)

forestunivar <- ggplot(univardf, aes(y = score_rename, x = HR, xmin=lower95, xmax=upper95, col = test))+
  geom_pointrange()+
  geom_vline(xintercept = 1, linetype = 'dashed')+
  theme_linedraw()+
  theme(axis.text.y = element_text(size = 7),
        axis.title.y = element_blank() )+
  labs(title = 'HR Forest plot, univariate continuous')+
  xlab('Hazard Ratio\nLeft = better survival ; Right = worse survival')






#hazard ratios for dich

### check the hazard ratios ###
dichlist <- list()
for(score in names(plotlist) ){
  summx <- summary(modellist[[score]]$dich)
  dich <- as.data.frame(summx$conf.int)
  dich$p <- summx$coefficients[,5]
  colnames(dich) <- c('HR', 'hrneg', 'lower95', 'upper95', 'p')
  dich<- data.frame(score=score,dich, row.names = NULL)
  dichlist[[score]]  <- dich
}
dichdf <- dplyr::bind_rows(dichlist)

dichdf$score_rename <- fixlongnames(dichdf$score)


#if more than 30, keep only 30 most significant on each side


# if(nrow(dichdf>30)){
#   dichdf <- dichdf[order(dichdf$p),]
#   up <- dichdf[dichdf$HR>1,] ; down <- dichdf[dichdf$HR<1,]
#   up <- up[up$p < 0.05,] ; down <- down[down$p < 0.05,]
#   
#   dichdf <- rbind(up,down)
# }


#sort by HR
dichdf <- dichdf[order(dichdf$HR),]

#if more than 30, keep only 30 strongest...
# if(nrow(dichdf) > 30){
#   dichdf <- rbind( head(dichdf, 15), tail(dichdf, 15))
# }


#significance...
dichdf$test <- ifelse(dichdf$p < 0.05, 'significant', 'non-sig')

#factorize..
dichdf$score_rename <- factor(dichdf$score_rename, levels=dichdf$score_rename)

forestdich <- ggplot(dichdf, aes(y = score_rename, x = HR, xmin=lower95, xmax=upper95, col=test))+
  geom_pointrange()+
  geom_vline(xintercept = 1, linetype = 'dashed')+
  theme_linedraw()+
  theme(axis.text.y = element_text(size = 7),
        axis.title.y = element_blank() )+
  labs(title = 'HR Forest plot, dichotomized')+
  xlab('Hazard Ratio\nLeft = better survival ; Right = worse survival')




#using significant pathways, save only worthwhile results...

#sigpways <- unique(c(univardf$score, dichdf$score))

#add any pathways we really want to see
# selectedpways <- c('OS_Stemness')
# for(i in selectedpways){
#   #check if in names and if not in sigpways
#   if(i %in% scorenames & !(i %in% sigpways)){
#     sigpways <- c(sigpways, i)
#   }
# }

sigpways <- scorenames

message('\n- Saving significant results')

sigresdir <- paste0(suboutdir, '/pathway-results-sigonly')
dir.create(sigresdir)

total = length(sigpways)
pb <- txtProgressBar(min = 0, max = total, style = 3)


#for significant results, make outputs
pdfplotlist <- list()
for(scoreidx in 1:length(sigpways) ){
  
  
  
  score <- sigpways[scoreidx]
  #message(score)
  
  #get models
  models <- modellist[[score]]
  df <- datalist[[score]] #load data into env or else, the diagnostics break
  
  #get plots
  plots <- plotlist[[score]]
  
  
  #get the coefficients, hazard ratios, confidence intervals, and P values
  summarytable <- bind_rows( lapply(models, function(x){as.data.frame( summary(x)$coefficients ) }) )
  summarytable <- bind_cols(summarytable,
                            bind_rows( lapply(models, function(x){as.data.frame( summary(x)$conf.int ) }) )[,3:4] )
  
  
  
  
  
  #only plot if significant...
  # ignore met
  # pvals <- summarytable$`Pr(>|z|)`[-(nrow(summarytable))]
  # if( !any(pvals < 0.10) ){
  #   next
  # }
  # this is defined by sigpways, not necessary
  
  #rename the continuous
  rownames <- rownames(summarytable)
  rownames[grepl('^continuous', rownames)] <- c('Univariable\nContinuous', 'Multivariable\nContinuous')
  rownames(summarytable) <- rownames
  
  #put conf int next to coef
  summarytable <- summarytable[,c(1,6,7,2,3,4,5)]
  
  #hr, lower/upperCI, P, and sig only...
  summarytable <- summarytable[,c(4,2,3,7)]
  
  
  #add significance marks
  summarytable$significance <- ''
  summarytable[summarytable$`Pr(>|z|)` < 0.1, "significance"] <- '.'
  summarytable[summarytable$`Pr(>|z|)` < 0.05, "significance"] <- '*'
  summarytable[summarytable$`Pr(>|z|)` < 0.01, "significance"] <- '**'
  summarytable[summarytable$`Pr(>|z|)` < 0.001, "significance"] <- '***'
  
  #round table
  summarytable[,-ncol(summarytable)] <- round( summarytable[,-ncol(summarytable)],3 )
  
  # plot it all together
  updatedplot <- (plots[[1]] + plots[[2]]) / ( gridExtra::tableGrob(summarytable)) + patchwork::plot_annotation(score) + 
    plot_layout(heights = c(1,2))
  
  
  pdfplotlist[[score]] <- updatedplot
  
  
  
  #save all outputs...
  resultdir <- paste0(sigresdir, '/', score)
  dir.create(resultdir)
  
  
  # saveRDS( file = paste0(resultdir, '/models.rds'), models  )
  
  
  #check model assumptions...
  modeldiagnosticsfile <- paste0(resultdir, '/modeldiagnostics.pdf')
  pdf(modeldiagnosticsfile, height = 10, width = 10)
  for(model in models){
    
    print( ggcoxzph( cox.zph(model) ) )
    
  }
  dev.off()
  
  
  #save plots
  
  
  gdich <- plots[[1]]
  gterc <- plots[[2]]
  distplots <- plots[[3]]
  suppressMessages(ggsave( paste0(resultdir, '/dist.jpg'), distplots , height = 7, width = 7))
  
  ggsave( paste0(resultdir, '/gdich.jpg'), gdich , height = 7, width = 7)
  ggsave( paste0(resultdir, '/gterc.jpg'), gterc , height = 7, width = 7)
  
  setTxtProgressBar(pb, scoreidx)
  
}




summaryfile <- paste0(suboutdir, '/summary-plots_models.pdf')
pdf(file = summaryfile, height = 10,width = 9)

print(forestunivar)
print(forestdich)
suppressMessages(print(pdfplotlist))

dev.off()







































### repeat with met free ###

### repeat with met free ###

### repeat with met free ###

### repeat with met free ###

### repeat with met free ###

### repeat with met free ###

### repeat with met free ###





#lets use MET
status <- 'First_event'
time <- 'First_event_time'



## suboutdir is defined above, gsub vital stat to metfree
suboutdir <- gsub('_vitalstatus', '_eventfree', suboutdir)

dir.create(suboutdir)


#get module score names
scorenames <- names(mm)




total = ncol(mm)
pb <- txtProgressBar(min = 0, max = total, style = 3)


plotlist <- list()
modellist <- list()
datalist <- list() 

for(scoreidx in 1:length(scorenames)) {
  
  #get the score name
  scorename <- colnames(mm)[scoreidx]
  
  #get the score valeues, numeric continuopus
  scorevec <- mm[,scorename]
  
  ## dichotomize the score ##
  highname <- paste0('High\nN=', table(na.omit(scorevec) >= median(na.omit(scorevec)) )[2] )
  lowname <- paste0('Low\nN=', table(na.omit(scorevec) >= median(na.omit(scorevec)) )[1] )
  
  dichscore <- ifelse(scorevec >= median(na.omit(scorevec)), yes = highname, no = lowname)
  dichscore <- factor(dichscore, levels = c(lowname, highname))
  
  ### split score to tertiles ###
  terciles <- factor(ntile(scorevec, 3))
  terciles <- plyr::mapvalues(x = terciles, from = levels(terciles),
                              to = c( paste0('Low\nN=', table(terciles)[1]),
                                      paste0('Med\nN=', table(terciles)[2]),
                                      paste0('High\nN=', table(terciles)[2]) )
  )
  
  
  
  ### set up survival df ###
  df <- data.frame(samps = clindat$Sample,
                   continuous = scorevec,
                   dichotomized = dichscore,
                   terciles = terciles,
                   surv = clindat[,time],
                   status = clindat[,status])
  
  
  
  #check dist
  #check dist over dead/alive
  hist <- ggplot(df, aes(continuous, fill=status))+
    geom_histogram(col='black',position='identity', alpha=0.5) +
    facet_wrap(~status, nrow=2)
  
  #check dist cor with survival?
  
  corCensored <- cor.test(df[df$status=='Censored', "surv"], df[df$status=='Censored', "continuous"])
  corNonCensored <- cor.test(df[df$status!='Censored', "surv"], df[df$status!='Censored', "continuous"])
  corlab <- paste0('Non-Censored cor: ', round(corNonCensored$estimate, 3), ', Non-Censored P:', round(corNonCensored$p.value, 3),
                   '\nCensored cor: ', round(corCensored$estimate, 3), ', Censored P: ', round(corCensored$p.value, 3))
  corp <- ggplot(df, aes(continuous, surv, col = status))+
    geom_point()+
    geom_smooth()+
    facet_wrap(~status, nrow=2)+
    labs(caption = corlab)
  
  distplots <- hist+corp
  
  # code status as numeric
  df$status_code <- ifelse(df$status == 'Censored', yes=0, no = 1)
  
  
  
  ### tertile analysis ###
  fit <- survfit(Surv(surv, status_code) ~ terciles, data = df)
  
  gterc <- ggsurvplot(fit,
                      conf.int = F,
                      pval = TRUE, 
                      pval.method = T, 
                      linetype = "strata",  
                      test.for.trend = F,
                      title=scorename,
                      ggtheme = theme_classic2(), 
                      palette =  c( "#8DA0CB", "#66C2A5", "#FC8D62" ) )+ 
    xlab('Time (Years)')
  
  
  ### dichotomized analysis ###
  fit <- survfit(Surv(surv, status_code) ~ dichotomized, data = df)
  
  gdich <- ggsurvplot(fit,
                      conf.int = F,
                      pval = TRUE, #pval.coord = c(0, 0.15),
                      pval.method = T, #pval.method.coord = c(0, 0.25),
                      linetype = "strata",  
                      test.for.trend = F,
                      title=scorename,
                      ggtheme = theme_classic2(), 
                      palette = c("#377EB8" , "#E41A1C") )+ 
    xlab('Time (Years)')
  
  
  gterc <- gterc$plot
  gdich <- gdich$plot
  
  ### cox modelling ###
  #dich model, will match plot
  dich <- coxph(Surv(surv, status_code) ~ dichotomized, data = df)
  
  #terc
  terc <- coxph(Surv(surv, status_code) ~ terciles, data = df)
  
  #continuous univariate cox with cibersort score
  continuous <- coxph(Surv(surv, status_code) ~ continuous, data = df)
  
  #prep for multivar
  df <- cbind(df, clindat)
  
  #multivar cox test 
  multivar <- coxph(Surv(surv, status_code) ~ continuous + Mestastasis_at_diagnosis , data = df)
  
  #bind all outputs
  models <- list(dich = dich,
                 terc = terc,
                 continuous = continuous,
                 multivar = multivar
  )
  
  plots <- list(gdich, gterc, distplots)
  
  #temporarily save outs to plot list...
  plotlist[[scorename]] <- plots
  modellist[[scorename]] <- models
  datalist[[scorename]] <- df
  
  
  setTxtProgressBar(pb, scoreidx)
  
}




### check the hazard ratios ###


#univar first
univarlist <- list()
for(score in names(plotlist) ){
  summx <- summary(modellist[[score]]$continuous)
  univar <- as.data.frame(summx$conf.int)
  univar$p <- summx$coefficients[,5]
  colnames(univar) <- c('HR', 'hrneg', 'lower95', 'upper95', 'p')
  univar<- data.frame(score=score,univar, row.names = NULL)
  univarlist[[score]]  <- univar
}
univardf <- dplyr::bind_rows(univarlist)

univardf$score_rename <- fixlongnames(univardf$score)


#if more than 30, keep only 30 most significant on each side
# no filtering for genes...
univardf <- univardf[order(univardf$p),]
# up <- univardf[univardf$HR>1,] ; down <- univardf[univardf$HR<1,]
# up <- up[up$p < 0.05,] ; down <- down[down$p < 0.05,]
# 
# univardf <- rbind(up,down)



#sort by HR
univardf <- univardf[order(univardf$HR),]

#if more than 30, keep only 30 strongest...
# if(nrow(univardf) > 30){
#   univardf <- rbind( head(univardf, 15), tail(univardf, 15))
# }

#significance...
univardf$test <- ifelse(univardf$p < 0.05, 'significant', 'non-sig')

#factorize..
univardf$score_rename <- factor(univardf$score_rename, levels=univardf$score_rename)

forestunivar <- ggplot(univardf, aes(y = score_rename, x = HR, xmin=lower95, xmax=upper95, col = test))+
  geom_pointrange()+
  geom_vline(xintercept = 1, linetype = 'dashed')+
  theme_linedraw()+
  theme(axis.text.y = element_text(size = 7),
        axis.title.y = element_blank() )+
  labs(title = 'HR Forest plot, univariate continuous')+
  xlab('Hazard Ratio\nLeft = better survival ; Right = worse survival')






#hazard ratios for dich

### check the hazard ratios ###
dichlist <- list()
for(score in names(plotlist) ){
  summx <- summary(modellist[[score]]$dich)
  dich <- as.data.frame(summx$conf.int)
  dich$p <- summx$coefficients[,5]
  colnames(dich) <- c('HR', 'hrneg', 'lower95', 'upper95', 'p')
  dich<- data.frame(score=score,dich, row.names = NULL)
  dichlist[[score]]  <- dich
}
dichdf <- dplyr::bind_rows(dichlist)

dichdf$score_rename <- fixlongnames(dichdf$score)


#if more than 30, keep only 30 most significant on each side


# if(nrow(dichdf>30)){
#   dichdf <- dichdf[order(dichdf$p),]
#   up <- dichdf[dichdf$HR>1,] ; down <- dichdf[dichdf$HR<1,]
#   up <- up[up$p < 0.05,] ; down <- down[down$p < 0.05,]
#   
#   dichdf <- rbind(up,down)
# }


#sort by HR
dichdf <- dichdf[order(dichdf$HR),]

#if more than 30, keep only 30 strongest...
# if(nrow(dichdf) > 30){
#   dichdf <- rbind( head(dichdf, 15), tail(dichdf, 15))
# }


#significance...
dichdf$test <- ifelse(dichdf$p < 0.05, 'significant', 'non-sig')

#factorize..
dichdf$score_rename <- factor(dichdf$score_rename, levels=dichdf$score_rename)

forestdich <- ggplot(dichdf, aes(y = score_rename, x = HR, xmin=lower95, xmax=upper95, col=test))+
  geom_pointrange()+
  geom_vline(xintercept = 1, linetype = 'dashed')+
  theme_linedraw()+
  theme(axis.text.y = element_text(size = 7),
        axis.title.y = element_blank() )+
  labs(title = 'HR Forest plot, dichotomized')+
  xlab('Hazard Ratio\nLeft = better survival ; Right = worse survival')




#using significant pathways, save only worthwhile results...

#sigpways <- unique(c(univardf$score, dichdf$score))

#add any pathways we really want to see
# selectedpways <- c('OS_Stemness')
# for(i in selectedpways){
#   #check if in names and if not in sigpways
#   if(i %in% scorenames & !(i %in% sigpways)){
#     sigpways <- c(sigpways, i)
#   }
# }

sigpways <- scorenames

message('\n- Saving significant results')

sigresdir <- paste0(suboutdir, '/pathway-results-sigonly')
dir.create(sigresdir)

total = length(sigpways)
pb <- txtProgressBar(min = 0, max = total, style = 3)


#for significant results, make outputs
pdfplotlist <- list()
for(scoreidx in 1:length(sigpways) ){
  
  
  
  score <- sigpways[scoreidx]
  #message(score)
  
  #get models
  models <- modellist[[score]]
  df <- datalist[[score]] #load data into env or else, the diagnostics break
  
  #get plots
  plots <- plotlist[[score]]
  
  
  #get the coefficients, hazard ratios, confidence intervals, and P values
  summarytable <- bind_rows( lapply(models, function(x){as.data.frame( summary(x)$coefficients ) }) )
  summarytable <- bind_cols(summarytable,
                            bind_rows( lapply(models, function(x){as.data.frame( summary(x)$conf.int ) }) )[,3:4] )
  
  
  
  
  
  #only plot if significant...
  # ignore met
  # pvals <- summarytable$`Pr(>|z|)`[-(nrow(summarytable))]
  # if( !any(pvals < 0.10) ){
  #   next
  # }
  # this is defined by sigpways, not necessary
  
  #rename the continuous
  rownames <- rownames(summarytable)
  rownames[grepl('^continuous', rownames)] <- c('Univariable\nContinuous', 'Multivariable\nContinuous')
  rownames(summarytable) <- rownames
  
  #put conf int next to coef
  summarytable <- summarytable[,c(1,6,7,2,3,4,5)]
  
  #hr, lower/upperCI, P, and sig only...
  summarytable <- summarytable[,c(4,2,3,7)]
  
  
  #add significance marks
  summarytable$significance <- ''
  summarytable[summarytable$`Pr(>|z|)` < 0.1, "significance"] <- '.'
  summarytable[summarytable$`Pr(>|z|)` < 0.05, "significance"] <- '*'
  summarytable[summarytable$`Pr(>|z|)` < 0.01, "significance"] <- '**'
  summarytable[summarytable$`Pr(>|z|)` < 0.001, "significance"] <- '***'
  
  #round table
  summarytable[,-ncol(summarytable)] <- round( summarytable[,-ncol(summarytable)],3 )
  
  # plot it all together
  updatedplot <- (plots[[1]] + plots[[2]]) / ( gridExtra::tableGrob(summarytable)) + patchwork::plot_annotation(score) + 
    plot_layout(heights = c(1,2))
  
  
  pdfplotlist[[score]] <- updatedplot
  
  
  
  #save all outputs...
  resultdir <- paste0(sigresdir, '/', score)
  dir.create(resultdir)
  
  
  # saveRDS( file = paste0(resultdir, '/models.rds'), models  )
  
  
  #check model assumptions...
  modeldiagnosticsfile <- paste0(resultdir, '/modeldiagnostics.pdf')
  pdf(modeldiagnosticsfile, height = 10, width = 10)
  for(model in models){
    
    print( ggcoxzph( cox.zph(model) ) )
    
  }
  dev.off()
  
  
  #save plots
  
  
  gdich <- plots[[1]]
  gterc <- plots[[2]]
  distplots <- plots[[3]]
  suppressMessages(ggsave( paste0(resultdir, '/dist.jpg'), distplots , height = 7, width = 7))
  
  ggsave( paste0(resultdir, '/gdich.jpg'), gdich , height = 7, width = 7)
  ggsave( paste0(resultdir, '/gterc.jpg'), gterc , height = 7, width = 7)
  
  setTxtProgressBar(pb, scoreidx)
  
}




summaryfile <- paste0(suboutdir, '/summary-plots_models.pdf')
pdf(file = summaryfile, height = 10,width = 9)

print(forestunivar)
print(forestdich)
suppressMessages(print(pdfplotlist))

dev.off()







beepr::beep()

