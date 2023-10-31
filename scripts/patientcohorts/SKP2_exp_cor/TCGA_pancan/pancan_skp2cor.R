### parse and read in...?
library(TCGAbiolinks)
library(tidyverse)
library(SummarizedExperiment)
library(maftools) # v‘2.2.10’

library(survival)
library(survminer)


setwd('~/Dropbox/data/bangdata/target/pancan/')
# setwd('/home/Alex/data/pancan')
set.seed(2022)


#2022.may.24
packageVersion('TCGAbiolinks') #‘2.25.0’
packageVersion("tidyverse") #‘1.3.1'
packageVersion("SummarizedExperiment") #‘1.22.0’
packageVersion("maftools") #‘2.8.5’



### loop thru rds files

projects <- list.files('parsed/')
projects <- str_split_fixed(projects, '\\.', 2)[,1]

proj = "TARGET-OS"

genes = c('SKP2')

reslist <- list()
for(proj in projects){
  
  samp <- list.files('parsed', pattern = proj, full.names = T)
  
  message('Reading ', proj)
  
  #read in data list
  dl <- readRDS(samp)
  
  #get the data
  gem <- dl[[2]]
  md <- dl[[3]]
  
  
  ### calculate MCP
  
  # log counts...
  mcp <- MCPcounter::MCPcounter.estimate(log(gem+1), featuresType = 'HUGO_symbols')
  
  
  #prep to calculate cor with skp2 (TPM)
  gem <- dl[[1]]
  mcp <- as.data.frame(t(mcp))
  
  #add genes
  genes <- genes[genes %in% rownames(gem)]
  genedf <- data.frame( gem[match(genes, rownames(gem)),])
  colnames(genedf) <- genes
  mcp <- cbind(mcp, genedf)
  
  #correlate
  cr <- cor.test(mcp$SKP2, mcp$`Monocytic lineage`)
  
  cr <- data.frame(proj=proj, r=cr$estimate, p=cr$p.value, meanMac=mean(mcp$`Monocytic lineage`), meanSKP2=mean(mcp$SKP2))
  
  
  
  
  
  reslist[[proj]] <- cr
  
}

res <- dplyr::bind_rows(reslist)

res$fdr <- p.adjust(res$p)



ggplot(res, aes(r,-log10(p), label=proj))+
  geom_point()+
  ggrepel::geom_text_repel()+
  geom_vline(xintercept = 0, linetype='dotted')+
  geom_hline(yintercept = -log10(0.05), linetype='dotted')
# scale_x_continuous(limits = c(-1,1))+
#scale_color_gradient2(low='blue',high = 'red', mid = 'gray')



ggplot(res, aes(meanSKP2,meanMac, label=proj, col = r, size = -log10(p)))+
  geom_point()+
  ggrepel::geom_text_repel(size=2.5)+
  scale_color_gradient2(low='blue',high = 'red', mid = 'gray')











### check sarc ###
proj = "TCGA-SARC"

samp <- list.files('parsed', pattern = proj, full.names = T)

message('Reading ', proj)

#read in data list
dl <- readRDS(samp)

#get the data
gem <- dl[[2]]
md <- dl[[3]]



### ignore subbtypes for now...?
#subtypes are a mess between two columns, 'primary_diagnosis' and 'paper_histology'
# lets also add in ICD10 code...?
st <- as.data.frame(md[,c('barcode', 'primary_diagnosis', 'paper_histology' ,'tissue_or_organ_of_origin')])
# st[is.na(st$paper_histology),"paper_histology"] <- st[is.na(st$paper_histology),"primary_diagnosis"]

write.csv('~/Desktop/tcgasarchistology.csv', x=st, quote = T, row.names = F)






### calculate MCP

# log counts...
mcp <- MCPcounter::MCPcounter.estimate(log(gem+1), featuresType = 'HUGO_symbols')


#prep to calculate cor with skp2 (TPM)
gem <- dl[[1]]
mcp <- as.data.frame(t(mcp))

#add genes
genes <- genes[genes %in% rownames(gem)]
genedf <- data.frame( gem[match(genes, rownames(gem)),])
colnames(genedf) <- genes
mcp <- cbind(mcp, genedf)

#correlate
cr <- cor.test(mcp$SKP2, mcp$`Monocytic lineage`)

cr <- data.frame(proj=proj, r=cr$estimate, p=cr$p.value, meanMac=mean(mcp$`Monocytic lineage`), meanSKP2=mean(mcp$SKP2))





#### survival analysis ######



### loop thru rds files

projects <- list.files('parsed/')
projects <- str_split_fixed(projects, '\\.', 2)[,1]

proj = "TARGET-OS"

genes = c('SKP2')

survlist <- list()
for(proj in projects){
  
  samp <- list.files('parsed', pattern = proj, full.names = T)
  
  message('Reading ', proj)
  
  #read in data list
  dl <- readRDS(samp)
  
  #get the data
  gem <- dl[[1]]
  md <- dl[[3]]
  
  
  ### calculate MCP
  
  # log counts...
  mcp <- MCPcounter::MCPcounter.estimate(log(gem+1), featuresType = 'HUGO_symbols')
  
  
  #prep to calculate cor with skp2 (TPM)
  gem <- dl[[1]]
  mcp <- as.data.frame(t(mcp))
  
  #add genes
  genes <- genes[genes %in% rownames(gem)]
  genedf <- data.frame( gem[match(genes, rownames(gem)),])
  colnames(genedf) <- genes
  mcp <- cbind(mcp, genedf)
  
  
  #correlate
  cr <- cor.test(mcp$SKP2, mcp$`Monocytic lineage`)
  
  cr <- data.frame(proj=proj, r=cr$estimate, p=cr$p.value, meanMac=mean(mcp$`Monocytic lineage`), meanSKP2=mean(mcp$SKP2))
  
  
  
  
  ### survival analysis... SKp2 and monocyte....
  surv <- md[,c('barcode', 'vital_status', 'days_to_death', 'days_to_last_follow_up')]
  
  #remove missing...
  surv <- surv[!(is.na(surv$vital_status)),]
  
  
  #make survival time: days to die for vital_status="dead", days to followup for alive
  surv$os_time <- surv$days_to_death
  
  surv[surv$vital_status=='Alive', "os_time"] <- surv[surv$vital_status=='Alive', "days_to_last_follow_up"]
  
  #exclude those not alive or dead... ie unknown...
  surv <- surv[,c('barcode','vital_status','os_time')]
  surv <- surv[complete.cases(surv),]
  
  # make sure same samples are in gem;mcp
  gem <- gem[,match(surv$barcode, colnames(gem))]
  mcp <- mcp[match(surv$barcode, rownames(mcp)),]
  
  
  #prep surv, mtch prior analysis exactly...
  clin <- surv
  clin[clin$vital_status=='Alive', "vital_status"] <- 'Censored'
  clin$vital_status <- factor( clin$vital_status, levels = c('Dead', 'Censored') )
  clin$vitalstatus_recoded <- factor(clin$vital_status, levels = c('Dead', 'Censored') )
  
  clin$vitalstatus_recoded_time <- clin$os_time / 365.25
  
  clindat <- clin
  rm(clin,surv)
  
  
  #keep only monocyte, genes...
  keep <- c('Monocytic lineage', genes)
  mm <- mcp[,keep]
  
  
  
  #lets use vital status
  status <- 'vitalstatus_recoded'
  time <- 'vitalstatus_recoded_time'
  
  
  
  #survival analysis
  scorenames <- names(mm)
  
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
    df <- data.frame(samps = clindat$barcode,
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

    #bind all outputs
    models <- list(dich = dich,
                   terc = terc,
                   continuous = continuous
    )
    
    plots <- list(gdich, gterc, distplots)
    
    #temporarily save outs to plot list...
    plotlist[[scorename]] <- plots
    modellist[[scorename]] <- models
    datalist[[scorename]] <- df
    
    
    # setTxtProgressBar(pb, scoreidx)
    
  }
  
  
  ### just get the HRs for now...
  
  modelreslist <- list()
  for(scoreidx in 1:length(scorenames)) {
    
    #get the score name
    scorename <- colnames(mm)[scoreidx]
    
    models <- modellist[[scoreidx]]
    
    univar <- models[[3]]
    
    univar <- summary(univar)
    
    modeldf <- data.frame(proj = proj,
                          scorename = scorename,
                          hr = univar$coefficients[2], 
                          lo = univar$conf.int[3],
                          hi = univar$conf.int[4],
                          p = univar$waldtest[3],
                          row.names = NULL
                          )
    
    modelreslist[[scorename]] <- modeldf
  }
  
  
  modeldf <- dplyr::bind_rows(modelreslist)
  
  survlist[[proj]] <- modeldf
  
  
}

survdf <- dplyr::bind_rows(survlist)



#just monocyte
mon <- survdf[survdf$scorename=='Monocytic lineage',]


res$hr <- mon$hr
res$hr_p <- mon$p


res$goodbad <- 'bad'
res[res$hr<1,"goodbad"] <- 'good'


ggplot(res, aes(r,-log10(p), label=proj, col = hr, size = -log10(hr_p), shape=goodbad))+
  geom_point()+
  ggrepel::geom_text_repel(size=2.5)+
  geom_vline(xintercept = 0, linetype='dotted')+
  geom_hline(yintercept = -log10(0.05), linetype='dotted')+
  scale_color_gradient2(low='blue',high = 'red', mid = 'gray', midpoint = 1)



