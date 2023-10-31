library(TCGAbiolinks)
library(tidyverse)
library(SummarizedExperiment)
library(maftools) # v‘2.2.10’


# setwd('~/Dropbox/data/bangdata/target/pancan/')
setwd('/home/Alex/data/pancan')
set.seed(2022)


#2022.may.24
packageVersion('TCGAbiolinks') #‘2.25.0’
packageVersion("tidyverse") #‘1.3.1'
packageVersion("SummarizedExperiment") #‘1.22.0’
packageVersion("maftools") #‘2.8.5’


#TCGA Biolinks has three steps
# 1 - query data, find it on the GDC portal
# 2 - actually download the data
# 3 - parse the data and read into R

#1 - query data

#queried on 2022.may.24
projects <- TCGAbiolinks:::getGDCprojects()$project_id

#keep only tcga and target data...
projects <- c( grep(projects, pattern = 'TCGA', value = T, ignore.case = T),
               grep(projects, pattern = 'TARGET', value = T, ignore.case = T)
)


project='TARGET-NBL'
project = 'TARGET-OS'

#download...

for(project in projects){
  
  #see if RNAseq is there, 
  deets <- TCGAbiolinks:::getProjectSummary(project)$data_categories
  
  if(project %in% list.files('GDCdata/')){
    
    message('\nAlready downloaded and skipping ', project)
    
    next()
    
  }
  
  
  if("Transcriptome Profiling" %in% deets$data_category){
    
    rdf <- data.frame(project = project,
                      downloaded = T)
    
    message('\n\n\nStarting ', project, '\n\n\n')
    
    
    tryCatch({
      #1 - query data
      query <- GDCquery(project = project,
                        data.category = "Transcriptome Profiling",
                        data.type = "Gene Expression Quantification", 
                        workflow.type = "STAR - Counts",
                        sample.type = "Primary Tumor")
      
      
      
      #2 - download data. downloads to current working directory
      GDCdownload(query)
      
      
    })
    
    #3 - parse and read in, as a RangedSummarizedExperiment-class object
    #rse <- GDCprepare(query = query)
    
    
    
  } 
  
  
}



# beepr::beep()















### above was done on threadripper

### parse and read in...?
library(TCGAbiolinks)
library(tidyverse)
library(SummarizedExperiment)
library(maftools) # v‘2.2.10’


# setwd('~/Dropbox/data/bangdata/target/pancan/')
setwd('/home/Alex/data/pancan')
set.seed(2022)


#2022.may.24
packageVersion('TCGAbiolinks') #‘2.25.0’
packageVersion("tidyverse") #‘1.3.1'
packageVersion("SummarizedExperiment") #‘1.22.0’
packageVersion("maftools") #‘2.8.5’


#TCGA Biolinks has three steps
# 1 - query data, find it on the GDC portal
# 2 - actually download the data
# 3 - parse the data and read into R

#1 - query data

#queried on 2022.may.24
projects <- TCGAbiolinks:::getGDCprojects()$project_id

#keep only tcga and target data...
projects <- c( grep(projects, pattern = 'TCGA', value = T, ignore.case = T),
               grep(projects, pattern = 'TARGET', value = T, ignore.case = T)
)



project = 'TARGET-OS'


#for now, exclude NBL, can't download properly
# https://github.com/BioinformaticsFMRP/TCGAbiolinks/issues/517
projects <- projects[projects != 'TARGET-NBL']



for(project in projects){
  
  
  
  
  #check if project is in gdc data...
  if( project %in% list.files('GDCdata/') ){
    
    
    #check if project is NOT ALREADY in parsed...
    savelistfilename <- paste0('parsed/', project, '.rds')
    
    if( !(savelistfilename %in% list.files('parsed', full.names = T))  ){
      
      
      
      message('\n\n\nStarting ', project, '\n\n\n')
      
      
      query <- NULL
      
      while(is.null(query) ){
        try(
          query <- GDCquery(project = project,
                            data.category = "Transcriptome Profiling",
                            data.type = "Gene Expression Quantification", 
                            workflow.type = "STAR - Counts",
                            sample.type = "Primary Tumor")
        )
        
      }
      
      
      #3 - parse and read in, as a RangedSummarizedExperiment-class object
      rse <- NULL
      
      while(is.null(rse) ){
        try(
          rse <- GDCprepare(query = query)
        )
        
      }
      
      
      ### prep to run ###
      
      #1 - get the gene expression matrix from the object
      gem_tpm <- assays(rse)$`tpm_unstrand`
      gem_counts <- assays(rse)$`unstranded`
      
      
      #2 - row ranges: stores the ensembl ID and the corresponding HUGO gene name
      geneIDs <- as.data.frame(rowRanges(rse))
      
      
      #3 - metadata, including sample IDs and some clinical / demographic data:
      md <- colData(rse)
      
      
      
      ### keep only protein-coding genes
      geneIDs <- geneIDs[geneIDs$gene_type=='protein_coding',]
      
      #remove dupicate symbols, apparently not many
      dups <- geneIDs$gene_name[duplicated(geneIDs$gene_name)]
      geneIDs <- geneIDs[!(geneIDs$gene_name %in% dups),]
      
      #subset gem
      gem_tpm <- gem_tpm[match(geneIDs$gene_id, rownames(gem_tpm)),] 
      gem_counts <- gem_counts[match(geneIDs$gene_id, rownames(gem_counts)),]
      
      #rename the gem
      rownames(gem_tpm) <- geneIDs$gene_name
      rownames(gem_counts) <- geneIDs$gene_name
      
      
      
      ### save everything nicely ###
      savelist <- list(gem_tpm=gem_tpm,
                       gem_counts=gem_counts,
                       md=md,
                       geneIDs=geneIDs)
      rm(gem_tpm,gem_counts,md,geneIDs, rse)
      
      savelistfilename <- paste0('parsed/', project, '.rds')
      saveRDS(savelist, savelistfilename)
      
      
      rm(gem,md)
      
    } #end if statement checking if project is NOT ALREADY in parsed
    
  } #end if statement checking if project is in GDCdata
  
  
  
} # end project loop

# beepr::beep()




