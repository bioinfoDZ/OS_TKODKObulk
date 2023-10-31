library(tidyverse)
library(clusterProfiler)
library(msigdbr)
library(enrichplot)
library(org.Mm.eg.db)

set.seed(2021)




## gsea for UPR and IFN-G, leading edge heatmaps
#read de res
save <- read.csv("results/comparative-de/TKO-vs-DKO/2.deresults/DEresults-normalized_count_matrix.csv")

#remove ID with dot
save$ensembl_gene_id <- str_split_fixed(save$ensembl_gene_id , pattern = '\\.', n = 2)[,1]

#make rownames as endembl ID
rownames(save) <- save$ensembl_gene_id 

#get pways
pathways <- msigdbr(species = 'Mus musculus')

cats <- unique(pathways$gs_cat)
cats <- str_sort(cats, numeric = T)




#run hallmarks

cat <- 'H'

cp <- pathways[pathways$gs_cat == cat,]


message('\nRunning ', cat)


#split df to list
cp  <- cp %>% split(x = .$ensembl_gene, f = .$gs_name)


#run GSEA and make main dotplot

names(cp) <- gsub('HALLMARK_', '', x = names(cp))
names(cp) <- gsub('_', ' ', x = names(cp))
names(cp) <- str_to_title(names(cp))

gseares <- FerrenaBulkRNAseq::gsea.results(save, pathways = cp)
gseadotplot <- FerrenaBulkRNAseq::gsea.dotplot.onecol(gseares = gseares[[1]], pathwayfontsize = 12) + labs(caption = '')


pdf("results/comparative-de/TKO-vs-DKO/3.gsea/H/BIGHALLMARKGSEA_updated.pdf", height = 8, width = 6)
gseadotplot
dev.off()


#save leading edge
le <- gseares[[2]]


### make GSEAplots for UPR and IFN

#get scores
tmp <- save
scores <- log(tmp$pvalue)
if (any(scores == -Inf)) {
  numuf = length(scores[scores == -Inf])
  adder = rev(1:numuf)
  for (uf in rev(1:numuf)) {
    scores[uf] <- scores[numuf + 1] + (adder[uf] * 
                                         -1)
  }
}
scores <- (-1 * scores) * sign(tmp$log2FoldChange)
names(scores) <- rownames(tmp)
rm(tmp)
scores <- sort(scores, decreasing = T)


res <- gseares[[1]]

res[res$pathway %in% c("Interferon Gamma Response", "Unfolded Protein Response"),]


pwayname <- "Interferon Gamma Response"
pway <- cp[[pwayname]]
ifngsea <- fgsea::plotEnrichment(pway, stats = scores) + labs(title = pwayname)


pwayname <- "Unfolded Protein Response"
pway <- cp[[pwayname]]
uprgsea <- fgsea::plotEnrichment(pway, stats = scores) + labs(title = pwayname)

pdf("results/comparative-de/TKO-vs-DKO/3.gsea/H/UPR_IFN_GSEAplots.pdf", width = 6, height = 4)
uprgsea
ifngsea
dev.off()






### leading egde heatmaps

#read de res
save <- read.csv("results/comparative-de/TKO-vs-DKO/2.deresults/DEresults-normalized_count_matrix.csv")

#remove ID with dot
save$ensembl_gene_id <- str_split_fixed(save$ensembl_gene_id , pattern = '\\.', n = 2)[,1]

#make rownames as endembl ID
rownames(save) <- save$ensembl_gene_id 

#rename columns
code <- data.frame(Sample = c("TL724M0T", "TL30F1T",  "TL744F2T", "DL761M2R", "DL579FOR", "DL608F1R"),
                   Code = c('TKO1', 'TKO2', 'TKO3', 'DKO1', 'DKO2', 'DKO3'))


colnames(save)[match(code$Sample, colnames(save))] <- code$Code




#adjust code so DKO first
code <- code[order(code$Code),]

#add in metadata
metadata <- as.data.frame(readxl::read_excel("data/metadata.xlsx"))
metadata <- metadata[match(code$Sample, metadata$Sample),]
metadata$Code <- code$Code


#get leading edge
le <- gseares[[2]]
le <- le[names(le) %in% c("Interferon Gamma Response", "Unfolded Protein Response")]


lehm <- lapply(1:length(le), function(i){
  
  #get genes and pway
  pway <- names(le)[i]
  ids <- le[[i]]
  message(pway)
  
  
  #remove any dups
  ids <- unique(ids)
  #make sure all in res
  ids <- ids[ids %in% save$ensembl_gene_id]
  
  #subset de res
  savesub <- save[match(ids, save$ensembl_gene_id),]
  
  
  mat <- savesub[match(code$Code, colnames(savesub))]
  rownames(mat) <- savesub$mgi_symbol
  
  #log and scale mat
  mat <- log1p(mat)
  mat <- t(scale(t(mat)))
  
  
  ### adjust 1-mar name.... ugh
  if('1-Mar' %in% rownames(mat)){
    rownames(mat)[rownames(mat) == '1-Mar'] <- "Marchf1"
  }
  
  
  #annotation df, has sample (column) name and color
  tmpgem <- mat
  
  annotdf <- data.frame(sample = colnames(tmpgem), 
                        sampcheck = metadata[match(colnames(tmpgem),  metadata$Code), "Sample"], 
                        cond = metadata[match(colnames(tmpgem), metadata$Code), "Condition"], 
                        color = metadata[match(colnames(tmpgem), metadata$Code), "Color"])
  
  hacol <- list(Condition = annotdf[, 4])
  names(hacol[[1]]) <- annotdf[, 3]
  
  #create the annotation object
  ha <- ComplexHeatmap::HeatmapAnnotation(Condition = annotdf[,3], col = hacol, name = 'Condition')
  
  hm <- ComplexHeatmap::Heatmap(mat, name = 'Scaled\nLog1p\nNorm.\nCounts', 
                                cluster_rows = F, cluster_columns = F,
                                column_title = pway,
                                top_annotation = ha)
  
  return(hm)
  
  
})

pdf("results/comparative-de/TKO-vs-DKO/3.gsea/H/UPR_IFN_leadingedge_heatmaps.pdf", height = 7, width = 5)
lehm
dev.off()



