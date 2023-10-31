library(tidyverse)
library(clusterProfiler)
library(msigdbr)
library(enrichplot)
library(org.Mm.eg.db)

set.seed(2021)



### check osteoclast vs monocyte ###


## calculate module score of osteoclast diff from TPM matrix

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




#get TPM
md <- as.data.frame(readxl::read_excel("data/metadata.xlsx"))



### read in TPM file ###
tpmfiles <- list.files("data/tpm/", full.names = T, recursive = T)

tpmfiles <- sapply(md$Sample, function(samp){tpmfiles[grepl(samp,tpmfiles)]})
names(tpmfiles) <- md$Code



### check ptprc ###
ptprc <- lapply(1:length(tpmfiles), function(i){
  
  tf <- tpmfiles[i]
  code <- names(tpmfiles)[i]
  
  ts <- read.table(tf, sep = '\t', header = T, fill = T)
  
  ts[ts$Gene_Id == 'ENSMUSG00000033016.16','TPM']
})

names(ptprc) <- names(tpmfiles)
ptprc



tpm <- lapply(1:length(tpmfiles), function(i){
  
  tf <- tpmfiles[i]
  code <- names(tpmfiles)[i]
  
  ts <- read.table(tf, sep = '\t', header = T, fill = T)
  
  
  #no NAs
  ts <- ts[complete.cases(ts),]
  
  #no duplicate IDs
  dups <- ts$Gene_Id[duplicated(ts$Gene_Id)]
  ts <- ts[!(ts$Gene_Id %in% dups),]
  
  tcol <- as.data.frame(ts$TPM)
  colnames(tcol) = code
  rownames(tcol) = ts$Gene_Id
  
  
  return(tcol)
  
})

names(tpm) <- names(tpmfiles)

#get shared genes...
sharedgenes <- Reduce(intersect, lapply(tpm, FUN = rownames))

#get all genes...
allgenes <- unique(unlist(sapply(tpm, rownames, simplify = T)))

#get missing genes
missing <- allgenes[!(allgenes %in% sharedgenes)]

#add missing genes as zeros
tpm <- lapply(tpm, function(tcol){
  
  missing_thissamp <- missing[!(missing %in% rownames(tcol))]
  missing_thissamp <- data.frame(rep(0, length(missing_thissamp)), 
                                 row.names = missing_thissamp)
  colnames(missing_thissamp) <- colnames(tcol)
  
  tcol <- rbind(tcol, missing_thissamp)
  
})


#match each with allgenes
tpm <- lapply(tpm, function(tcol){
  tcol <- tcol[match(allgenes, rownames(tcol)),,drop=F]
  tcol
})

tpm <- dplyr::bind_cols(tpm)



### convert ensembl to symbol

gencode <- read.csv("data/gencode.vM23.ids_names_types.csv")

#convert ensembl to symbol
gencode <- gencode[gencode$gene_id %in% rownames(tpm),]
gencode <- gencode[match(rownames(tpm), gencode$gene_id),]
# rownames(tpm) <- gencode$gene_name

#remove duplicate symbols
dups <- gencode$gene_name[duplicated(gencode$gene_name)]
gencode <- gencode[!(gencode$gene_name %in% dups),]
tpm <- tpm[match(gencode$gene_id, rownames(tpm)),]
rownames(tpm) <- gencode$gene_name

#use codes
colnames(tpm) <- md$Code
















## UPDATED ANALYSIS, OLD IS BELOW ##

# OSTEOCLAST / OSTEOBLAST: CHECK MARKERS FROM SCRNASEQ
os_mouse_markers <- read.csv("~/Dropbox/data/bangdata/scrnaseq-TKO-DKOAA-DKO/results_finalized/integration_res/risc_fromraw_samp4ref_RISC1.6/markers/Seurat_allmarkers_intcelltype.csv")

oscl <- os_mouse_markers[os_mouse_markers$cluster=='Osteoclast',]


#ref markers from stroma dataset
ref_stroma_markers <- readRDS('~/Dropbox/Result_from_Alex/deyoudata/stromadata_Jan2022/ref/triple/m_reference.rds')

osbl <- ref_stroma_markers[ref_stroma_markers$cluster=='Osteo',]

osbl$cluster <- 'Osteoblast'

#remove ctsk from osteoblast..
osbl <- osbl[osbl$gene != 'Ctsk',]


m <- rbind(oscl,osbl)




### top 20 genes ###
m <- m[m$gene %in% rownames(tpm),]
n <- 20
top <- m %>% group_by(cluster) %>% top_n(n = n, wt = avg_log2FC)
# top <- rbind(top, m)



mat <- tpm[match(top$gene, rownames(tpm)),]


mat <- log1p(mat)
mat <- t(scale(t(mat)))
# mat <- scale(mat)


#just tko
columnsplitnames <- md$Condition
# mat <- mat[,grepl('TKO', colnames(mat))]
# columnsplitnames <- rep('TKO', 3)


hm_osclbl <- ComplexHeatmap::Heatmap(mat, 
                                     # row_names_side = 'left',
                                     # row_title_gp = grid::gpar(fontsize = 5),
                                     width = ncol(mat)*unit(20, "mm"),
                                     name = 'Scaled\nTPM',
                                     column_title = "Osteoclast and Osteoblast markers",
                                     row_split = top$cluster,
                                     column_split = columnsplitnames
)
hm_osclbl














## try m0 m1 m2 markers


#two papers say tnfA is high in m0; though other papers consder it a marker of m1.
# pdac m0 paper
# https://pubmed.ncbi.nlm.nih.gov/33267818/
# hepatocellular carcinoma cell line m1 vs m2 etoposide paper
# https://bmccancer.biomedcentral.com/articles/10.1186/s12885-015-1546-9


# also , cd68 is used as pga-differentiated m0 macrophage (as opposed to nondiff monocyte) in HCC paper; 
# cd14 and cd16 (Fcgr3) are used as "mo monocyte" in pdac paper

#add in monocyte markers from scrnaseq OS paper 2:  https://www.frontiersin.org/journals/oncology/articles/10.3389/fonc.2021.709210/full

#m0?
m0 <- 'Tnf, Cd68, Cd14, Fcgr3, Cxcl3, Cxcl2'



#m1
# from Pdac paper, HCC paper, nd this review fig 1:
# https://www.frontiersin.org/journals/oncology/articles/10.3389/fonc.2019.01512/full

# also from OS scRNAseq paper 2:  https://www.frontiersin.org/journals/oncology/articles/10.3389/fonc.2021.709210/full

m1 <- 'Il6, Il1b, Cxcl10, Cd80, Cd86, Nos2, H2-Eb1, Axl, Ccl2, Cxcl9, Ifit1'

#m2 
# from Pdac paper, HCC paper, nd this review fig 1:
# https://www.frontiersin.org/journals/oncology/articles/10.3389/fonc.2019.01512/full

# also from OS scRNAseq paper 2:  https://www.frontiersin.org/journals/oncology/articles/10.3389/fonc.2021.709210/full


m2 <- 'Msr1, Mrc1, Vegf1, Maf, Fn1, Ccl3, Il10, Gpnmb, Stab1, F13a1, Mcm5, Fabp5, Txnip, Folr2, Mertk, Scart1'


## conflicting markers
# cd163 is conflicting, mostly m2, pdac paper shows higher in m1
# scpaper, it is enriched in the "inflammatory mac"s Ifit1+



m0 <- str_split_fixed(m0, pattern = ', ', n = Inf)[1,]
m1 <- str_split_fixed(m1, pattern = ', ', n = Inf)[1,]
m2 <- str_split_fixed(m2, pattern = ', ', n = Inf)[1,]



m0m1m2 <- data.frame(gene = c(m0,m1,m2),
                     group = c(
                       rep('M0', length(m0)),
                       rep('M1', length(m1)),
                       rep('M2', length(m2))
                     ))




m0m1m2 <- m0m1m2[m0m1m2$gene %in% rownames(tpm),]

mat <- tpm[match(m0m1m2$gene, rownames(tpm)),]

### try removing low DKO hi genes ###
rem <- as.data.frame(t(mat))

rem$Code <- md$Condition

avgs <- aggregate(. ~ Code, rem, mean)

mat <- mat[avgs[1,-1] < 5,]

group <- m0m1m2
group <- group[group$gene %in% rownames(mat),]


mat <- log1p(mat)
# mat <- t(scale(t(mat)))
# mat <- scale(mat)

hm_m0m1m2 <- ComplexHeatmap::Heatmap(mat, 
                                     # row_names_side = 'left',
                                     # row_title_gp = grid::gpar(fontsize = 5),
                                     width = ncol(mat)*unit(20, "mm"),
                                     name = 'Log1p\nTPM',
                                     column_title = "M0 vs M1 vs M2",
                                     row_split = group$group,
                                     column_split = md$Condition
)
hm_m0m1m2



m0m1m2 <- m0m1m2[m0m1m2$gene %in% rownames(mat),]





### now try geom mean of these



types <- unique(m0m1m2$group)
sl <- lapply(types, function(mt){
  genes <- m0m1m2[m0m1m2$group==mt,'gene']
  
  mat <- tpm[rownames(tpm) %in% genes,]
  
  #first log
  mat <- log1p(mat)
  
  #next get sample-wise mean
  ss <- Matrix::colMeans(mat)
  
  outdf <- as.data.frame(t(data.frame(score = ss, row.names = names(ss))))
  rownames(outdf) <- mt
  
  return(outdf)
  
})


scores = dplyr::bind_rows(sl)


#make barplot
ggs <- reshape2::melt(t(scores))

colnames(ggs) <- c('Code', 'Type', 'Score')


ggplot(ggs, aes(Code, Score, fill = Type))+
  geom_bar(stat="identity")

ggplot(ggs, aes(Code, Score, fill = Type))+
  geom_bar(stat="identity", position = 'fill')+
  scale_y_continuous(labels = scales::percent)



### get means
sm <- as.data.frame(t(scores))
sm$Condition <- md$Condition

avgs <- aggregate(. ~ Condition, sm, mean)




### check if these genes are DE
m0m1m2
deres <- read.csv('results/comparative-de/TKO-vs-DKO/2.deresults/DEresults-normalized_count_matrix.csv')


deres <- deres[match(m0m1m2$gene, deres$mgi_symbol),]

m0m1m2$log2FoldChange = deres$log2FoldChange
m0m1m2$pvalue = deres$pvalue
m0m1m2$padj = deres$padj



deres <- read.csv('results/comparative-de/TKO-vs-DKO/2.deresults/DEresults-normalized_count_matrix.csv')
deres <- deres[,-1]

deres <- deres[!duplicated(deres$mgi_symbol),]
rownames(deres) <- deres$mgi_symbol

# vp <- FerrenaBulkRNAseq::volcanoplot(deres, colors = unique(md$Color))
# 
# vpdat <- vp$data
# minires <- vpdat[match(m0m1m2$gene, vpdat$mgi_symbol),]
# 
# vp + geom_point(data = minires, size = 30, col = 'black', inherit.aes = F,
#                 aes(log2FoldChange, neglogtenp))



deres$neglogtenp <- -log10(deres$pvalue)
minires <- deres[match(m0m1m2$gene, deres$mgi_symbol),]
minires$Type <- m0m1m2$group
remgenevp <- ggplot(deres, aes(log2FoldChange, neglogtenp))+
  geom_point(size = 0.1)+
  geom_vline(xintercept = 0, linetype = 'dotted', col = 'red')+
  geom_hline(yintercept = -log10(0.05), linetype = 'dotted', col = 'red')+
  geom_point(data = minires, size = 2, aes(col = Type))+
  ggrepel::geom_label_repel(data = minires,  aes(col = Type, label = mgi_symbol), size = 3, seed = 20223)
remgenevp



deres$neglogtenp <- -log10(deres$pvalue)
minires <- deres[match(m0m1m2$gene, deres$mgi_symbol),]
minires$Type <- m0m1m2$group
genevp <- ggplot(deres, aes(log2FoldChange, neglogtenp))+
  geom_point(size = 0.1)+
  geom_vline(xintercept = 0, linetype = 'dotted', col = 'red')+
  geom_hline(yintercept = -log10(0.05), linetype = 'dotted', col = 'red')+
  geom_point(data = minires, size = 2, aes(col = Type))+
  ggrepel::geom_label_repel(data = minires,  aes(col = Type, label = mgi_symbol), size = 3, seed = 20223)
genevp


# also get mMCP markers
# library(mMCPcounter)
# data("mMCPcounter_signatures_GCRm39")
# saveRDS(mMCPcounter_signatures_GCRm39, "data/MCPCONTERSIGNATURES_JUL25_2023.RDS")
mMCPcounter_signatures_GCRm39 <- readRDS("data/MCPCONTERSIGNATURES_JUL25_2023.RDS")
mac <- mMCPcounter_signatures_GCRm39[mMCPcounter_signatures_GCRm39$Denomination == 'Monocytes / macrophages',]

#just one gene....
mcpdf <- data.frame(gene = mac$Gene.Symbol,
                    group = 'mMCPCounter_Macrophage')

m0m1m2_mcp <- rbind(m0m1m2, mcpdf)





types <- unique(m0m1m2_mcp$group)
sl <- lapply(types, function(mt){
  genes <- m0m1m2_mcp[m0m1m2_mcp$group==mt,'gene']
  
  mat <- tpm[rownames(tpm) %in% genes,]
  
  #first log
  mat <- log1p(mat)
  
  #next get sample-wise mean
  ss <- Matrix::colMeans(mat)
  
  outdf <- as.data.frame(t(data.frame(score = ss, row.names = names(ss))))
  rownames(outdf) <- mt
  
  return(outdf)
  
})


scores = dplyr::bind_rows(sl)


#make barplot
ggs <- reshape2::melt(t(scores))

colnames(ggs) <- c('Code', 'Type', 'Score')

ggplot(ggs, aes(Code, Score, fill = Type))+
  geom_bar(stat="identity")


















#get pways
# pathways <- as.data.frame(msigdbr(species = 'Mus musculus'))
# dir.create('results/comparative-de/TKO-vs-DKO/4.downstream/review/')
# saveRDS(pathways, "results/comparative-de/TKO-vs-DKO/4.downstream/review/msigdbr.rds")
pathways <- readRDS("results/comparative-de/TKO-vs-DKO/4.downstream/review/msigdbr.rds")

#subset just go bp
pathways <- pathways[pathways$gs_subcat == 'GO:BP',]

# pwaynames <- c("GOBP_OSTEOCLAST_DIFFERENTIATION", "GOBP_MACROPHAGE_ACTIVATION")

allpwaynames <- unique(pathways$gs_name)
pwaynames <- allpwaynames[grepl('OSTEOCLAST', allpwaynames)]
pwaynames <- allpwaynames[grepl('OSTEOBLAST', allpwaynames)]
pwaynames <- allpwaynames[grepl('MACROPHAGE', allpwaynames)]
pwaynames <- allpwaynames[grepl('MONOCYT', allpwaynames)]
pwaynames <- pwaynames[!grepl('FACTOR', pwaynames)]
pwaynames <- pwaynames[!grepl('REGULATION', pwaynames)]
pwaynames <- pwaynames[!grepl('PRODUCTION', pwaynames)]


# pwaynames <- c(pwaynames, allpwaynames[grepl('MONOCYT', allpwaynames)] )
# pwaynames <- c(pwaynames, allpwaynames[grepl('MACROPHAGE', allpwaynames)] )


pways_of_interest <- lapply(pwaynames, function(pway){
  
  pwaysub <- pathways[pathways$gs_name == pway,]
  
  return(pwaysub$gene_symbol)
  
})
names(pways_of_interest) <- pwaynames

# rm(pathways)





#calculate module scores

mmlist <- lapply(1:length(pways_of_interest), function(i){
  
  gl <- pways_of_interest[[i]]
  modulescore(tpm, gl)
})


mm <- as.data.frame(dplyr::bind_rows(mmlist))
rownames(mm) <- pwaynames
colnames(mm) <- colnames(tpm)


# mm <- modulescore(tpm, pways_of_interest$GOBP_OSTEOCLAST_DIFFERENTIATION)
mm
mm <- t(scale(t(mm)))


library(ComplexHeatmap)
hm <- ComplexHeatmap::Heatmap(mm, 
                              # row_names_side = 'left',
                              row_title_gp = grid::gpar(fontsize = 5),
                              width = ncol(mm)*unit(20, "mm"),
                              name = 'Scaled\nTPM'
)






pdf("~/Desktop/R3_Osteoblast_module_heatmap.pdf", height = 12, width = 20)
draw(hm, heatmap_legend_side = "left")
dev.off()




### heatmap of genes
pway_genes_to_plot = "GOBP_OSTEOCLAST_DIFFERENTIATION"
pway_genes_to_plot = "GOBP_OSTEOBLAST_DIFFERENTIATION"
pwaygenes <- pways_of_interest[[pway_genes_to_plot]]
pwaygenes <- pwaygenes[pwaygenes %in% rownames(tpm)]
mat <- tpm[rownames(tpm) %in% pwaygenes,]

mat <- t(scale(t(mat)))

hmg <- ComplexHeatmap::Heatmap(mat, 
                               # row_names_side = 'left',
                               row_title_gp = grid::gpar(fontsize = 3),
                               width = ncol(mm)*unit(20, "mm"),
                               name = 'Scaled\nTPM',
                               column_title = pway_genes_to_plot
)
hmg



pdf("~/Desktop/R3_geneheatmap_GOBP_OSTEOBLAST_DIFFERENTIATION.pdf", height = 25, width = 20)
draw(hmg, heatmap_legend_side = "left")
dev.off()

# library(GSVA)
# 
# mmgsva <- GSVA::gsva(expr = as.matrix(tpm),
#                      gset.idx.list = pways_of_interest)
# 
# mmgsea <- GSVA::gsva(expr = as.matrix(tpm),
#                      gset.idx.list = pways_of_interest,
#                      method='ssgsea'
#                      )













### m0 vs m1 vs m2


m1 <- 'Ccl5, Ccr7, Cd40, Cd86, Cxcl9, Cxcl10, Cxcl11, Ido1, Il1a, Il1b, Il6, Irf1, Irf5, Kynu'
m1 <- str_split_fixed(m1, ', ', Inf)[1,]



m2 <- "Ccl4, Ccl13, Ccl18, Ccl20, Ccl22, Cd276, Clec7a, Ctsa, Ctsb, Ctsc, Ctsd, Fn1, Il4r, Irf4, Lyve1, Mmp9, Mmp14, Mmp19, Msr1, Tgfb1, Tgfb2, Tgfb3, Tnfsf8, Tnfsf12, Vegfa, Vegfb, Vegfc"
m2 <- str_split_fixed(m2, ', ', Inf)[1,]



m1m2 <- data.frame(gene = c(m1,m2),
                   group =  c(rep('M1', length(m1)), rep('M2', length(m2)))
)

m1m2 <- m1m2[m1m2$gene %in% rownames(tpm),]

mat <- tpm[match(m1m2$gene, rownames(tpm)),]


mat <- log1p(mat)
# mat <- t(scale(t(mat)))
mat <- scale(mat)

hmm1m2 <- ComplexHeatmap::Heatmap(mat, 
                                  # row_names_side = 'left',
                                  # row_title_gp = grid::gpar(fontsize = 5),
                                  width = ncol(mat)*unit(20, "mm"),
                                  name = 'Within\nSample\nScaled\nTPM',
                                  column_title = "M1 vs M2",
                                  row_split = m1m2$group,
                                  column_split = md$Condition
)
hmm1m2









### try heatmap of m0, m1, m2 markers


markerfile <- '~/Dropbox/data/deyou/macspectrum/finalized_m0m1m2_markers.csv'
m <- read.csv(markerfile)


### top 20 genes ###
m <- m[m$gene %in% rownames(tpm),]
n <- 20
top <- m %>% group_by(cluster) %>% top_n(n = n, wt = avg_log2FC)




mat <- tpm[match(top$gene, rownames(tpm)),]


mat <- log1p(mat)
# mat <- t(scale(t(mat)))
# mat <- scale(mat)


#just tko
columnsplitnames <- md$Condition
mat <- mat[,grepl('TKO', colnames(mat))]
columnsplitnames <- rep('TKO', 3)


hm_ms <- ComplexHeatmap::Heatmap(mat, 
                                 # row_names_side = 'left',
                                 # row_title_gp = grid::gpar(fontsize = 5),
                                 width = ncol(mat)*unit(20, "mm"),
                                 name = 'Log1p\nTPM',
                                 column_title = "MacSpectrum M0, M1, M2 markers",
                                 row_split = top$cluster,
                                 column_split = columnsplitnames
)
hm_ms





## all ms genes
top <- m
mat <- tpm[match(top$gene, rownames(tpm)),]


mat <- log1p(mat)
# mat <- t(scale(t(mat)))
# mat <- scale(mat)


#just tko
columnsplitnames <- md$Condition
mat <- mat[,grepl('TKO', colnames(mat))]
columnsplitnames <- rep('TKO', 3)


hm_ms_all <- ComplexHeatmap::Heatmap(mat, 
                                     # row_names_side = 'left',
                                     # row_title_gp = grid::gpar(fontsize = 5),
                                     width = ncol(mat)*unit(20, "mm"),
                                     name = 'Log1p\nTPM',
                                     column_title = "MacSpectrum M0, M1, M2 markers",
                                     row_split = top$cluster,
                                     show_row_names = F,
                                     column_split = columnsplitnames
)
hm_ms_all


clusters <- unique(m$cluster)
mmlist <- lapply(clusters, function(clust){
  modulescore(tpm, m[m$cluster==clust,'gene'])
})


mm <- as.data.frame(dplyr::bind_rows(mmlist))
rownames(mm) <- clusters
colnames(mm) <- colnames(tpm)



mm <- mm[,grepl('TKO', colnames(mm))]
#mm <- t(scale(t(mm)))
mm <- scale(mm)

hm_ms_module <- hm <- ComplexHeatmap::Heatmap(mm, 
                                              # row_names_side = 'left',
                                              row_title_gp = grid::gpar(fontsize = 5),
                                              width = ncol(mm)*unit(20, "mm"),
                                              name = 'Scaled\nTPM'
)
hm_ms_module
