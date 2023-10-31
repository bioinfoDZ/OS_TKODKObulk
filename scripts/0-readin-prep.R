library(tidyverse)


### we will use the star reference geneInfo.tab file
gencode <- read.table('data/geneInfo.tab', sep = '\t', skip = 1, header = F)
colnames(gencode) <- c('gene_id', 'gene_name', 'gene_type')


#readin data
md <- as.data.frame(readxl::read_excel('data/metadata.xlsx'))
files <- list.files('data/counts/', recursive = T, full.names = T, pattern = '.tab')

samplist <- list()
for(file in files){
  
  sampname <- str_split_fixed(file, '/', Inf)[,4]
  # basename <- str_sub(sampname, 10)
  basename <- sampname
  
  message(sampname)
  
  # take second column;
  # just pick col with highest colsums
  samp <- read.table(file, sep = '\t', skip = 4)
  samp <- data.frame(row.names = samp[,1],
                     counts = samp[,2])
  
  colnames(samp)[1] <- basename
  
  
  
  samplist[[sampname]] <- samp
  
  rm(samp, sampname, basename)
  
  
}

gem <- dplyr::bind_cols(samplist)
rm(samplist)


#match with metadata
gem <- gem[,match(md$Sample, colnames(gem))]


### check sex ###

#normalize
sizefactors <- DESeq2::estimateSizeFactorsForMatrix(counts = gem)
gemnorm <- as.data.frame(t( t(gem) / sizefactors ))




#read in md, with sex labels
xist <- t(gemnorm['ENSG00000229807',])

df <- data.frame(xist = xist,
                 sample = rownames(xist))

xistplot <- ggplot(df, aes(sample, xist))+
  geom_point()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  ylab('Normalized Xist expression')

xistplot

ggsave(xistplot, filename = 'results/allsamples/xist-expression.jpg', height = 5, width = 5, dpi = 300)


rm(gemnorm, df, xist, xistplot)






#### plot the lib size of all samples #####
pdf <- data.frame(samp = colnames(gem), 
                  numreadsaligned = colSums(gem),
                  condition = md$Treatment,
                  cellline = md$CellLine,
                  color=md$TreatmentColor,
                  stringsAsFactors = F
)


#order them from hi to low
pdf$samp <- factor(pdf$samp, levels = pdf[order(pdf$numreadsaligned, decreasing = T),"samp"])
pdf$condition <- factor(pdf$condition, levels = unique(pdf[order(pdf$numreadsaligned, decreasing = T),"condition"]))
cols <- unique(pdf[order(pdf$numreadsaligned, decreasing = T),"color"])

data.frame(cols = cols, condition = levels(pdf$condition))

libsize_rawall <- ggplot(pdf, aes(x = samp, y = numreadsaligned, fill = condition))+
  geom_bar(stat = 'identity')+
  scale_fill_brewer(palette = 'Set1',direction = -1)+
  scale_y_continuous(labels = scales::comma)+
  theme_light()+
  #scale_fill_brewer(palette = 'Set2')+
  scale_fill_manual(values = cols)+
  #scale_y_continuous(limits = c(0,25000000), labels = scales::comma)+
  coord_flip()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  labs(title = 'Number of aligned reads', 
       subtitle = 'Non-normalized',
       y = 'Number of reads aligned', x = 'Sample')

libsize_rawall

ggsave(libsize_rawall, filename = 'results/allsamples/libsize-whollelib.jpg', height = 5, width = 5, dpi = 300)




### plot lib sizes by aligned type ###

rnatypes <- unique(gencode$gene_type)
typelist <- list()
for(type in rnatypes){
  tgc <- gencode[gencode$gene_type == type,]
  typegem <- gem[rownames(gem) %in% tgc$gene_id,]
  typerow <- data.frame(samp = names(typegem),
                        SpeciesCount = colSums(typegem),
                        type = type)
  typelist[[type]] <- typerow
}
typedf <- dplyr::bind_rows(typelist)

#sort by how many and remove zeros
agg <- aggregate(SpeciesCount ~ type, typedf,sum)
agg <- agg[agg$SpeciesCount > 0,]
agg <- agg[order(agg$SpeciesCount),]

typedf <- typedf[typedf$type %in% agg$type,]
typedf$type <- factor(typedf$type, levels = agg$type)


#set custom color palette
palette <- c(RColorBrewer::brewer.pal(Inf, 'Dark2'),
             RColorBrewer::brewer.pal(Inf, 'Set2'),
             RColorBrewer::brewer.pal(Inf, 'Set1'),
             RColorBrewer::brewer.pal(Inf, 'Paired'),
             RColorBrewer::brewer.pal(Inf, 'Accent')
             
)


set.seed(2021)
pal <- sample(palette, nrow(agg))

#order low to high
typedf$samp <- factor(typedf$samp, levels = levels(pdf$samp))

libsize_species <- ggplot(typedf, aes(fill = type, y = SpeciesCount, x  = samp))+
  geom_bar(position="stack", stat="identity")+
  # scale_y_continuous(limits = c(0,25000000), labels = scales::comma)+
  coord_flip()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  labs(title = 'Number of aligned reads', 
       subtitle = 'Non-normalized library size with RNA species',
       y = 'Number of reads aligned', x = 'Sample')+
  scale_fill_manual(values = pal)

libsize_species

ggsave(libsize_species, filename = 'results/allsamples/libsize-species.jpg', height = 5, width = 10, dpi = 300)


# save and go

saveRDS(file = 'data/gem.rds', gem)

rm(list=ls())




