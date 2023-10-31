library(tidyverse)

#format:
# first row with # is the sample names
# next 20 rows are clinical information and sample metadata
# rest of data is the actual gene expression


### DOWNLOADED 2022.04.04
# APRIL 4 2022

raw <- read.table('data/kuijjer/ps_avgpres_osteosarcoma127_ilmnhwg6v2_box1649104359-datagrabber.txt',
           sep = '\t',  comment.char = '', na.strings = 'nd')

#take a look
raw[1:30,1:6]


#separate the gem vs the metadata
md <- raw[1:21,]

# get the gene exp matrix
gem <- raw[22:nrow(raw),]

#remove raw from mem
rm(raw)

###fix it up ###
#remove hash (#) from col names
md$V1 <- gsub('#', '', md$V1)

#remove colon form colnames, just the first one
md$V1[1] <- 'hugo'


# grab sample names
colnames(gem) <- md[1,]




### separate out gene names
genes <- gem[,1:2]

#fix matrix
gem <-gem[,3:ncol(gem)]

#make matrix numeric
gem <- as.data.frame(apply(gem, 2, as.numeric))

#make rownames of gem to the illumina probes
rownames(gem) <- genes$probeset



## fix up the md ##

#remove extraneous first column
md <- md[,-1]

#transpose the md
md <- t(md)

#get appropriate colnames, then remove extraneous column
colnames(md) <- md[1,]
md <- md[-1,]
md <- as.data.frame(md)
rownames(md) <- NULL



#make appropriate columns numeric
# to check:
apply(md, 2, head)

# fix numeric cols
md$age <- as.numeric(md$age)
md$tod <- as.numeric(md$tod)
md$tof <- as.numeric(md$tof)
md$tom1 <- as.numeric(md$tom1)
md$tom2 <- as.numeric(md$tom2)


#finally, rename "probest" column to "sample"
colnames(md)[1] <- 'sample'



# at this point, it's all good to go i think




### test how to deal with the duplicate genes...

#default
x <- colSums(gem)


#with non-dup genes...
dups <- names( table(genes$hugo)[table(genes$hugo)>1] )
nondupgenes <- genes[!(genes$hugo %in% dups),]

gem_nondups <- gem[match(nondupgenes$probeset, rownames(gem)),]

y <- colSums(gem_nondups)


xx <- data.frame(x,y)

plot(x,y)


### it seems they make up a decent amount of the library...

rm(x,y,dups,nondupgenes,gem_nondups, xx)


### check the clinical stuff:
# survival status, time vars...

# get only tumor
df <- md[md$type == 'biopsy',]



dead <- df[df$deceased=='true',c('deceased', 'tod', 'tof')]
alive <- df[df$deceased=='false',c('deceased', 'tod', 'tof')]

#it looks like tof = time of followup
#tod = time of death
# for dead patients, tod = tof
# for alive patients tof is there, tod is missing

#tom1 and tom2 -= time of mets? multiple mts?

#the r2 website indicates these variables are in MONTHS

#the age info seems to be in months also; 
# dividing by 12 to get years does seem to match expected age...




rm(alive,dead,df)










###save:
# raw gem
# raw md...

savegem <- gem
savemd <- md
savegenes <- genes



write.csv(savegem, file = 'data/kuijjer/parsed-gem.csv', quote = F, row.names = T)
write.csv(savemd, file = 'data/kuijjer/parsed-md.csv', quote = F, row.names = F)
write.csv(savegenes, file = 'data/kuijjer/parsed-genesymbols.csv', quote = F, row.names = F)










### run pca ###

md <- md[md$type %in% c('biopsy', 'resection'),]
gem <- gem[,match(md$sample, colnames(gem))]

#PCA

numhvgs <- 2000

#transform; sqrt or log. can also try VST or rlog...
gemt <- sqrt(gem)


#get HVGs
rowvars <- apply(gemt, 1, var)

#select top hvgs
rowvars <- sort(rowvars, decreasing = T)[1:numhvgs]

#get a small mat
gem_hvgs <- gemt[rownames(gemt) %in% names(rowvars),]

#scale small mat
gem_scale <- t(scale(t(gem_hvgs)))




#pca
pca <- prcomp(t(gem_scale))

rm(gemt, gem_hvgs, rowvars, gem_scale)

plot(pca$sdev)

#get embeddings
emb <- as.data.frame(pca$x)

#get genes
pcagenes <- pca$rotation
head(sort(pcagenes[,1],decreasing = T))
head(sort(pcagenes[,2],decreasing = T))


emb$libsize <- colSums(gem)

#plot PCs
libsize <- ggplot(emb, aes(PC1, PC2, col=libsize))+
  geom_point()+
  scale_color_distiller(palette = 'Reds', direction = 1)

#plot lab
emb <- emb <- cbind(emb,md)

lab <- ggplot(emb, aes(PC1, PC2, col=lab))+
  geom_point()



ggplot(emb, aes(PC1, PC2, col=type))+
  geom_point()


lab / libsize

  