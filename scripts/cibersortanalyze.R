library(tidyverse)



cs <- read.csv('data/cibersort/CIBERSORTx_Job14_Results.csv')


rownames(cs) <- cs$Mixture

cs <- cs[,c(-1, -13, -14, -15, -16)]

md <- readxl::read_excel('data/metadata.xlsx')
cs <- cs[match(md$Sample, rownames(cs)),]


md <- cbind(md, cs)


cts <- colnames(cs)

plotlist <- list()

for(ct in cts){
  
  pdf <- data.frame(samp  = md$Sample,
                    cond = md$Condition,
                    score = md[,ct])
  
  plotlist[[ct]] <- ggplot(pdf, aes(x = samp, y = score, col = cond))+
    geom_point()+ 
    scale_color_manual(values =md$Color)+
    ggtitle(ct)+
    Seurat::NoLegend()
}

patchwork::wrap_plots(plotlist)
