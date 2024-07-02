---
title: "Notebook5_fulltermP2_and_pretermP3"
author: "Laura Wolbeck"
date: "2023-02-28"
---

In this notebook DEG and GSEA analysis is performed for samples fullterm P2 and preterm P3. Moreover, this script includes code for all figures of fullterm P2 and preterm P3 samples. 
For DEG and GSEA the package "cacoa" was used: https://github.com/kharchenkolab/cacoa/tree/main

# Setup
```{r setup, results='hide', warning=FALSE, message = FALSE}
#Load helper functions
source("scRNA_helper_functions_preterm.R")

#Load libraries
library(conos) #version 1.4.0
library(pagoda2)
library(magrittr)
library(ggplot2)
library(cacoa)
library(xlsx)
library(enrichplot)
```

# 1. DEG analysis
Read in high resolution annotation file "anno"
Con2 from notebook1 is read in as "con"

Create condition factor
```{r,eval=FALSE }
condition <- setNames(c("fullterm", "fullterm", "preterm", "preterm"), names(con$samples))
```

Create Cacoa object
```{r,eval=FALSE}
cao <- Cacoa$new(con, cell.groups = anno, sample.groups=condition, n.cores = 10, ref.level = "fullterm", target.level = "preterm")
```

estimate DE genes
```{r,eval=FALSE}
cao$estimatePerCellTypeDE(n.cores=50)  #DEs not calculated for Microglia
```

save DEG of RG as excel file 
```{r}
DEG <- cao$test.results$de %>% lapply("[[", "res") 
write.xlsx(x = DEG$RG, file= "DEG_late_RG.xlsx")
```

# 2. GSEA
```{r}
cao$estimateOntology(type="GSEA", org.db=org.Mm.eg.db)
```

save GSEA results of RG as excel file
```{r}
cao$saveOntologyAsTable("late_GSEA.tsv", name="GSEA")
GSEA <- read.table("late_GSEA.tsv", sep="\t")
write.xlsx(x = GSEA[GSEA$V1=="RG",], file= "GSEA_late_RG.xlsx")
```

# 3. Figures

## Fig 1H lower plot
```{r}
enrichplot::dotplot(cao$test.results$GSEA$res$RG$BP,showCategory=c("ribonucleoprotein complex biogenesis","ribonucleoprotein complex assembly","ribonucleoprotein complex subunit organization"))
```
## Fig S5C
```{r}
cao$plotEmbedding(color.by='cell.groups', alpha=0.1, size=0.2, title='', 
                    plot.na=FALSE, show.legend=F, font.size=c(2,3))
```
```{r}
cao$estimateCellLoadings()
```
```{r, fig.height=6}
cao$plotCellLoadings(signif.threshold=0.05, show.pvals = F) + xlab("separation coefficient") 
ggsave("Compositionshifts_fullterm.pdf",width=6, height=4.5)
```

## Fig S7G
```{r}
#defining colours of clusters
colours_high <- c("#e1af32","#FFAF32","#3232FF","#32FF32","#FFFF00","#CDCD32","#FF9b9b","#EB8C32","#C8961E","#FF3232","#AF3232","#B49191","#694646","#820000","#965050","#550505","#C80F0F","#AF32FF","#E632E6","#FFC8FF","#DBFF00") %>% setNames(c("Immature ependymal cells-2","Immature ependymal cells-2","RG","Intermediate progenitor cells","Dividing cells-2","Dividing cells-1","Neuroblasts-1", "Dividing cells-4", "Dividing cells-3","Neuroblasts-2","Neurons 1-1","Neurons 1-2","Neurons 2","Neurons 3-3","Neurons 3-4","Neurons 3-1","Neurons 3-2","OPC","Microglia", "Endothelial/Pericytes/VSMC","Erythrocyte"))

con$plotGraph(groups=anno, label="") + scale_colour_manual(values=colours_high)+  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
```

## Fig S7H

```{r}
#define words to exclude from collapsed GO terms
ex_words <- c('regulation', 'process', 'cell', 'a', 'an')
RG <- cao$plotOntologyHeatmapCollapsed(
  name="GSEA", genes="down", n=300, size.range=c(1, 4), subtype="BP", exclude.words=ex_words, clust.method="ward.D")

#filter for Radial glia
RG$data %<>% filter(G2=="RG", value > 0)

df <- RG$data
df %<>% arrange(value) %>% mutate(G1 = G1 %>% factor(., levels = .))

#define color legend
color.guide <- guide_colorbar(title = "-log10(p-value)", title.position = "left", 
        title.theme = element_text(angle = 90, hjust = 0.5))

ggplot(df, aes(x=G2, y = G1, fill = value)) + geom_tile( colour="white", size=0.2) +  scale_fill_distiller(direction = 1 , palette = "Blues") + guides(fill = color.guide) + labs(title="Top 20 down of RG ", y = " ", x= " ") + theme_grey(base_size = 11) + theme(axis.ticks.x = element_blank(),  axis.text.x = element_blank())
```

## Fig S7I
```{r, fig.width=7, fig.height=7}
#define words to exclude from collapsed GO terms
ex_words <- c('regulation', 'process', 'cell')
RG <- cao$plotOntologyHeatmapCollapsed(
  name="GSEA", genes="up", n=118, size.range=c(1, 4), subtype="BP", exclude.words=ex_words, clust.method="ward.D")

RG$data %<>% filter(G2=="RG", value > 0)

#filter for Radial glia
df <- RG$data
df %<>% arrange(value) %>% mutate(G1 = G1 %>% factor(., levels = .))

#define color legend
color.guide <- guide_colorbar(title = "-log10(p-value)", title.position = "left", 
        title.theme = element_text(angle = 90, hjust = 0.5))

ggplot(df, aes(x=G2, y = G1, fill = value)) + geom_tile( colour="white", size=0.2) +  scale_fill_distiller(direction = 1 , palette = "Reds") + guides(fill = color.guide) + labs(title="Top 30 up of RG", y = " ", x= " ") + theme_grey(base_size = 11) + theme(axis.ticks.x = element_blank(),  axis.text.x = element_blank())
```
# Fig S8A
