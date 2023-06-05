---
title: "Notebook4_E18.5_and_preterm_P0"
author: "Laura Wolbeck"
date: "2023-02-28"
---

In this notebook DEG and GSEA analysis is performed for samples E 18.5 and preterm P0. Moreover, this script includes code for all figures of E 18.5 and preterm P0 samples. 
For DEG and GSEA the package "cacoa" was used: https://github.com/kharchenkolab/cacoa/tree/main.


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
```

# 1. DEG analysis
Read in high resolution annotation file "anno"
Read in conos object from notebook1 as "con"

Create condition factor for cacoa object
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
write.xlsx(x = DEG$RG, file= "DEG_early_RG.xlsx")
```

# 2. GSEA
```{r}
cao$estimateOntology(type="GSEA", org.db=org.Mm.eg.db)
```

save GSEA results of RG as excel file
```{r}
cao$saveOntologyAsTable("early_GSEA.tsv", name="GSEA")
GSEA <- read.table("early_GSEA.tsv", sep="\t")
write.xlsx(x = GSEA[GSEA$V1=="RG",], file= "GSEA_early_RG.xlsx")
```

# 3. Figures
## Fig S5 D
```{r}
#defining colours of clusters
colours_high <- c("#e1af32","#FFAF32","#3232FF","#32FF32","#FFFF00","#CDCD32","#FF9b9b","#EB8C32","#C8961E","#FF3232","#AF3232","#B49191","#694646","#820000","#965050","#550505","#C80F0F","#AF32FF","#E632E6","#FFC8FF","#DBFF00", "#0096FF") %>% setNames(c("Immature ependymal cells-2","Immature ependymal cells-2","RG","Intermediate progenitor cells","Dividing cells-2","Dividing cells-1","Neuroblasts-1", "Dividing cells-4", "Dividing cells-3","Neuroblasts-2","Neurons 1-1","Neurons 1-2","Neurons 2","Neurons 3-3","Neurons 3-4","Neurons 3-1","Neurons 3-2","OPC","Microglia", "Endothelial/Pericytes/VSMC","Erythrocyte", "MKi67+_RG"))

con$plotGraph(groups=anno, label="") + scale_colour_manual(values=colours_high)+  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
```

## Fig. S5 E
```{r}
#define words to exclude from collapsed GO terms
ex_words <- c('regulation', 'process', 'cell')
RG <- cao$plotOntologyHeatmapCollapsed(
  name="GSEA", genes="down", n=70, size.range=c(1, 4), subtype="BP", exclude.words=ex_words, clust.method="ward.D")

#filter for Radial glia
RG$data %<>% filter(G2=="RG", value > 0)

df <- RG$data
df %<>% arrange(value) %>% mutate(G1 = G1 %>% factor(., levels = .))

#define color legend
color.guide <- guide_colorbar(title = "-log10(p-value)", title.position = "left", 
        title.theme = element_text(angle = 90, hjust = 0.5))


ggplot(df, aes(x=G2, y = G1, fill = value)) + geom_tile( colour="white", size=0.2) +  scale_fill_distiller(direction = 1 , palette = "Blues") + guides(fill = color.guide) + labs(title="Top 30 down of RG ", y = " ", x= " ") + theme_grey(base_size = 11) + theme(axis.ticks.x = element_blank(),  axis.text.x = element_blank())
```

## Fig. S5 F
```{r, fig.width=7, fig.height=7}
#define words to exclude from collapsed GO terms
ex_words <- c('regulation', 'process', 'cell')
RG <- cao$plotOntologyHeatmapCollapsed(
  name="GSEA", genes="up", n=100, size.range=c(1, 4), subtype="BP", exclude.words=ex_words, clust.method="ward.D")

RG$data %<>% filter(G2=="RG", value > 0)

#filter for Radial glia
df <- RG$data
df %<>% arrange(value) %>% mutate(G1 = G1 %>% factor(., levels = .))

#define color legend
color.guide <- guide_colorbar(title = "-log10(p-value)", title.position = "left", 
        title.theme = element_text(angle = 90, hjust = 0.5))

ggplot(df, aes(x=G2, y = G1, fill = value)) + geom_tile( colour="white", size=0.2) +  scale_fill_distiller(direction = 1 , palette = "Reds") + guides(fill = color.guide) + labs(title="Top 30 up of RG", y = " ", x= " ") + theme_grey(base_size = 11) + theme(axis.ticks.x = element_blank(),  axis.text.x = element_blank())
```
