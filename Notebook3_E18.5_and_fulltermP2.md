---
title: "Notebook3_E18.5_and_fullterm_P2"
author: "Laura Wolbeck"
date: "2023-02-23"
---

In this notebook UMAP embedding, DEG and GSEA analysis will be performed for full-term samples E18.5 and P2. Moreover, this script includes code for all figures of full-term samples. 
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

# 1. Read in count matrices (cms)
Build cms from filtered cms from notebook 1

```{r, eval=F, echo=F}
#load filtered cms
#extract fullterm samples
cms <- append(cms1[c(1, 2)], cms2[c(1,2)])
```

```{r}
names = c("E18_5_1_",
           "E18_5_2_",
           "full_P2_1_",
           "full_P2_2_")
```

# 2. UMAP embedding and clustering
pre-processing of each dataset with pagoda2, followed by using conos to build the joint UMAP graph with forced alignment (alignment.strength = 0.2) to integrate samples well
```{r, eval=FALSE}
con <- quickConos(cms,
                  names,
                  n.cores.p2=10,
                  n.cores.con=20, get.tsne = TRUE, alignment.strength=0.2)
con <- con$con
```

## Additional clusters
Rerun leiden clustering to find additional clusters
```{r, eval=FALSE}
con$findCommunities(method = leiden.community, resolution = 5, min.group.size = 15)
leiden30_full <- con$clusters$leiden$groups %>% factor
```
The final clusters were generated using a combination of manual selection and findSubcommunities() function of conos to split the clusters further. The clusters were annotated using known cell-type specific marker genes from the literature.

# 3. DEG analysis
Read in high resolution annotation file "anno"

Create condition factor for cacoa object
```{r,eval=FALSE }
condition <- setNames(c("early", "early", "late", "late"), names(con$samples))
```

Create Cacoa object
```{r,eval=FALSE}
cao <- Cacoa$new(con, cell.groups = anno, sample.groups=condition, n.cores = 10, ref.level = "early", target.level = "late")
```

estimate DE genes
```{r,eval=FALSE}
cao$estimatePerCellTypeDE(n.cores=50) 
```

save DEG of RG as excel file 
```{r}
DEG <- cao$test.results$de %>% lapply("[[", "res") 
write.xlsx(x = DEG$RG, file= "DEG_fullterm_RG.xlsx")
```

# 4. GSEA
```{r}
cao$estimateOntology(type="GSEA", org.db=org.Mm.eg.db)
```

save GSEA results of RG as excel file
```{r}
cao$saveOntologyAsTable("fullterm_GSEA.tsv", name="GSEA")
GSEA <- read.table("fullterm_GSEA.tsv", sep="\t")
write.xlsx(x = GSEA[GSEA$V1=="RG",], file= "GSEA_fullterm_RG.xlsx")
```

# 5. Figures
## Fig 1H upper plot
```{r}
enrichplot::dotplot(cao$test.results$GSEA$res$RG$BP,showCategory=c("ribosome biogenesis" , "ribonucleoprotein complex biogenesis", "ribonucleoprotein complex subunit organization", "ribonucleoprotein complex assembly"))
```
## Fig S5A
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

## Fig S7A
```{r}
#defining colours of high resolution clusters
colours_high <- c("#e1af32","#FFAF32","#3232FF","#32FF32","#FFFF00","#CDCD32","#FF9b9b","#EB8C32","#C8961E","#FF3232","#AF3232","#B49191","#694646","#820000","#965050","#550505","#C80F0F","#AF32FF","#E632E6","#FFC8FF","#DBFF00") %>% setNames(c("Immature ependymal cells-2","Immature ependymal cells-2","RG","Intermediate progenitor cells","Dividing cells-2","Dividing cells-1","Neuroblasts-1", "Dividing cells-4", "Dividing cells-3","Neuroblasts-2","Neurons 1-1","Neurons 1-2","Neurons 2","Neurons 3-3","Neurons 3-4","Neurons 3-1","Neurons 3-2","OPC","Microglia", "Endothelial/Pericytes/VSMC","Erythrocytes"))
```

```{r}
con$plotGraph(groups=anno, embedding="90", label="") + scale_colour_manual(values=colours_high) +  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
```

## Fig. S7B
```{r}
#define words to exclude from colapsed GO terms
ex_words <- c('regulation', 'process', 'cell')
RG <- cao$plotOntologyHeatmapCollapsed(
  name="GSEA", genes="down", n=50, size.range=c(1, 4), subtype="BP", exclude.words=ex_words, clust.method="ward.D")

#filter for Radial glia
RG$data %<>% filter(G2=="RG", value > 0)

df <- RG$data
df %<>% arrange(value) %>% mutate(G1 = G1 %>% factor(., levels = .))

#define color legend
color.guide <- guide_colorbar(title = "-log10(p-value)", title.position = "left", 
        title.theme = element_text(angle = 90, hjust = 0.5))


ggplot(df, aes(x=G2, y = G1, fill = value)) + geom_tile( colour="white", size=0.2) +  scale_fill_distiller(direction = 1 , palette = "Blues") + guides(fill = color.guide) + labs(title="Top 30 down of RG fullterm", y = " ", x= " ") + theme_grey(base_size = 11) + theme(axis.ticks.x = element_blank(),  axis.text.x = element_blank())
```

## Fig. S7C
```{r, fig.width=7, fig.height=7}
#define words to exclude from colapsed GO terms
ex_words <- c('regulation', 'process', 'cell')
RG <- cao$plotOntologyHeatmapCollapsed(
  name="GSEA", genes="up", n=46, size.range=c(1, 4), subtype="BP", exclude.words=ex_words, clust.method="ward.D")

#filter for Radial glia
RG$data %<>% filter(G2=="RG", value > 0)

df <- RG$data
df %<>% arrange(value) %>% mutate(G1 = G1 %>% factor(., levels = .))

#define color legend
color.guide <- guide_colorbar(title = "-log10(p-value)", title.position = "left", 
        title.theme = element_text(angle = 90, hjust = 0.5))

ggplot(df, aes(x=G2, y = G1, fill = value)) + geom_tile( colour="white", size=0.2) +  scale_fill_distiller(direction = 1 , palette = "Reds") + guides(fill = color.guide) + labs(title="Top 30 up of RG fullterm", y = " ", x= " ") + theme_grey(base_size = 11) + theme(axis.ticks.x = element_blank(),  axis.text.x = element_blank())
```
