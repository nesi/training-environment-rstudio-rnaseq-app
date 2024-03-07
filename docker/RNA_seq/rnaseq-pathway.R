## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
knitr::opts_chunk$set(warning = FALSE)
knitr::opts_chunk$set(message = FALSE)


## -----------------------------------------------------------------------------
library(dplyr)
library(goseq)
supportedOrganisms() %>% head()


## ---- eval=FALSE--------------------------------------------------------------
## supportedOrganisms() %>% View()


## -----------------------------------------------------------------------------
load('topTable.RData')


## -----------------------------------------------------------------------------
# Note: If you want to use tt from DESeq, replace $adj.P.Val with $padj below

genes <- ifelse(tt$adj.P.Val < 0.05, 1, 0) 
names(genes) <- rownames(tt)
head(genes)
table(genes)


## -----------------------------------------------------------------------------
pwf=nullp(genes, "sacCer1", "ensGene")


## ---- eval=FALSE, echo=FALSE--------------------------------------------------
## library(GenomicFeatures)
## txdb <- makeTxDbFromGFF(file="Saccharomyces_cerevisiae.R64-1-1.99.gtf", format = "gtf", dataSource = "SGD", organism = "Saccharomyces cerevisiae")
## txsByGene=transcriptsBy(txdb,"gene")
## lengthData=median(width(txsByGene))


## -----------------------------------------------------------------------------
head(pwf)


## -----------------------------------------------------------------------------
par(mfrow=c(1,2))
hist(pwf$bias.data,30)
hist(pwf$pwf,30)


## -----------------------------------------------------------------------------
library(ggplot2)
data.frame(logGeneLength = log2(pwf$bias.data), 
           avgExpr = tt$AveExpr) %>% 
  ggplot(., aes(x=logGeneLength, y=avgExpr)) + 
  geom_point(size=0.2) + 
  geom_smooth(method='lm')


## -----------------------------------------------------------------------------
data.frame(logGeneLength = log2(pwf$bias.data), 
           negLogAdjP = -log10(tt$adj.P.Val)) %>% 
  ggplot(., aes(x=logGeneLength, negLogAdjP)) + 
  geom_point(size=0.2) + 
  geom_smooth(method='lm')


## -----------------------------------------------------------------------------
GO.wall=goseq(pwf, "sacCer1", "ensGene")


## -----------------------------------------------------------------------------
head(GO.wall)


## -----------------------------------------------------------------------------
GO.wall.padj <- p.adjust(GO.wall$over_represented_pvalue, method="fdr")
sum(GO.wall.padj < 0.05)
GO.wall.sig <- GO.wall$category[GO.wall.padj < 0.05]
length(GO.wall.sig)
head(GO.wall.sig)


## -----------------------------------------------------------------------------
GO.wall[GO.wall.padj < 0.05, ] %>% filter(numInCat < 500)


## -----------------------------------------------------------------------------
library(GO.db)

GOTERM[[GO.wall.sig[1]]]


## -----------------------------------------------------------------------------
GO.nobias=goseq(pwf, "sacCer1", "ensGene", method="Hypergeometric")


## -----------------------------------------------------------------------------
head(GO.nobias)


## -----------------------------------------------------------------------------
GO.nobias.padj <- p.adjust(GO.nobias$over_represented_pvalue, method="fdr")
sum(GO.nobias.padj < 0.05)
GO.nobias.sig <- GO.nobias$category[GO.nobias.padj < 0.05]
length(GO.nobias.sig)
head(GO.nobias.sig)


## -----------------------------------------------------------------------------
library(gplots)
venn(list(GO.wall=GO.wall.sig, GO.nobias=GO.nobias.sig))


## -----------------------------------------------------------------------------
## Only significant in Hypergeomtric analysis
onlySig.nobias <- setdiff(GO.nobias.sig, GO.wall.sig)

## Only significant in Wallenius analysis
onlySig.wall <- setdiff(GO.wall.sig, GO.nobias.sig)

## Significant in both
sig.wall.nobias <- intersect(GO.wall.sig, GO.nobias.sig)


## -----------------------------------------------------------------------------
len=getlength(names(genes),"sacCer1","ensGene")

head(len)

go = getgo(names(genes),"sacCer1","ensGene")

names(go) %>% head()
class(go)

head(go[[1]])
head(go[[2]])


## -----------------------------------------------------------------------------
lengths.onlySig.nobias <- list()

for(i in 1:length(onlySig.nobias)){
  inGo <- lapply(go, function(x)  onlySig.nobias[i] %in% x) %>% unlist()
  lengths.onlySig.nobias[[i]] <- len[inGo]
}

lengths.onlySig.wall <- list()

for(i in 1:length(onlySig.wall)){
  inGo <- lapply(go, function(x)  onlySig.wall[i] %in% x) %>% unlist()
  lengths.onlySig.wall[[i]] <- len[inGo]
}


## -----------------------------------------------------------------------------
cols <- rep(c("lightpink", "lightblue"), c(10,7))
boxplot(c(lengths.onlySig.nobias, lengths.onlySig.wall), col=cols)


## -----------------------------------------------------------------------------
lengths.sig.wall.nobias <- list()

for(i in 1:length(sig.wall.nobias)){
  inGo <- lapply(go, function(x)  sig.wall.nobias[i] %in% x) %>% unlist()
  lengths.sig.wall.nobias[[i]] <- len[inGo]
}


## -----------------------------------------------------------------------------
cols <- rep(c("lightpink", grey(0.7), "lightblue"), c(10,37,7))

avgLength <- lapply(c(lengths.onlySig.nobias, lengths.sig.wall.nobias, lengths.onlySig.wall),
                    median) %>% unlist()

oo <- order(avgLength, decreasing=TRUE)


## -----------------------------------------------------------------------------
boxplot(c(lengths.onlySig.nobias, lengths.sig.wall.nobias, lengths.onlySig.wall)[oo],
        col=cols[oo], ylab="Gene Length", xlab = "GO term")


## -----------------------------------------------------------------------------
avgLength.wall <- lapply(c(lengths.onlySig.wall, lengths.sig.wall.nobias), median)

avgLength.nobias <- lapply(c(lengths.onlySig.nobias, lengths.sig.wall.nobias), median)

cols <- rep(c("blue", "lightblue", "red","lightpink"),
            c(length(lengths.onlySig.wall), length(lengths.sig.wall.nobias),
              length(lengths.onlySig.nobias), length(lengths.sig.wall.nobias)))

plot(c(avgLength.wall, avgLength.nobias), 
     -log(c(GO.nobias.padj[GO.nobias.padj < 0.05], GO.wall.padj[GO.wall.padj < 0.05])), 
     col=cols, pch=16, xlab="Median Gene Length", ylab ="-log(FDR adj-pval)")
legend('topright', c("Only sig in NoBias", "Sig in both (nobias adjp)", 
                     "Sig in both (wal adjp)", "Only sig in Wall"), 
       fill=c("red", "pink", "lightblue", "blue"))

