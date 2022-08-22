getwd()
library(DESeq2)
library(parathyroidSE)
data("parathyroidGenesSE")

head(assay(parathyroidGenesSE))
colData(parathyroidGenesSE)
rowRanges(parathyroidGenesSE)

as.data.frame(colData(parathyroidGenesSE)[, c("sample","patient","treatment", "time")])


allColSamples <- colData(parathyroidGenesSE)$sample
sp <- split(seq(along = allColSamples), allColSamples)

countdata <- sapply(sp, function(columns) {
  rowSums(assay(parathyroidGenesSE)[, columns, drop = FALSE])
})
head(countdata)

coldata <- colData(parathyroidGenesSE)[sapply(sp, `[`, 1), ]
rownames(coldata) <- coldata$sample
coldata

coldata <- coldata[, c("patient", "treatment", "time")]
head(coldata)

rowdata <-  rowRanges(parathyroidGenesSE)
rowdata

ddsFull <- DESeqDataSetFromMatrix(countData = countdata, colData = coldata, 
                                  design = ~ patient + treatment, rowData = rowdata)
ddsFull

dds <- ddsFull[, colData(ddsFull)$treatment %in% c("Control", "DPN") & colData(ddsFull)$time == "48h"]

dds$patient <- factor(dds$patient)
dds$treatment <- factor(dds$treatment)

dds$treatment <- relevel(dds$treatment, "Control")
colData(dds)

dds <- DESeq(dds)
res <- results(dds)
res
mcols(res)

sum(res$pvalue < 0.05, na.rm = TRUE)

sum(res$padj < 0.1, na.rm = TRUE)

resSig <- res[which(res$padj < 0.1), ]
head(resSig[order(resSig$log2FoldChange), ])

tail(resSig[order(resSig$log2FoldChange), ])

plotMA(dds, ylim = c(-2.5,2.5))

plotDispEsts(dds)

filterThreshold <- 2
keep <- rowMeans(counts(dds, normalized = TRUE)) > filterThreshold
table(keep)       

table(p.adjust(res$pvalue, method = "BH") < 0.1)
table(p.adjust(res$pvalue[keep], method = "BH") < 0.1)
hist(res$pvalue[keep], breaks = 100)

library(readr)
library("org.Hs.eg.db")
columns(org.Hs.eg.db)

convertIDs <- function(ids, fromKey, toKey, db, ifMultiple = c("putNA", "useFirst")) { 
  stopifnot(inherits(db, "AnnotationDb"))
  ifMultiple <- match.arg(ifMultiple)
  suppressWarnings(selRes <- AnnotationDbi::select(db, keys = ids, keytype = fromKey, columns = c(fromKey, toKey))) 
  if (ifMultiple == "putNA") {
    duplicatedIds <- selRes[duplicated(selRes[, 1]), 1]
    selRes <- selRes[!selRes[, 1] %in% duplicatedIds, ] 
    }
  return(selRes[match(ids, selRes[, 1]), 2]) }

res$symbol <- convertIDs(row.names(res), "ENSEMBL", "SYMBOL", org.Hs.eg.db)
res


library("reactome.db")
res$entrez <- convertIDs(row.names(res), "ENSEMBL", "ENTREZID", org.Hs.eg.db)

res2 <- res[res$entrez %in% keys(reactome.db, "ENTREZID") & !is.na(res$pvalue),]
head(res2)

reactomeTable <- AnnotationDbi::select(reactome.db, keys = res2$entrez, keytype = "ENTREZID",
                                       columns = c("ENTREZID", "REACTOMEID"))
head(reactomeTable)

incm <- do.call(rbind, with(reactomeTable, tapply(ENTREZID, factor(REACTOMEID),
                                                  function(x) res2$entrez %in% x)))
colnames(incm) <- res2$entrez
str(incm)

incm <- incm[rowSums(incm) >= 5, ]

testCategory <- function(reactomeid) {
  isMember <- incm[reactomeid, ]
  data.frame(reactomeid = reactomeid, numGenes = sum(isMember), 
             avgLFC = mean(res2$log2FoldChange[isMember]), 
             strength = sum(res2$log2FoldChange[isMember])/sqrt(sum(isMember)), 
             pvalue = t.test(res2$log2FoldChange[isMember])$p.value, 
             reactomeName = reactomePATHID2NAME[[reactomeid]])
}
testCategory("R-HSA-109581")

reactomeResult <- do.call(rbind, lapply(rownames(incm), testCategory))
reactomeResult$padjust <- p.adjust(reactomeResult$pvalue, "BH")

reactomeResultSignif <- reactomeResult[reactomeResult$padjust < 0.05, ]
reactomeResultSignif[order(reactomeResultSignif$strength), ]

rld <- rlogTransformation(dds)
head(assay(rld))

par(mfrow = c(1, 2))
plot(log2(1 + counts(dds, normalized = TRUE)[, 1:2]), col = "#00000020", pch = 20,
     cex = 0.3)
plot(assay(rld)[, 1:2], col = "#00000020", pch = 20, cex = 0.3)

sampleDists <- dist(t(assay(rld)))
sampleDists

sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(colData(rld)$treatment, colData(rld)$patient,
                                    sep = "-")
library("gplots")
heatmap.2(sampleDistMatrix, trace = "none")

library("RColorBrewer")
colours = colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
heatmap.2(sampleDistMatrix, trace = "none", col = colours)

print(plotPCA(rld, intgroup = c("patient", "treatment")))

library("genefilter")
topVarGenes <- head(order(rowVars(assay(rld)), decreasing = TRUE), 35)

heatmap.2(assay(rld)[topVarGenes, ], scale = "row", trace = "none", dendrogram = "column",
          col = colorRampPalette(rev(brewer.pal(9, "RdBu")))(255))
