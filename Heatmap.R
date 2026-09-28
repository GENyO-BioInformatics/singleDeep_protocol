
# Load and install missing packages ---------------------------------------

packages <- c("Hmisc", "pheatmap", "paletteer")
package.check <- lapply(
    packages,
    FUN = function(x) {
        if (!require(x, character.only = TRUE)) {
            install.packages(x, dependencies = TRUE)
            library(x, character.only = TRUE)
        }
    }
)

MCCClust <- read.delim("SjS/results_SjS/Status_clusterResults.tsv", row.names = 1)[,"MCC",drop=F]
cellTypes <- rownames(MCCClust)[order(MCCClust$MCC, decreasing = T)]
MCCClust <- MCCClust[cellTypes,,drop=F]

pseudobulkSjS <- read.delim("SjS/data/training/pseudobulk.tsv", check.names = F, row.names = 1)

phenoSjS <- read.delim("SjS/data/SjS_training/Phenodata.tsv")
cellTypesTop <- cellTypes[1:6] 
clusterContributions <- list()
for (cluster in cellTypesTop) {
    clustContrib <- read.delim(paste0("SjS/results_SjS/gene_contributions/geneContributions_cluster_", cluster, ".tsv"),
                               row.names = 1, check.names = F)
    clusterContributions[[cluster]] <- clustContrib
}

genesSjS <- lapply(clusterContributions, rowMeans)
genesSjS <- do.call(cbind, genesSjS)
rownames(genesSjS) <- rownames(clustContrib)
genesSjS <- genesSjS[,cellTypesTop]

genesSjSRank <- apply(-genesSjS, 2, rank)
rownames(genesSjSRank) = rownames(genesSjS)

selectedGenes <- c()

for (cluster in cellTypesTop) {
    topGenes <- rownames(genesSjS)[order(rank(-abs(genesSjS[,cluster])), decreasing = F)[1:5]]
    selectedGenes <- c(selectedGenes, topGenes)
}

selectedGenes <- names(sort(table(selectedGenes), decreasing = T))
genesSjSRank[(genesSjSRank > 100 & genesSjSRank < (nrow(genesSjSRank) - 100))] <- NA

for (cellType in cellTypesTop) {
    samplesCellType <- grep(paste0("_", cellType), rownames(pseudobulkSjS), value = T)
    exprCellType <- colMeans(pseudobulkSjS[samplesCellType,])
    cutValues = cut2(exprCellType, g=20, onlycuts = T)
    topExprGenes <- names(exprCellType)[exprCellType >= cutValues[20]]
    for (gene in topExprGenes) {
        if (gene %in% rownames(genesSjSRank)) {
            if(!is.na(genesSjSRank[gene, cellType])) {
                if (genesSjSRank[gene, cellType] <= 100) {
                    genesSjSRank[gene, cellType] <- 350
                }
                else if (genesSjSRank[gene, cellType] > 1000){
                    genesSjSRank[gene, cellType] <- (nrow(genesSjSRank) - 350)
                }
            }
        }
    }
}

siglas <- c("proT", "E", "Myo", "SM", "CD4", "IgA")
pheatmap(-genesSjSRank[selectedGenes,], cluster_cols = F, cluster_rows = F,
         color = paletteer_d("beyonce::X39"), border_color = "black", labels_col = siglas,
         angle_col = 90, width = 3.5, height = 5, legend = F, fontsize = 9,
         filename = "SjS/heatmap.png")