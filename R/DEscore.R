# copyright: Xi Wang (sunlightwang@gmail.com)
## DESeq pipeline w/ permutation for RNASeq_GSEA
#require(DESeq2)
#require(locfit)

nbinomTestForMatrices <- function (countsA, countsB, sizeFactorsA, sizeFactorsB, dispsA, dispsB) {
  kAs <- rowSums(cbind(countsA))
  kBs <- rowSums(cbind(countsB))
  mus <- rowMeans(cbind(t(t(countsA)/sizeFactorsA), t(t(countsB)/sizeFactorsB)))
  fullVarsA <- pmax(mus * sum(sizeFactorsA) + dispsA * mus^2 * 
                    sum(sizeFactorsA^2), mus * sum(sizeFactorsA) * (1 + 1e-08))
  fullVarsB <- pmax(mus * sum(sizeFactorsB) + dispsB * mus^2 * 
                    sum(sizeFactorsB^2), mus * sum(sizeFactorsB) * (1 + 1e-08))
  sumDispsA <- (fullVarsA - mus * sum(sizeFactorsA))/(mus * 
                                                      sum(sizeFactorsA))^2
  sumDispsB <- (fullVarsB - mus * sum(sizeFactorsB))/(mus * 
                                                      sum(sizeFactorsB))^2
  sapply(seq(along = kAs), function(i) {
           if (kAs[i] == 0 & kBs[i] == 0) 
             return(NA)
           ks <- 0:(kAs[i] + kBs[i])
           ps <- dnbinom(ks, mu = mus[i] * sum(sizeFactorsA), size = 1/sumDispsA[i]) * 
             dnbinom(kAs[i] + kBs[i] - ks, mu = mus[i] * sum(sizeFactorsB), 
                     size = 1/sumDispsB[i])
           pobs <- dnbinom(kAs[i], mu = mus[i] * sum(sizeFactorsA), 
                           size = 1/sumDispsA[i]) * dnbinom(kBs[i], mu = mus[i] * 
                           sum(sizeFactorsB), size = 1/sumDispsB[i])
           stopifnot(pobs == ps[kAs[i] + 1])
           if (kAs[i] * sum(sizeFactorsB) < kBs[i] * sum(sizeFactorsA)) 
             numer <- ps[1:(kAs[i] + 1)]
           else numer <- ps[(kAs[i] + 1):length(ps)]
           min(1, 2 * sum(numer)/sum(ps))
                                                      })
}

runDESeq <- function(geneCounts, label) {
  dds <- DESeqDataSetFromMatrix(geneCounts, DataFrame(label), ~ label)
  dds <- estimateSizeFactors(dds)
  dds <- estimateDispersions(dds)
  dds
}

# modified DESeq function to report variances for each group
DENBTest <- function (dds) 
{
stopifnot(is(dds, "DESeqDataSet"))
if (length(levels(colData(dds)$label)) != 2)
  stop("For DESeqDataSet with multivariate conditions, only the GLM-based test can be used.")
if (is.null(dispersions(dds)))
  stop("Call 'estimateDispersions' first.")

colA <- colData(dds)$label == levels(colData(dds)$label)[2]
colB <- colData(dds)$label == levels(colData(dds)$label)[1]
baseMean <- rowMeans( counts( dds, normalized=TRUE )[, colA | colB] ) 
sizeFactorsA <- sizeFactors(dds)[colA]
sizeFactorsB <- sizeFactors(dds)[colB]
dispsA <- pmax(dispersions(dds), 1e-8)
dispsB <- pmax(dispersions(dds), 1e-8)
musA <- rowMeans(t(t(counts(dds)[, colA])/sizeFactorsA))
musB <- rowMeans(t(t(counts(dds)[, colB])/sizeFactorsB))
VarsA <- musA * sum(1/sizeFactorsA) / sum(colA)^2 + dispsA * musA^2 / sum(colA)
VarsB <- musB * sum(1/sizeFactorsB) / sum(colB)^2 + dispsB * musB^2 / sum(colB)
pval <- nbinomTestForMatrices(counts(dds)[, colA],
                              counts(dds)[, colB],
                              sizeFactorsA,
                              sizeFactorsB,
                              dispsA, dispsB)
data.frame(id = rownames(counts(dds)), baseMean = baseMean,
           baseMeanA = musA, VarA = VarsA,
           baseMeanB = musB, VarB = VarsB,
           NBstat = (musA - musB) ^ 2 / (VarsA + VarsB),
           foldChange = musB/musA, log2FoldChange = log2(musB/musA),
           pval = pval, padj = p.adjust(pval, method = "BH"),
           stringsAsFactors = FALSE)
}


DEpermutePval <- function(DEGres, permuteNBstat) {
  times <- ncol(permuteNBstat)
  permutePval <- rowSums(DEGres$NBstat <= permuteNBstat) / times
  data.frame(cbind(DEGres, perm.pval = permutePval, perm.padj = p.adjust(permutePval)), 
             stringsAsFactors = FALSE)  
}

## Added on Sep 12: new version for var 
DENBStat4GSEA <- function(dds) {
  stopifnot(is(dds, "DESeqDataSet"))
  stopifnot(length(levels(colData(dds)$label)) == 2)
  colA <- colData(dds)$label == levels(colData(dds)$label)[2]
  colB <- colData(dds)$label == levels(colData(dds)$label)[1]
  countsA <- counts(dds)[, colA]
  countsB <- counts(dds)[, colB]
  sizeFactorsA <- sizeFactors(dds)[colA]
  sizeFactorsB <- sizeFactors(dds)[colB]

  dispsA <- pmax(dispersions(dds), 1e-8)
  dispsB <- pmax(dispersions(dds), 1e-8)
  musA <- rowMeans(t(t(countsA)/sizeFactorsA))
  musB <- rowMeans(t(t(countsB)/sizeFactorsB))
  VarsA <- musA * sum(1/sizeFactorsA) / sum(colA)^2 + dispsA * musA^2 / sum(colA)
  VarsB <- musB * sum(1/sizeFactorsB) / sum(colB)^2 + dispsB * musB^2 / sum(colB)

  data.frame(id = rownames(counts(dds)),
             baseMeanA = musA, VarA = VarsA,
             baseMeanB = musB, VarB = VarsB,
             NBstat = (musA - musB) ^ 2 / (VarsA + VarsB),
             stringsAsFactors = FALSE)
}

## Added on Sep 12: new version for var
## permutate for GSEA
DENBStatPermut4GSEA <- function(dds, permuteMat) {
  stopifnot( is( dds, "DESeqDataSet" ) )
  times <- ncol(permuteMat)
  n_gene <- nrow(counts(dds))
  #permuteNBstatGene <- matrix(NA_real_, n_gene, times)
  #for(i in 1:times) {
  foreach(i = 1:times, .combine='cbind', .packages = c("DESeq2", "SeqGSEA"))  %dopar% {
    dds@colData$label <- as.factor(permuteMat[,i])
    dds <- estimateDispersions(dds)
    DEGresPerm <- DENBStat4GSEA( dds ) 
    DEGresPerm$NBstat
  }
  #permuteNBstatGene
}

topDEGenes <- function(DEGres, n = 20, sortBy = c("padj", "pval", "perm.pval", "perm.padj", "NBstat", "foldChange")) {
## DEGres is an output from DENBtest and DEpermutePval sequentially. 
  sortBy <- match.arg(sortBy, c("padj", "pval", "perm.pval", "perm.padj", "NBstat", "foldChange"))
  stopifnot(sortBy %in% colnames(DEGres))
  switch(sortBy, padj = {
    res <- DEGres[order(DEGres$padj)[1:n],] 
  }, pval = {
    res <- DEGres[order(DEGres$pval)[1:n],] 
  },  perm.padj = {
    res <- DEGres[order(DEGres$perm.padj)[1:n],] 
  }, perm.pval = {
    res <- DEGres[order(DEGres$perm.pval)[1:n],] 
  }, NBstat = {
    res <- DEGres[order(DEGres$NBstat, decreasing = TRUE)[1:n],]
  }, foldChange = {
    res <- DEGres[order(DEGres$foldChange, decreasing = TRUE)[1:n],]
  })
  res 
}

