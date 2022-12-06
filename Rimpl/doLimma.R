suppressMessages(library(limma))
args <- commandArgs(trailingOnly = TRUE)
exprsMat_path <- args[1]
cellTypes_path <- args[2]
topN <- 50
pSig <- 0.001

exprsMat <- read.csv(exprsMat_path, row.names=1)
cellTypes <- readLines(cellTypes_path)

exprs_pct <- 0.05

## doLimma
cellTypes <- droplevels(as.factor(cellTypes))
tt <- list()
for (i in seq_len(nlevels(cellTypes))) {
    tmp_celltype <- (ifelse(cellTypes == levels(cellTypes)[i], 1, 0))
    design <- stats::model.matrix(~tmp_celltype)


    meanExprs <- do.call(cbind, lapply(c(0,1), function(i){
        Matrix::rowMeans(exprsMat[, tmp_celltype == i, drop = FALSE])
    }))

    meanPct <- do.call(cbind, lapply(c(0,1), function(i){
        Matrix::rowSums(exprsMat[, tmp_celltype == i,
                                 drop = FALSE] > 0)/sum(tmp_celltype == i)
    }))

    keep <- meanPct[,2] > exprs_pct

    y <- methods::new("EList")
    y$E <- exprsMat[keep, ]
    fit <- limma::lmFit(y, design = design)
    fit <- limma::eBayes(fit, trend = TRUE, robust = TRUE)
    tt[[i]] <- limma::topTable(fit, n = Inf, adjust.method = "BH", coef = 2)



    if (!is.null(tt[[i]]$ID)) {
        tt[[i]] <- tt[[i]][!duplicated(tt[[i]]$ID),]
        rownames(tt[[i]]) <- tt[[i]]$ID
    }

    tt[[i]]$meanExprs.1 <- meanExprs[rownames(tt[[i]]), 1]
    tt[[i]]$meanExprs.2 <- meanExprs[rownames(tt[[i]]), 2]
    tt[[i]]$meanPct.1 <- meanPct[rownames(tt[[i]]), 1]
    tt[[i]]$meanPct.2 <- meanPct[rownames(tt[[i]]), 2]
}
print(tt)
res <- Reduce(union, lapply(tt, function(t)
            rownames(t[t$logFC > 0 & (t$meanPct.2 - t$meanPct.1) > 0.05 &
                           t$adj.P.Val < pSig,])[seq_len(topN)]))

write(res, file.path(dirname(exprsMat_path), "sub_features.txt"))

