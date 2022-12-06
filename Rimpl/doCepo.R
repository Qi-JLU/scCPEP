
args <- commandArgs(trailingOnly = TRUE)
exprsMat_path <- args[1]
cellTypes_path <- args[2]
topN <- 50
pSig <- 0.001

exprsMat <- read.csv(exprsMat_path, row.names=1)
cellTypes <- readLines(cellTypes_path)

tt <- Cepo::Cepo(as.matrix(exprsMat), cellTypes, exprsPct = 0.05)
print(typeof(tt))
res <- Reduce(union, Cepo::topGenes(tt, n=topN))

write(res, file.path(dirname(exprsMat_path), "sub_features.txt"))