### R code from vignette source 'viper.Rnw'

###################################################
### code chunk number 1: viper.Rnw:68-72 (eval = FALSE)
###################################################
## source("http://bioconductor.org/biocLite.R")
## biocLite("mixtools")
## biocLite("bcellViper")
## biocLite("viper")


###################################################
### code chunk number 2: viper.Rnw:78-79
###################################################
library(viper)


###################################################
### code chunk number 3: viper.Rnw:94-97
###################################################
data(bcellViper, package="bcellViper")
adjfile <- system.file("aracne", "bcellaracne.adj", package = "bcellViper")
regul <- aracne2regulon(adjfile, dset, verbose = FALSE)


###################################################
### code chunk number 4: viper.Rnw:99-100
###################################################
print(regul)


###################################################
### code chunk number 5: viper.Rnw:117-118
###################################################
signature <- rowTtest(dset, "description", c("CB", "CC"), "N")


###################################################
### code chunk number 6: viper.Rnw:124-126
###################################################
signature <- (qnorm(signature$p.value/2, lower.tail = FALSE) * 
                sign(signature$statistic))[, 1]


###################################################
### code chunk number 7: viper.Rnw:134-136
###################################################
nullmodel <- ttestNull(dset, "description", c("CB", "CC"), "N", per = 1000,
                       repos = TRUE, verbose = FALSE)


###################################################
### code chunk number 8: viper.Rnw:144-145
###################################################
regulon


###################################################
### code chunk number 9: viper.Rnw:150-151
###################################################
mrs <- msviper(signature, regulon, nullmodel, verbose = FALSE)


###################################################
### code chunk number 10: viper.Rnw:156-157
###################################################
summary(mrs)


###################################################
### code chunk number 11: msviper
###################################################
plot(mrs, cex = .7)


###################################################
### code chunk number 12: viper.Rnw:179-181
###################################################
mrs <- ledge(mrs)
summary(mrs)


###################################################
### code chunk number 13: viper.Rnw:191-193
###################################################
signature <- bootstrapTtest(dset, "description", c("CB", "CC"), "N", verbose = FALSE)
mrs <- msviper(signature, regulon, nullmodel, verbose = FALSE)


###################################################
### code chunk number 14: viper.Rnw:197-198
###################################################
mrs <- bootstrapmsviper(mrs, "mode")


###################################################
### code chunk number 15: bsmsviper
###################################################
plot(mrs, cex = .7)


###################################################
### code chunk number 16: viper.Rnw:219-220
###################################################
mrshadow <- shadow(mrs, regulators = 25, verbose = FALSE)


###################################################
### code chunk number 17: viper.Rnw:225-226
###################################################
summary(mrshadow)


###################################################
### code chunk number 18: viper.Rnw:235-236
###################################################
mrs <- msviperCombinatorial(mrs, regulators = 25, verbose = FALSE)


###################################################
### code chunk number 19: viper.Rnw:240-241
###################################################
mrs <- msviperSynergy(mrs, verbose = FALSE)


###################################################
### code chunk number 20: synmsviper
###################################################
summary(mrs)
plot(mrs, 25, cex = .7)


###################################################
### code chunk number 21: viper.Rnw:263-264
###################################################
vpres <- viper(dset, regulon, verbose = FALSE)


###################################################
### code chunk number 22: viper.Rnw:268-269
###################################################
dim(vpres)


###################################################
### code chunk number 23: viper.Rnw:273-276
###################################################
tmp <- rowTtest(vpres, "description", c("CB", "CC"), "N")
data.frame(Gene = rownames(tmp$p.value), t = round(tmp$statistic, 2),
"p-value" = signif(tmp$p.value, 3))[order(tmp$p.value)[1:10], ]


###################################################
### code chunk number 24: viper.Rnw:287-289
###################################################
vpsig <- viperSignature(dset, "description", "N", verbose = FALSE)
vpres <- viper(vpsig, regulon, verbose = FALSE)


###################################################
### code chunk number 25: euviper
###################################################
pos <- pData(vpres)[["description"]] %in% c("M", "CB", "CC")
d1 <- exprs(vpres)[, pos]
colnames(d1) <- pData(vpres)[["description"]][pos]
dd <- dist(t(d1), method = "euclidean")
heatmap(as.matrix(dd), Rowv = as.dendrogram(hclust(dd, method = "average")), symm = T)


###################################################
### code chunk number 26: viper.Rnw:313-314
###################################################
dd <- viperSimilarity(d1)


###################################################
### code chunk number 27: sigviper
###################################################
heatmap(as.matrix(as.dist(dd)), Rowv = as.dendrogram(hclust(as.dist(dd),
method = "average")), symm = T)


