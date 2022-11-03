# Version info: R 3.2.3, Biobase 2.30.0, GEOquery 2.40.0, limma 3.26.8
################################################################
#   Differential expression analysis with limma
library(GEOquery)
library(limma)
library(umap)

adjustment_list <- list("fdr","BY","bonferroni", "hochberg", "holm", "hommel", "none")
log_transform_list <- list(TRUE,FALSE)
normalization_list <- list(TRUE,FALSE)
limma_precision_weights_list <- list(TRUE,FALSE)
count <- 1

for (adjust_value in adjustment_list)
{
  for (log_transform_value in log_transform_list)
  {
    for (normalization_value in normalization_list)
    {
      for (limma_precision_weights_value in limma_precision_weights_list)
      {
        print(adjust_value)
        print(log_transform_value)
        print(normalization_value)
        print(limma_precision_weights_value)
        # load series and platform data from GEO
        gset <- getGEO("GSE121239", GSEMatrix =TRUE, AnnotGPL=TRUE)
        if (length(gset) > 1) idx <- grep("GPL13158", attr(gset, "names")) else idx <- 1
        gset <- gset[[idx]]
        
        # make proper column names to match toptable 
        fvarLabels(gset) <- make.names(fvarLabels(gset))
        
        # group membership for all samples
        gsms <- paste0("00000000000000000000111111111111111111111111111111",
                       "11111111111111111111111111111111111111111111111111",
                       "11111111111111111111111111111111111111111111111111",
                       "11111111111111111111111111111111111111111111111111",
                       "11111111111111111111111111111111111111111111111111",
                       "11111111111111111111111111111111111111111111111111",
                       "111111111111")
        sml <- strsplit(gsms, split="")[[1]]
        
        meta_row <- data.frame(count, log_transform_value, normalization_value, limma_precision_weights_value, adjust_value)
        write.table(meta_row, file="/Users/gaskell/Dropbox/Mac/Documents/Cirrone_lab/GRN/gea/meta_data.csv", sep = ",",
                    append = TRUE, quote = FALSE, col.names=FALSE, row.names=FALSE)
        
        if (log_transform_value)
        {
          ex <- exprs(gset)
          ex[which(ex <= 0)] <- NaN
          exprs(gset) <- log2(ex) # log2 transform
        }
        
        if (normalization_value)
        {
          exprs(gset) <- normalizeBetweenArrays(exprs(gset)) # normalize data
        }
        
        # assign samples to groups and set up design matrix
        gs <- factor(sml)
        groups <- make.names(c("healthy","SLE"))
        levels(gs) <- groups
        gset$group <- gs
        design <- model.matrix(~group + 0, gset)
        colnames(design) <- levels(gs)
        
        nall <- nrow(gset)
        gset <- gset[complete.cases(exprs(gset)), ]
        
        if (limma_precision_weights_value)
        {
          # calculate precision weights and show plot of mean-variance trend
          v <- vooma(gset, design, plot=T)
          # OR weights by group
          # v <- voomaByGroup(gset, group=groups, design, plot=T, cex=0.1, pch=".", col=1:nlevels(gs))
          v$genes <- fData(gset) # attach gene annotations
        }
        
        # fit linear model
        fit  <- lmFit(v)
        
        # set up contrasts of interest and recalculate model coefficients
        cts <- c(paste(groups[1],"-",groups[2],sep=""))
        cont.matrix <- makeContrasts(contrasts=cts, levels=design)
        fit2 <- contrasts.fit(fit, cont.matrix)
        
        # compute statistics and table of top significant genes
        fit2 <- eBayes(fit2, 0.01)
        tT <- topTable(fit2, adjust=adjust_value, sort.by="P", number=Inf)
        
        tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","t","B","logFC","Gene.symbol","Gene.title"))
        # write.table(tT, file=stdout(), row.names=F, sep="\t")
        save_path <- paste0("/Users/gaskell/Dropbox/Mac/Documents/Cirrone_lab/GRN/gea/","dataset_",count,".csv")
        write.csv(tT,save_path, row.names = TRUE)
        count = count + 1
      }
    }
  }
}

# if (log_transfrom)
# {
#   exprs(gset) <- log2(ex) # log2 transform
# }
# 
# if (normalization)
# {
#   exprs(gset) <- normalizeBetweenArrays(exprs(gset)) # normalize data
# }
# 
# 
# # assign samples to groups and set up design matrix
# gs <- factor(sml)
# groups <- make.names(c("healthy","SLE"))
# levels(gs) <- groups
# gset$group <- gs
# design <- model.matrix(~group + 0, gset)
# colnames(design) <- levels(gs)
# 
# nall <- nrow(gset)
# gset <- gset[complete.cases(exprs(gset)), ]
# 
# if (limma_precision_weights)
# {
#   # calculate precision weights and show plot of mean-variance trend
#   v <- vooma(gset, design, plot=T)
#   # OR weights by group
#   # v <- voomaByGroup(gset, group=groups, design, plot=T, cex=0.1, pch=".", col=1:nlevels(gs))
#   v$genes <- fData(gset) # attach gene annotations
# }
# 
# # fit linear model
# fit  <- lmFit(v)
# 
# # set up contrasts of interest and recalculate model coefficients
# cts <- c(paste(groups[1],"-",groups[2],sep=""))
# cont.matrix <- makeContrasts(contrasts=cts, levels=design)
# fit2 <- contrasts.fit(fit, cont.matrix)
# 
# # compute statistics and table of top significant genes
# fit2 <- eBayes(fit2, 0.01)
# tT <- topTable(fit2, adjust="hommel", sort.by="P", number=Inf)
# 
# tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","t","B","logFC","Gene.symbol","Gene.title"))
# # write.table(tT, file=stdout(), row.names=F, sep="\t")
# write.csv(tT,"/Users/gaskell/Dropbox/Mac/Documents/Cirrone_lab/GRN/gea/dataset_2.csv", row.names = TRUE)

# # Visualize and quality control test results.
# # Build histogram of P-values for all genes. Normal test
# # assumption is that most genes are not differentially expressed.
# tT2 <- topTable(fit2, adjust="hommel", sort.by="B", number=Inf)
# hist(tT2$adj.P.Val, col = "grey", border = "white", xlab = "P-adj",
#      ylab = "Number of genes", main = "P-adj value distribution")
# 
# # summarize test results as "up", "down" or "not expressed"
# dT <- decideTests(fit2, adjust.method="hommel", p.value=0.05)
# 
# # Venn diagram of results
# vennDiagram(dT, circle.col=palette())
# 
# # create Q-Q plot for t-statistic
# t.good <- which(!is.na(fit2$F)) # filter out bad probes
# qqt(fit2$t[t.good], fit2$df.total[t.good], main="Moderated t statistic")
# 
# # volcano plot (log P-value vs log fold change)
# colnames(fit2) # list contrast names
# ct <- 1        # choose contrast of interest
# volcanoplot(fit2, coef=ct, main=colnames(fit2)[ct], pch=20,
#             highlight=length(which(dT[,ct]!=0)), names=rep('+', nrow(fit2)))
# 
# # MD plot (log fold change vs mean log expression)
# # highlight statistically significant (p-adj < 0.05) probes
# plotMD(fit2, column=ct, status=dT[,ct], legend=F, pch=20, cex=1)
# abline(h=0)
# 
# ################################################################
# # General expression data analysis
# ex <- exprs(gset)
# 
# # box-and-whisker plot
# dev.new(width=3+ncol(gset)/6, height=5)
# ord <- order(gs)  # order samples by group
# palette(c("#1B9E77", "#7570B3", "#E7298A", "#E6AB02", "#D95F02",
#           "#66A61E", "#A6761D", "#B32424", "#B324B3", "#666666"))
# par(mar=c(7,4,2,1))
# title <- paste ("GSE121239", "/", annotation(gset), sep ="")
# boxplot(ex[,ord], boxwex=0.6, notch=T, main=title, outline=FALSE, las=2, col=gs[ord])
# legend("topleft", groups, fill=palette(), bty="n")
# dev.off()
# 
# # expression value distribution
# par(mar=c(4,4,2,1))
# title <- paste ("GSE121239", "/", annotation(gset), " value distribution", sep ="")
# plotDensities(ex, group=gs, main=title, legend ="topright")
# 
# # UMAP plot (dimensionality reduction)
# ex <- na.omit(ex) # eliminate rows with NAs
# ex <- ex[!duplicated(ex), ]  # remove duplicates
# ump <- umap(t(ex), n_neighbors = 15, random_state = 123)
# par(mar=c(3,3,2,6), xpd=TRUE)
# plot(ump$layout, main="UMAP plot, nbrs=15", xlab="", ylab="", col=gs, pch=20, cex=1.5)
# legend("topright", inset=c(-0.15,0), legend=levels(gs), pch=20,
#        col=1:nlevels(gs), title="Group", pt.cex=1.5)
# library("maptools")  # point labels without overlaps
# pointLabel(ump$layout, labels = rownames(ump$layout), method="SANN", cex=0.6)

