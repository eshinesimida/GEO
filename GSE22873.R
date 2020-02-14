library(limma)
library(GEOquery)


#================================================================
#======================Local processing
library(affy)
#perform mas5 normalization
affy_data = ReadAffy(celfile.path='.')
eset.mas5 = mas5(affy_data)
exprSet.nologs = exprs(eset.mas5)
exprSet = log(exprSet.nologs, 2)  #transform to Log_2 if needed


a=eset.mas5 #
#dat=exprs(eset) #a现在是一个对象，取a这个对象通过看说明书知道要用exprs这个函数
dat <- exprSet
dim(dat)#看一下dat这个矩阵的维度
# GPL6244

dat[1:4,1:4] #查看dat这个矩阵的1至4行和1至4列，逗号前为行，逗号后为列
boxplot(dat,las=2)
pd=pData(a) #通过查看说明书知道取对象a里的临床信息用pData


#==============================================



library(mouse4302.db)
#library(mouse430a2.db)
x <- mouse4302ENTREZID# Get the probe identifiers that are mapped to an ENTREZ Gene ID
mapped_probes <- mappedkeys(x)# Convert to a list
xx <- as.list(x[mapped_probes])

process <- function(x){
  y <- x[[1]]
  y
}
entrez <- sapply(xx,process)

idx <- data.frame(ID = names(xx), ENTREZID = entrez, stringsAsFactors = F)
tail(sort(table(idx$ENTREZID)))

#idx[which(idx$ENTREZID=='170791'),]

#dat <- log(dat, 2)
dat <- data.frame(dat)

dat1 <- dat[rownames(dat) %in% idx$ID,]

ids <- idx[match(rownames(dat), idx$ID),]
ids <- ids[!is.na(ids$ENTREZID),]
tail(sort(table(ids$ENTREZID)))

tmp <- by(dat1,
          ids$ENTREZID,
          function(x) rownames(x)[which.max(rowMeans(x))])

probes <- as.character(tmp)

exprSet <- dat1[rownames(dat1) %in% probes,]

head(ids[match(rownames(exprSet),ids$ID),])

rownames(exprSet) <- ids[match(rownames(exprSet),ids$ID),]$ENTREZID

library(stringr)
####=======================differential expression=====================
colnames(exprSet) <- sub(".CEL.gz", "", colnames(exprSet))

#eset <- gset[[1]]
#pDat <- pd
gpCase <- colnames(exprSet)[12:19]
gpControl <- colnames(exprSet)[1:11]

#names <- colnames(exprs(gset[[1]]))[colnames(exprs(gset[[1]]))%in%c(gpCase, gpControl)]
#exprSet <- exprSet[,  c(gpControl,gpCase)  ]
#pDat <- pData(gset[[1]])
#rownames(pData(gset[[1]]))%in%names

f1 <- factor( rep("case",ncol(exprSet)),  levels = c("case", "control"))
a$description <- f1


#aDFrame <- new("AnnotatedDataFrame",  data = pDat)
#eset <- new('ExpressionSet')
#eset <- new("ExpressionSet", exprs = exprSet,phenoData = aDFrame)

design <- model.matrix(~ description + 0, a)
colnames(design) <- levels(f1)
row.names(design) <- sub(".CEL.gz", "", row.names(design))

design[row.names(design)%in%gpCase, "case"] <- 1
design[row.names(design)%in%gpCase, "control"] <- 0
design[row.names(design)%in%gpControl, "control"] <- 1
design[row.names(design)%in%gpControl, "case"] <- 0

fit <- lmFit(exprSet, design)
#expr_gse <- exprs(eset)
#expr_gse <- data.frame(expr_gse)

cont.matrix <- makeContrasts(case-control, levels = design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)

tT <- topTable(fit2, adjust = "fdr", sort.by = "B", number = 35000)

tT$ID<-rownames(tT)

library("org.Mm.eg.db")
x <- org.Mm.egSYMBOL

mapped_genes <- mappedkeys(x)
xx <- as.list(x[mapped_genes])
myAnnot <- unlist(xx)

myAnnot <- myAnnot[which(names(myAnnot)%in%tT$ID)]
iMat <- match(tT$ID, names(myAnnot))
symbols <- myAnnot[iMat]
names(symbols) <- NULL
head(tT[order(abs(tT$logFC)), ])
tT1 <- cbind(symbols, tT)


#===============================================

tmp = tT1[str_detect(tT1$symbols, "Myd88"), ]
tmp[!is.na(tmp[,1]), ]


KOgeneSymbol <- "Myd88" # gene symbol
KOgeneEntrez <- 17874# gene ID


#=======================================
library(SPIA)
p <- 0.05

tT_m <- tT1[tT1$P.Value<p,]
de <- tT_m$logFC
names(de) <- rownames(tT_m)
all <- rownames(tT1)

spiaRES=spia(de = de, all = all, organism = "mmu",nB = 2000,verbose = T,combine = "fisher",
             data.dir = "./")


load("mmuSPIA.RData")
KO_Flag_Table = sapply(path.info, function(X){any(X$nodes%in%KOgeneEntrez)})
tmpID = spiaRES$ID
KO_Flag_Table = KO_Flag_Table[names(KO_Flag_Table)%in%tmpID]
names(tmpID) = tmpID
KO_Flag_Table = (KO_Flag_Table)[tmpID]
Rank = 1:length(KO_Flag_Table)
spiaRES <- cbind(spiaRES, p.FDR = spiaRES$pGFdr, KO_Flag_Table, Rank)
#TPR_FPR <- getTPR_FPR(spiaRES)


#load pathways
library(HighEdgeS)
toKO = "Myd88"

assign(x = "toKO",toKO,envir = .GlobalEnv)
loadDataObject(KO.gene = toKO)
Organism <- DataObject$Organism
loadKEGG.Objects()

#alpha <- 0.1




getTPR_FPR1 <- function (pathwayesRank, alpha){
  numberTruePW=sum(sapply(kpg, function(X,y){any(nodes(X)%in%y)},y=paste("mmu:",KOgeneEntrez,sep = "")))
  sigPathways = as.vector(pathwayesRank$p.FDR)<alpha
  TP = sum(pathwayesRank$KO_Flag_Table[sigPathways]=="TRUE")
  FP = sum(pathwayesRank$KO_Flag_Table[sigPathways]=="FALSE")
  nonSigPathways = length(kpg) - sum(sigPathways)
  FN = numberTruePW - TP
  TN = length(kpg) - TP - FP - FN
  TPR = (TP) / (TP + FN)
  FPR = (FP) / (FP + TN)
  return(list(TPR=TPR, FPR=FPR))
}
alphas <- seq(0,1,0.05)


TPRs <- c()
FPRs <- c()
for(i in alphas){
  TPR_FPR <- getTPR_FPR1(spiaRES, i)
  TPR <- TPR_FPR$TPR
  FPR <- TPR_FPR$FPR

  TPRs <- c(TPRs,TPR)
  FPRs <- c(FPRs,FPR)
}
AUC_dataframe <- data.frame(p=alphas,TPR=TPRs,FPR=FPRs)
write.csv(AUC_dataframe,'AUC_spia_GSE22873.csv')





