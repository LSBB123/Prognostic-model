library(survival)
library(survminer)
library(dplyr)
library(glmnet)
library(GEOquery)
library(hgu133plus2.db)
library(dplyr)
library(tibble)
library(limma)
library(data.table)
library(devtools)
library(GEOmirror)
library(stringr)
library(limma)
library(ggthemes)
#library('rtracklayer')
setwd('C:/Users/rog/Desktop/2022031')
windowsFonts(arial=windowsFont("arial"))

#gtf <- rtracklayer::import('GSE202203_UCSC_hg38_knownGenes_22sep2014.gtf')
#gtf <- as.data.frame(gtf)

data <- read.table(file = 'GSE202203_RawCounts_gene_3207.txt', header = T, sep = '\t', check.names = F, row.names = 1)
gene <- read.table(file = 'geneCoef1.txt', header = T, sep = '\t')
data <- data[rownames(data)%in%gene$Gene, ]
data <- log2(data+1)
anno <- read.table(file = 'cli.txt', sep = '\t', check.names = F, row.names = 1, header = T)
anno_1 <- read.table(file = 'cli1.txt', sep = '\t', check.names = F, row.names = 1, header = T)
anno_z <- cbind(anno, anno_1)
anno_OS <- anno_z[c(3,4,5), ]
anno_OS <- as.data.frame(t(anno_OS))
anno_OS <- anno_OS[anno_OS$`!Sample_characteristics_ch5` == 'clinical groups: TNBC', ]
names(anno_OS)[names(anno_OS) == '!Sample_characteristics_ch3'] <- 'fustat'
names(anno_OS)[names(anno_OS) == '!Sample_characteristics_ch4'] <- 'futime'
anno_OS$fustat <- str_sub(anno_OS$fustat,25,-1)
anno_OS$futime <- str_sub(anno_OS$futime,24,-1)
anno_OS <- anno_OS[, c(1,2)]
data <- as.data.frame(t(data))
samSample=intersect(row.names(data), row.names(anno_OS))
data=data[samSample,,drop=F]
anno_OS=anno_OS[samSample,,drop=F]
rt=cbind(data, anno_OS)
rt <- as.matrix(rt)
dimnames <- list(rownames(rt), colnames(rt))
rt <- matrix(as.numeric(as.matrix(rt))
                  , nrow=nrow(rt)
                  , ncol = ncol(rt)
                  , dimnames = dimnames
)
rt <- rt[, c(7,6,4,3,5,1,2)]
exprSet <- rt
outdata <- data.frame(matrix(nrow=320,ncol=1))
for (i in 1: 5) {
  a <- exprSet[ ,c(2+i)]
  b <- gene[i, 2]
  b <- as.numeric(b)
  c <- a*b
  outdata <- cbind(outdata, c)
}
outdata <- outdata[, -1]
outdata <- cbind(outdata,Total=rowSums(outdata))
rt <- cbind(exprSet, outdata$Total)
rt <- as.data.frame(rt)
names(rt)[names(rt) == "V8"] <- 'riskScore'
write.table(rt, file = 'GSE202203的风险文件.txt', sep="\t", quote=F)

rt <- read.table(file = 'GSE202203的风险文件.txt', check.names = F, sep = '\t', row.names = 1)
rt$futime <- rt$futime/30
# res.cut=surv_cutpoint(rt, time='futime', event='fustat', variables='riskScore')
# cutoff=as.numeric(res.cut$cutpoint[1])
# print(cutoff)
data=rt
Type=ifelse(data[,'riskScore']<= median(data$riskScore), "Low", "High")
table(Type)
data$group=Type
data$group=factor(data$group, levels=c("Low", "High"))
diff=survdiff(Surv(futime, fustat) ~ group, data = data)
length=length(levels(factor(data[,"group"])))
pValue=1-pchisq(diff$chisq, df=length-1)

fit <- survfit(Surv(futime, fustat) ~ group, data = data)
bioCol=c('#1E7DC1','#CC4D65',"#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=bioCol[1:length]

p=ggsurvplot(fit, 
             data=data,
             conf.int=F,
             pval=TRUE,
             pval.size=5,
             legend.labs=levels(factor(data[,"group"])),
             legend = c(0.4, 0.27),
             font.legend=14,
             xlab="Time (Months)",
             palette = bioCol,
             pval.method = TRUE, 
             cumevents=F
             , break.time.by = 50
             , ggtheme = theme_few()
             ,font.tickslab = c(16, "black")
             , xlim=c(0,115)
)
p[["plot"]][["theme"]][["axis.title.x"]][["family"]] = 'arial'  #x轴标题
p[["plot"]][["layers"]][[4]][["computed_geom_params"]][["family"]] = 'arial'
p[["plot"]][["layers"]][[5]][["aes_params"]][["size"]] = 5
p[["plot"]][["layers"]][[4]][["computed_geom_params"]][["family"]] = 'italic'
p[["plot"]][["theme"]][["axis.title.y"]][["family"]] = 'arial'
p[["plot"]][["layers"]][[4]][["aes_params"]][["label"]] <- expression(italic('P = 0.0032'))
p[["plot"]][["theme"]][["axis.text.x"]][["size"]] = 14
p[["plot"]][["theme"]][["axis.text.y"]][["size"]] = 14
p[["table"]][["layers"]][[1]][["aes_params"]][['family']] <- 'sans'
p[["table"]][["theme"]][["axis.text"]][["colour"]] <- 'black'
p[["plot"]][["theme"]][["panel.background"]][["fill"]] <- NA
p[["plot"]][["theme"]][["axis.title.y"]][["size"]] <- 16
p[["plot"]][["theme"]][["axis.title.x"]][["size"]] <- 16
p[["plot"]][["theme"]][["axis.ticks"]][["size"]] <- 0.5
p
ggsave(filename = 'C:/Users/rog/Desktop/2022031/GSE202203生存曲线1.pdf',width = 5,height = 4.5)



rt <- read.table(file = 'GSE202203的风险文件.txt', check.names = F, sep = '\t', row.names = 1)
rt$futime <- rt$futime/30
# res.cut=surv_cutpoint(rt, time='futime', event='fustat', variables='riskScore')
# cutoff=as.numeric(res.cut$cutpoint[1])
# print(cutoff)
data=rt
Type=ifelse(data[,'GNLY']<= median(data$GNLY), "Low", "High")
table(Type)
data$group=Type
data$group=factor(data$group, levels=c("Low", "High"))
diff=survdiff(Surv(futime, fustat) ~ group, data = data)
length=length(levels(factor(data[,"group"])))
pValue=1-pchisq(diff$chisq, df=length-1)

fit <- survfit(Surv(futime, fustat) ~ group, data = data)
bioCol=c('#1E7DC1','#CC4D65',"#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=bioCol[1:length]

p=ggsurvplot(fit, 
             data=data,
             conf.int=F,
             pval=TRUE,
             pval.size=5,
             legend.labs=levels(factor(data[,"group"])),
             legend = c(0.4, 0.27),
             font.legend=14,
             xlab="Time (Months)",
             palette = bioCol,
             pval.method = TRUE, 
             cumevents=F
             , break.time.by = 50
             , ggtheme = theme_few()
             ,font.tickslab = c(16, "black")
             , xlim=c(0,115)
)
p[["plot"]][["theme"]][["axis.title.x"]][["family"]] = 'arial'  #x轴标题
p[["plot"]][["layers"]][[4]][["computed_geom_params"]][["family"]] = 'arial'
p[["plot"]][["layers"]][[5]][["aes_params"]][["size"]] = 5
p[["plot"]][["layers"]][[4]][["computed_geom_params"]][["family"]] = 'italic'
p[["plot"]][["theme"]][["axis.title.y"]][["family"]] = 'arial'
p[["plot"]][["layers"]][[4]][["aes_params"]][["label"]] <- expression(italic('P = 0.0255'))
p[["plot"]][["theme"]][["axis.text.x"]][["size"]] = 14
p[["plot"]][["theme"]][["axis.text.y"]][["size"]] = 14
p[["table"]][["layers"]][[1]][["aes_params"]][['family']] <- 'sans'
p[["table"]][["theme"]][["axis.text"]][["colour"]] <- 'black'
p[["plot"]][["theme"]][["panel.background"]][["fill"]] <- NA
p[["plot"]][["theme"]][["axis.title.y"]][["size"]] <- 16
p[["plot"]][["theme"]][["axis.title.x"]][["size"]] <- 16
p[["plot"]][["theme"]][["axis.ticks"]][["size"]] <- 0.5
p
ggsave(filename = 'C:/Users/rog/Desktop/2022031/GNLYGSE202203生存曲线1.pdf',width = 5,height = 4.5)













