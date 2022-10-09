library(survival)
library(survminer)
library(timeROC)
setwd('C:/Users/rog/Desktop/2022031')

riskFile="GSE202203的风险文件.txt"    
risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
risk$riskScore <- scale(risk$riskScore)
cox <- coxph(Surv(futime, fustat) ~ risk$riskScore, data = risk)
coxSummary = summary(cox)
HR=coxSummary$conf.int[,"exp(coef)"]
HR.95L=coxSummary$conf.int[,"lower .95"]
HR.95H=coxSummary$conf.int[,"upper .95"]
pvalue=coxSummary$coefficients[,"Pr(>|z|)"]
C_index <- coxSummary$concordance

uniTab=data.frame()
for(i in colnames(risk[,3:ncol(risk)])){
  cox <- coxph(Surv(futime, fustat) ~ risk[,i], data = risk)
  coxSummary = summary(cox)
  uniTab=rbind(uniTab,
               cbind(id=i,
                     HR=coxSummary$conf.int[,"exp(coef)"],
                     HR.95L=coxSummary$conf.int[,"lower .95"],
                     HR.95H=coxSummary$conf.int[,"upper .95"],
                     pvalue=coxSummary$coefficients[,"Pr(>|z|)"]
                     #,C_index <- coxSummary$concordance[1]
                     )
  )
}

uniTab=uniTab[as.numeric(uniTab[,"pvalue"])<1,]
risk1=risk[,c("futime", "fustat", as.vector(uniTab[,"id"]))]
multiCox=coxph(Surv(futime, fustat) ~ ., data = risk1)
multiCoxSum=summary(multiCox)
multiTab=data.frame()
multiTab=cbind(
  HR=multiCoxSum$conf.int[,"exp(coef)"],
  HR.95L=multiCoxSum$conf.int[,"lower .95"],
  HR.95H=multiCoxSum$conf.int[,"upper .95"],
  pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
multiTab=cbind(id=row.names(multiTab),multiTab)


risk$riskScore <- scale(risk$riskScore)
cox <- coxph(Surv(futime, fustat) ~ risk$GNLY, data = risk)
coxSummary = summary(cox)
HR=coxSummary$conf.int[,"exp(coef)"]
HR.95L=coxSummary$conf.int[,"lower .95"]
HR.95H=coxSummary$conf.int[,"upper .95"]
pvalue=coxSummary$coefficients[,"Pr(>|z|)"]
C_index <- coxSummary$concordance









