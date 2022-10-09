library(survival)
library(survminer)
library(timeROC)
setwd('C:/Users/rog/Desktop/2022031')

riskFile="GSE202203的风险文件.txt"    


#读取风险文件
risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
risk=risk[,c("futime", "fustat", "riskScore")]


bioCol=c('#1E7DC1', '#CC4D65', '#F4A758',"#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")


# ######绘制1 3 5年的ROC曲线######
# ROC_rt=timeROC(T=risk$futime, delta=risk$fustat,
#                marker=risk$riskScore, cause=1,
#                weighting='aalen',
#                times=c(1095,1825，2555), ROC=TRUE)

result <-with(risk, timeROC(T=risk$futime,
                            delta=risk$fustat,
                            marker=risk$riskScore,
                            cause=1,
                            times=c(365,1095,1825)))
#identical(c(result$TP[,1],result$TP[,2],result$TP[,3]),as.numeric(result$TP))
dat = data.frame(fpr = as.numeric(result$FP),
                 tpr = as.numeric(result$TP),
                 time = rep(as.factor(c(365,1095,1825)),each = nrow(result$TP)))

library(ggplot2)
ggplot() + 
  geom_line(data = dat,aes(x = fpr, y = tpr,color = time),size = 1) + 
  scale_color_manual(name = NULL,values = c('#F4A758', '#CC4D65', '#1E7DC1'),
                     labels = paste0("AUC of ",c(3,5,7),"-y survival: ",
                                     format(round(result$AUC,2),nsmall = 2)))+
  geom_line(aes(x=c(0,1),y=c(0,1)),color = "grey50",linetype=5)+
  theme_bw()+
  theme(panel.grid = element_blank(),
        legend.background = element_rect(linetype = 1, size = 0.2, colour = "black"),
        legend.position = c(0.765,0.125))+
  scale_x_continuous(expand = c(0.005,0.005))+
  scale_y_continuous(expand = c(0.005,0.005))+
  labs(x = "1 - Specificity",
       y = "Sensitivity")+
  coord_fixed() + theme(axis.title = element_text(size = 14),
                        axis.text = element_text(size = 12, colour = "black"),
                        plot.title = element_text(size = 14),
                        legend.title = element_text(size = 12),
                        panel.background = element_rect(fill = NA,
                                                        linetype = "dotdash"), plot.background = element_rect(linetype = "solid"),
                        legend.key = element_rect(fill = NA),
                        legend.position = c(0.74, 0.15)) + theme(legend.background = element_rect(colour = NA),
                                                                 legend.position = c(0.76, 0.12)) + theme(legend.text = element_text(size = 12))+ 
  theme(legend.position = c(0.78, 0.1)) + theme(legend.position = c(0.71, 0.13))+ 
  theme(axis.title = element_text(size = 12.5))




