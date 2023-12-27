# load linear mixed effect model / stats/ plotting librarires
library(lme4)
library(lmerTest)
library(readxl)
library(emmeans)
library(MuMIn)
library(sjPlot)
library(sjmisc)
library(sjstats)
library(ggplot2)
library(effects)
library(shiny)
library(relaimpo)

# directories
filesep <- .Platform$file.sep
d0 <- getwd()
input.dir <- paste(d0,filesep,"..",filesep,"data",sep = "")
output.dir <- paste(d0,filesep,"..",filesep,"results",sep = "")

signed.data.root <- "ROIs_bilat_QSM_data_PD_v1+v3"
abs.data.root <- "ROIs_bilat_abs_QSM_data_PD_v1+v3"


# do some housekeeping

signed.QSMdata <- read.csv(paste(input.dir,filesep,signed.data.root,".csv",sep = ""))
signed.QSMdata$Sex <- as.factor(signed.QSMdata$Sex)
abs.QSMdata <- read.csv(paste(input.dir,filesep,abs.data.root,".csv",sep = ""))
abs.QSMdata$Sex <- as.factor(abs.QSMdata$Sex)
ROI.names <- colnames(abs.QSMdata)[-1:-15]

a_thr <- 0.05
mc_thr <- a_thr/length(ROI.names)

# get data subsets required
PD.V1.abs.QSMdata <- subset(abs.QSMdata, Timepoint == 1)
PD.V1.signed.QSMdata <- subset(signed.QSMdata, Timepoint == 1)
PD.V3.abs.QSMdata <- subset(abs.QSMdata, Timepoint == 2)
PD.V3.signed.QSMdata <- subset(signed.QSMdata, Timepoint == 2)

##### PD only linear regression V1 QSM against bl CCog ----

for (ROI in ROI.names) {
  # define a linear regression model for each ROI
  # absolute QSM
  lmdl.PD.abs.V1.bl.CCog <- lm(PD.V1.abs.QSMdata[,ROI] ~ Age.baseline + Sex + CCog.baseline, data = PD.V1.abs.QSMdata)
  abs.relimp <- calc.relimp.lm(lmdl.PD.abs.V1.bl.CCog)
  abs.mdl.name <- paste("PD.abs.V1.lmdl.bl.CCog.",ROI,sep="")
  abs.res <- summary(lmdl.PD.abs.V1.bl.CCog)
  abs.mdl.pval <- pf(abs.res$fstatistic[1],abs.res$fstatistic[2],abs.res$fstatistic[3], lower.tail = FALSE)
  abs.eff.pval <- anova(lmdl.PD.abs.V1.bl.CCog)$`Pr(>F)`[3]
  abs.eff.slope <- formatC(coef(lmdl.PD.abs.V1.bl.CCog)[4], format = "e", digits = 2)
  abs.aRsq = abs.res$adj.r.squared
  abs.mRsq <- abs.res$r.squared
  abs.rel.Rsq <- abs.relimp@lmg[4]
  # signed QSM
  lmdl.PD.signed.V1.bl.CCog <- lm(PD.V1.signed.QSMdata[,ROI] ~ Age.baseline + Sex + CCog.baseline, data = PD.V1.signed.QSMdata)
  signed.relimp <- calc.relimp.lm(lmdl.PD.signed.V1.bl.CCog)
  signed.mdl.name <- paste("PD.signed.V1.lmdl.bl.CCog.",ROI,sep="")
  signed.res <- summary(lmdl.PD.signed.V1.bl.CCog)
  signed.mdl.pval <- pf(signed.res$fstatistic[1],signed.res$fstatistic[2],signed.res$fstatistic[3], lower.tail = FALSE)
  signed.eff.pval <- anova(lmdl.PD.signed.V1.bl.CCog)$`Pr(>F)`[3]
  signed.eff.slope <- formatC(coef(lmdl.PD.signed.V1.bl.CCog)[4], format = "e", digits = 2)
  signed.aRsq = signed.res$adj.r.squared
  signed.mRsq <- signed.res$r.squared
  signed.rel.Rsq <- signed.relimp@lmg[4]
  if (signed.eff.pval < 1 && signed.mdl.pval < 1){
    heading <- paste0("abs_p = ", format.pval(abs.eff.pval, digits = 3), ", abs_beta = ", abs.eff.slope, ", signed_p = ", format.pval(signed.eff.pval, digits = 3), ", signed_beta = ", signed.eff.slope)
    print(abs.mdl.name)
    print(anova(lmdl.PD.abs.V1.bl.CCog))
    print(signed.mdl.name)
    print(anova(lmdl.PD.signed.V1.bl.CCog))
    abs.mdl.effects <- effect(term = "CCog.baseline", mod = lmdl.PD.abs.V1.bl.CCog, xlevels = 20)
    abs.df.mdl <- as.data.frame(abs.mdl.effects)
    signed.mdl.effects <- effect(term = "CCog.baseline", mod = lmdl.PD.signed.V1.bl.CCog, xlevels = 20)
    signed.df.mdl <- as.data.frame(signed.mdl.effects)
    png.name <- paste("PD.abs.signed.V1.lmdl.bl.CCog.",ROI,sep="")
    png(filename = paste(output.dir,filesep, "plot_",png.name,".png",sep=""), 
        width = 650, height = 600, units = "px")
    plot.mdl <- ggplot() + theme_minimal() +
      geom_point(data = PD.V1.abs.QSMdata, aes(CCog.baseline, PD.V1.abs.QSMdata[,ROI]), color = "red", size = 3.5, stroke=0,  shape = 16, alpha=0.5) +
      geom_line(data = abs.df.mdl, aes(x=CCog.baseline, y=fit), color = "red", linewidth=0.8) +
      geom_ribbon(data = abs.df.mdl, aes(x=CCog.baseline, ymin=lower, ymax=upper), alpha = 0.2, fill="red") +
      geom_point(data = PD.V1.signed.QSMdata, aes(CCog.baseline, PD.V1.signed.QSMdata[,ROI]), color = "blue", size = 3.5, stroke=0,  shape = 16, alpha=0.5) +
      geom_line(data = signed.df.mdl, aes(x=CCog.baseline, y=fit), color = "blue", linewidth=0.8) +
      geom_ribbon(data = signed.df.mdl, aes(x=CCog.baseline, ymin=lower, ymax=upper), alpha = 0.2, fill="blue") +
      labs(x="combined cognitive score (baseline)", y=paste(ROI, "baseline")) +
      theme(axis.title = element_text(size = 20), axis.text = element_text(size=25)) +
      ggtitle(ggtitle(heading))
    print(plot.mdl)
    dev.off()
    
  }
  assign(abs.mdl.name,lmdl.PD.abs.V1.bl.CCog)
  assign(signed.mdl.name,lmdl.PD.signed.V1.bl.CCog)
}

##### PD only linear regression V1 QSM against bl UPDRS.III ----

for (ROI in ROI.names) {
  # define a linear regression model for each ROI
  # absolute QSM
  lmdl.PD.abs.V1.bl.UPDRS.III <- lm(PD.V1.abs.QSMdata[,ROI] ~ Age.baseline + Sex + UPDRS.III.baseline, data = PD.V1.abs.QSMdata)
  abs.relimp <- calc.relimp.lm(lmdl.PD.abs.V1.bl.UPDRS.III)
  abs.mdl.name <- paste("PD.abs.V1.lmdl.bl.UPDRS.III.",ROI,sep="")
  abs.res <- summary(lmdl.PD.abs.V1.bl.UPDRS.III)
  abs.mdl.pval <- pf(abs.res$fstatistic[1],abs.res$fstatistic[2],abs.res$fstatistic[3], lower.tail = FALSE)
  abs.eff.pval <- anova(lmdl.PD.abs.V1.bl.UPDRS.III)$`Pr(>F)`[3]
  abs.eff.slope <- formatC(coef(lmdl.PD.abs.V1.bl.UPDRS.III)[4], format = "e", digits = 2)
  abs.aRsq = abs.res$adj.r.squared
  abs.mRsq <- abs.res$r.squared
  abs.rel.Rsq <- abs.relimp@lmg[4]
  # signed QSM
  lmdl.PD.signed.V1.bl.UPDRS.III <- lm(PD.V1.signed.QSMdata[,ROI] ~ Age.baseline + Sex + UPDRS.III.baseline, data = PD.V1.signed.QSMdata)
  signed.relimp <- calc.relimp.lm(lmdl.PD.signed.V1.bl.UPDRS.III)
  signed.mdl.name <- paste("PD.signed.V1.lmdl.bl.UPDRS.III.",ROI,sep="")
  signed.res <- summary(lmdl.PD.signed.V1.bl.UPDRS.III)
  signed.mdl.pval <- pf(signed.res$fstatistic[1],signed.res$fstatistic[2],signed.res$fstatistic[3], lower.tail = FALSE)
  signed.eff.pval <- anova(lmdl.PD.signed.V1.bl.UPDRS.III)$`Pr(>F)`[3]
  signed.eff.slope <- formatC(coef(lmdl.PD.signed.V1.bl.UPDRS.III)[4], format = "e", digits = 2)
  signed.aRsq = signed.res$adj.r.squared
  signed.mRsq <- signed.res$r.squared
  signed.rel.Rsq <- signed.relimp@lmg[4]
  if (signed.eff.pval < 1 && signed.mdl.pval < 1){
    heading <- paste0("abs_p = ", format.pval(abs.eff.pval, digits = 3), ", abs_beta = ", abs.eff.slope, ", signed_p = ", format.pval(signed.eff.pval, digits = 3), ", signed_beta = ", signed.eff.slope)
    print(abs.mdl.name)
    print(anova(lmdl.PD.abs.V1.bl.UPDRS.III))
    print(signed.mdl.name)
    print(anova(lmdl.PD.signed.V1.bl.UPDRS.III))
    abs.mdl.effects <- effect(term = "UPDRS.III.baseline", mod = lmdl.PD.abs.V1.bl.UPDRS.III, xlevels = 20)
    abs.df.mdl <- as.data.frame(abs.mdl.effects)
    signed.mdl.effects <- effect(term = "UPDRS.III.baseline", mod = lmdl.PD.signed.V1.bl.UPDRS.III, xlevels = 20)
    signed.df.mdl <- as.data.frame(signed.mdl.effects)
    png.name <- paste("PD.abs.signed.V1.lmdl.bl.UPDRS.III.",ROI,sep="")
    png(filename = paste(output.dir,filesep, "plot_",png.name,".png",sep=""), 
        width = 650, height = 600, units = "px")
    plot.mdl <- ggplot() + theme_minimal() +
      geom_point(data = PD.V1.abs.QSMdata, aes(UPDRS.III.baseline, PD.V1.abs.QSMdata[,ROI]), color = "red", size = 3.5, stroke=0,  shape = 16, alpha=0.5) +
      geom_line(data = abs.df.mdl, aes(x=UPDRS.III.baseline, y=fit), color = "red", linewidth=0.8) +
      geom_ribbon(data = abs.df.mdl, aes(x=UPDRS.III.baseline, ymin=lower, ymax=upper), alpha = 0.2, fill="red") +
      geom_point(data = PD.V1.signed.QSMdata, aes(UPDRS.III.baseline, PD.V1.signed.QSMdata[,ROI]), color = "blue", size = 3.5, stroke=0,  shape = 16, alpha=0.5) +
      geom_line(data = signed.df.mdl, aes(x=UPDRS.III.baseline, y=fit), color = "blue", linewidth=0.8) +
      geom_ribbon(data = signed.df.mdl, aes(x=UPDRS.III.baseline, ymin=lower, ymax=upper), alpha = 0.2, fill="blue") +
      labs(x="UPDRS-III (baseline)", y=paste(ROI, "baseline")) +
      theme(axis.title = element_text(size = 20), axis.text = element_text(size=25)) +
      ggtitle(ggtitle(heading))
    print(plot.mdl)
    dev.off()
    
  }
  assign(abs.mdl.name,lmdl.PD.abs.V1.bl.UPDRS.III)
  assign(signed.mdl.name,lmdl.PD.signed.V1.bl.UPDRS.III)
}

##### PD only linear regression V1 QSM against fup CCog ----

for (ROI in ROI.names) {
  # define a linear regression model for each ROI
  # absolute QSM
  lmdl.PD.abs.V1.CCog <- lm(PD.V1.abs.QSMdata[,ROI] ~ Age.baseline + Sex + fup.interval.const + CCog.fup, data = PD.V1.abs.QSMdata)
  abs.relimp <- calc.relimp.lm(lmdl.PD.abs.V1.CCog)
  abs.mdl.name <- paste("PD.abs.V1.lmdl.fup.CCog.",ROI,sep="")
  abs.res <- summary(lmdl.PD.abs.V1.CCog)
  abs.mdl.pval <- pf(abs.res$fstatistic[1],abs.res$fstatistic[2],abs.res$fstatistic[3], lower.tail = FALSE)
  abs.eff.pval <- anova(lmdl.PD.abs.V1.CCog)$`Pr(>F)`[4]
  abs.eff.slope <- formatC(coef(lmdl.PD.abs.V1.CCog)[5], format = "e", digits = 2)
  abs.aRsq = abs.res$adj.r.squared
  abs.mRsq <- abs.res$r.squared
  abs.rel.Rsq <- abs.relimp@lmg[5]
  # signed QSM
  lmdl.PD.signed.V1.CCog <- lm(PD.V1.signed.QSMdata[,ROI] ~ Age.baseline + Sex + fup.interval.const + CCog.fup, data = PD.V1.signed.QSMdata)
  signed.relimp <- calc.relimp.lm(lmdl.PD.signed.V1.CCog)
  signed.mdl.name <- paste("PD.signed.V1.lmdl.fup.CCog.",ROI,sep="")
  signed.res <- summary(lmdl.PD.signed.V1.CCog)
  signed.mdl.pval <- pf(signed.res$fstatistic[1],signed.res$fstatistic[2],signed.res$fstatistic[3], lower.tail = FALSE)
  signed.eff.pval <- anova(lmdl.PD.signed.V1.CCog)$`Pr(>F)`[4]
  signed.eff.slope <- formatC(coef(lmdl.PD.signed.V1.CCog)[5], format = "e", digits = 2)
  signed.aRsq = signed.res$adj.r.squared
  signed.mRsq <- signed.res$r.squared
  signed.rel.Rsq <- signed.relimp@lmg[5]
  if (signed.eff.pval < 1 && signed.mdl.pval < 1){
    heading <- paste0("abs_p = ", format.pval(abs.eff.pval, digits = 3), ", abs_beta = ", abs.eff.slope, ", signed_p = ", format.pval(signed.eff.pval, digits = 3), ", signed_beta = ", signed.eff.slope)
    print(abs.mdl.name)
    print(anova(lmdl.PD.abs.V1.CCog))
    print(signed.mdl.name)
    print(anova(lmdl.PD.signed.V1.CCog))
    abs.mdl.effects <- effect(term = "CCog.fup", mod = lmdl.PD.abs.V1.CCog, xlevels = 20)
    abs.df.mdl <- as.data.frame(abs.mdl.effects)
    signed.mdl.effects <- effect(term = "CCog.fup", mod = lmdl.PD.signed.V1.CCog, xlevels = 20)
    signed.df.mdl <- as.data.frame(signed.mdl.effects)
    png.name <- paste("PD.abs.signed.V1.lmdl.fup.CCog.",ROI,sep="")
    png(filename = paste(output.dir,filesep, "plot_",png.name,".png",sep=""), 
                         width = 650, height = 600, units = "px")
    plot.mdl <- ggplot() + theme_minimal() +
      geom_point(data = PD.V1.abs.QSMdata, aes(CCog.fup, PD.V1.abs.QSMdata[,ROI]), color = "red", size = 3.5, stroke=0,  shape = 16, alpha=0.5) +
      geom_line(data = abs.df.mdl, aes(x=CCog.fup, y=fit), color = "red", linewidth=0.8) +
      geom_ribbon(data = abs.df.mdl, aes(x=CCog.fup, ymin=lower, ymax=upper), alpha = 0.2, fill="red") +
      geom_point(data = PD.V1.signed.QSMdata, aes(CCog.fup, PD.V1.signed.QSMdata[,ROI]), color = "blue", size = 3.5, stroke=0,  shape = 16, alpha=0.5) +
      geom_line(data = signed.df.mdl, aes(x=CCog.fup, y=fit), color = "blue", linewidth=0.8) +
      geom_ribbon(data = signed.df.mdl, aes(x=CCog.fup, ymin=lower, ymax=upper), alpha = 0.2, fill="blue") +
      labs(x="combined cognitive score (follow-up)", y=paste(ROI, "baseline")) +
      theme(axis.title = element_text(size = 20), axis.text = element_text(size=25)) +
      ggtitle(ggtitle(heading))
    print(plot.mdl)
    dev.off()

  }
  assign(abs.mdl.name,lmdl.PD.abs.V1.CCog)
  assign(signed.mdl.name,lmdl.PD.signed.V1.CCog)
}

##### PD only linear regression V1 QSM against fup UPDRS-III ----

for (ROI in ROI.names) {
  # define a linear regression model for each ROI
  # absolute QSM
  lmdl.PD.abs.V1.UPDRS.III <- lm(PD.V1.abs.QSMdata[,ROI] ~ Age.baseline + Sex + fup.interval.const + UPDRS.III.fup, data = PD.V1.abs.QSMdata)
  abs.relimp <- calc.relimp.lm(lmdl.PD.abs.V1.UPDRS.III)
  abs.mdl.name <- paste("PD.abs.V1.lmdl.fup.UPDRS.III.",ROI,sep="")
  abs.res <- summary(lmdl.PD.abs.V1.UPDRS.III)
  abs.mdl.pval <- pf(abs.res$fstatistic[1],abs.res$fstatistic[2],abs.res$fstatistic[3], lower.tail = FALSE)
  abs.eff.pval <- anova(lmdl.PD.abs.V1.UPDRS.III)$`Pr(>F)`[4]
  abs.eff.slope <- formatC(coef(lmdl.PD.abs.V1.UPDRS.III)[5], format = "e", digits = 2)
  abs.aRsq = abs.res$adj.r.squared
  abs.mRsq <- abs.res$r.squared
  abs.rel.Rsq <- abs.relimp@lmg[5]
  # signed QSM
  lmdl.PD.signed.V1.UPDRS.III <- lm(PD.V1.signed.QSMdata[,ROI] ~ Age.baseline + Sex + fup.interval.const + UPDRS.III.fup, data = PD.V1.signed.QSMdata)
  signed.relimp <- calc.relimp.lm(lmdl.PD.signed.V1.UPDRS.III)
  signed.mdl.name <- paste("PD.signed.V1.lmdl.fup.UPDRS.III.",ROI,sep="")
  signed.res <- summary(lmdl.PD.signed.V1.UPDRS.III)
  signed.mdl.pval <- pf(signed.res$fstatistic[1],signed.res$fstatistic[2],signed.res$fstatistic[3], lower.tail = FALSE)
  signed.eff.pval <- anova(lmdl.PD.signed.V1.UPDRS.III)$`Pr(>F)`[4]
  signed.eff.slope <- formatC(coef(lmdl.PD.signed.V1.UPDRS.III)[5], format = "e", digits = 2)
  signed.aRsq = signed.res$adj.r.squared
  signed.mRsq <- signed.res$r.squared
  signed.rel.Rsq <- signed.relimp@lmg[5]
  if (signed.eff.pval < 1 && signed.mdl.pval < 1){
    heading <- paste0("abs_p = ", format.pval(abs.eff.pval, digits = 3), ", abs_beta = ", abs.eff.slope, ", signed_p = ", format.pval(signed.eff.pval, digits = 3), ", signed_beta = ", signed.eff.slope)
    print(abs.mdl.name)
    print(anova(lmdl.PD.abs.V1.UPDRS.III))
    print(signed.mdl.name)
    print(anova(lmdl.PD.signed.V1.UPDRS.III))
    abs.mdl.effects <- effect(term = "UPDRS.III.fup", mod = lmdl.PD.abs.V1.UPDRS.III, xlevels = 20)
    abs.df.mdl <- as.data.frame(abs.mdl.effects)
    signed.mdl.effects <- effect(term = "UPDRS.III.fup", mod = lmdl.PD.signed.V1.UPDRS.III, xlevels = 20)
    signed.df.mdl <- as.data.frame(signed.mdl.effects)
    png.name <- paste("PD.abs.signed.V1.lmdl.fup.UPDRS.III.",ROI,sep="")
    png(filename = paste(output.dir,filesep, "plot_",png.name,".png",sep=""), 
        width = 650, height = 600, units = "px")
    plot.mdl <- ggplot() + theme_minimal() +
      geom_point(data = PD.V1.abs.QSMdata, aes(UPDRS.III.fup, PD.V1.abs.QSMdata[,ROI]), color = "red",size = 3.5, stroke=0,  shape = 16, alpha=0.5) +
      geom_line(data = abs.df.mdl, aes(x=UPDRS.III.fup, y=fit), color = "red",linewidth=0.8) +
      geom_ribbon(data = abs.df.mdl, aes(x=UPDRS.III.fup, ymin=lower, ymax=upper), alpha = 0.2, fill="red") +
      geom_point(data = PD.V1.signed.QSMdata, aes(UPDRS.III.fup, PD.V1.signed.QSMdata[,ROI]), color = "blue",size = 3.5, stroke=0,  shape = 16, alpha=0.5) +
      geom_line(data = signed.df.mdl, aes(x=UPDRS.III.fup, y=fit), color = "blue",linewidth=0.8) +
      geom_ribbon(data = signed.df.mdl, aes(x=UPDRS.III.fup, ymin=lower, ymax=upper), alpha = 0.2, fill="blue") +
      labs(x="UPDRS-III (follow-up)", y=paste(ROI, "baseline")) +
      theme(axis.title = element_text(size = 20), axis.text = element_text(size=25)) +
      ggtitle(ggtitle(heading))
    print(plot.mdl)
    dev.off()
    
  }
  assign(abs.mdl.name,lmdl.PD.abs.V1.UPDRS.III)
  assign(signed.mdl.name,lmdl.PD.signed.V1.UPDRS.III)
}


##### PD only linear regression V3 QSM against fup CCog ----

for (ROI in ROI.names) {
  # define a linear regression model for each ROI
  # absolute QSM
  lmdl.PD.abs.V3.fup.CCog <- lm(PD.V3.abs.QSMdata[,ROI] ~ Age.fup + Sex + CCog.fup, data = PD.V3.abs.QSMdata)
  abs.relimp <- calc.relimp.lm(lmdl.PD.abs.V3.fup.CCog)
  abs.mdl.name <- paste("PD.abs.V3.lmdl.fup.CCog.",ROI,sep="")
  abs.res <- summary(lmdl.PD.abs.V3.fup.CCog)
  abs.mdl.pval <- pf(abs.res$fstatistic[1],abs.res$fstatistic[2],abs.res$fstatistic[3], lower.tail = FALSE)
  abs.eff.pval <- anova(lmdl.PD.abs.V3.fup.CCog)$`Pr(>F)`[3]
  abs.eff.slope <- formatC(coef(lmdl.PD.abs.V3.fup.CCog)[4], format = "e", digits = 2)
  abs.aRsq = abs.res$adj.r.squared
  abs.mRsq <- abs.res$r.squared
  abs.rel.Rsq <- abs.relimp@lmg[4]
  # signed QSM
  lmdl.PD.signed.V3.fup.CCog <- lm(PD.V3.signed.QSMdata[,ROI] ~ Age.fup + Sex + CCog.fup, data = PD.V3.signed.QSMdata)
  signed.relimp <- calc.relimp.lm(lmdl.PD.signed.V3.fup.CCog)
  signed.mdl.name <- paste("PD.signed.V3.lmdl.fup.CCog.",ROI,sep="")
  signed.res <- summary(lmdl.PD.signed.V3.fup.CCog)
  signed.mdl.pval <- pf(signed.res$fstatistic[1],signed.res$fstatistic[2],signed.res$fstatistic[3], lower.tail = FALSE)
  signed.eff.pval <- anova(lmdl.PD.signed.V3.fup.CCog)$`Pr(>F)`[3]
  signed.eff.slope <- formatC(coef(lmdl.PD.signed.V3.fup.CCog)[4], format = "e", digits = 2)
  signed.aRsq = signed.res$adj.r.squared
  signed.mRsq <- signed.res$r.squared
  signed.rel.Rsq <- signed.relimp@lmg[4]
  if (signed.eff.pval < 1 && signed.mdl.pval < 1){
    heading <- paste0("abs_p = ", format.pval(abs.eff.pval, digits = 3), ", abs_beta = ", abs.eff.slope, ", signed_p = ", format.pval(signed.eff.pval, digits = 3), ", signed_beta = ", signed.eff.slope)
    print(abs.mdl.name)
    print(anova(lmdl.PD.abs.V3.fup.CCog))
    print(signed.mdl.name)
    print(anova(lmdl.PD.signed.V3.fup.CCog))
    abs.mdl.effects <- effect(term = "CCog.fup", mod = lmdl.PD.abs.V3.fup.CCog, xlevels = 20)
    abs.df.mdl <- as.data.frame(abs.mdl.effects)
    signed.mdl.effects <- effect(term = "CCog.fup", mod = lmdl.PD.signed.V3.fup.CCog, xlevels = 20)
    signed.df.mdl <- as.data.frame(signed.mdl.effects)
    png.name <- paste("PD.abs.signed.V3.lmdl.fup.CCog.",ROI,sep="")
    png(filename = paste(output.dir,filesep, "plot_",png.name,".png",sep=""), 
        width = 650, height = 600, units = "px")
    plot.mdl <- ggplot() + theme_minimal() +
      geom_point(data = PD.V3.abs.QSMdata, aes(CCog.fup, PD.V1.abs.QSMdata[,ROI]), color = "red", size = 3.5, stroke=0,  shape = 16, alpha=0.5) +
      geom_line(data = abs.df.mdl, aes(x=CCog.fup, y=fit), color = "red", linewidth=0.8) +
      geom_ribbon(data = abs.df.mdl, aes(x=CCog.fup, ymin=lower, ymax=upper), alpha = 0.2, fill="red") +
      geom_point(data = PD.V1.signed.QSMdata, aes(CCog.fup, PD.V1.signed.QSMdata[,ROI]), color = "blue", size = 3.5, stroke=0,  shape = 16, alpha=0.5) +
      geom_line(data = signed.df.mdl, aes(x=CCog.fup, y=fit), color = "blue", linewidth=0.8) +
      geom_ribbon(data = signed.df.mdl, aes(x=CCog.fup, ymin=lower, ymax=upper), alpha = 0.2, fill="blue") +
      labs(x="combined cognitive score (follow-up)", y=paste(ROI, "follow-up")) +
      theme(axis.title = element_text(size = 20), axis.text = element_text(size=25)) +
      ggtitle(ggtitle(heading))
    print(plot.mdl)
    dev.off()
    
  }
  assign(abs.mdl.name,lmdl.PD.abs.V3.fup.CCog)
  assign(signed.mdl.name,lmdl.PD.signed.V3.fup.CCog)
}

##### PD linear regression V3 QSM against fup UPDRS.III ----

for (ROI in ROI.names) {
  # define a linear regression model for each ROI
  # absolute QSM
  lmdl.PD.abs.V3.fup.UPDRS.III <- lm(PD.V3.abs.QSMdata[,ROI] ~ Age.fup + Sex + UPDRS.III.fup, data = PD.V3.abs.QSMdata)
  abs.relimp <- calc.relimp.lm(lmdl.PD.abs.V3.fup.UPDRS.III)
  abs.mdl.name <- paste("PD.abs.V3.lmdl.fup.UPDRS.III.",ROI,sep="")
  abs.res <- summary(lmdl.PD.abs.V3.fup.UPDRS.III)
  abs.mdl.pval <- pf(abs.res$fstatistic[1],abs.res$fstatistic[2],abs.res$fstatistic[3], lower.tail = FALSE)
  abs.eff.pval <- anova(lmdl.PD.abs.V3.fup.UPDRS.III)$`Pr(>F)`[3]
  abs.eff.slope <- formatC(coef(lmdl.PD.abs.V3.fup.UPDRS.III)[4], format = "e", digits = 2)
  abs.aRsq = abs.res$adj.r.squared
  abs.mRsq <- abs.res$r.squared
  abs.rel.Rsq <- abs.relimp@lmg[4]
  # signed QSM
  lmdl.PD.signed.V3.fup.UPDRS.III <- lm(PD.V3.signed.QSMdata[,ROI] ~ Age.fup + Sex + UPDRS.III.fup, data = PD.V3.signed.QSMdata)
  signed.relimp <- calc.relimp.lm(lmdl.PD.signed.V3.fup.UPDRS.III)
  signed.mdl.name <- paste("PD.signed.V3.lmdl.fup.UPDRS.III.",ROI,sep="")
  signed.res <- summary(lmdl.PD.signed.V3.fup.UPDRS.III)
  signed.mdl.pval <- pf(signed.res$fstatistic[1],signed.res$fstatistic[2],signed.res$fstatistic[3], lower.tail = FALSE)
  signed.eff.pval <- anova(lmdl.PD.signed.V3.fup.UPDRS.III)$`Pr(>F)`[3]
  signed.eff.slope <- formatC(coef(lmdl.PD.signed.V3.fup.UPDRS.III)[4], format = "e", digits = 2)
  signed.aRsq = signed.res$adj.r.squared
  signed.mRsq <- signed.res$r.squared
  signed.rel.Rsq <- signed.relimp@lmg[4]
  if (signed.eff.pval < 1 && signed.mdl.pval < 1){
    heading <- paste0("abs_p = ", format.pval(abs.eff.pval, digits = 3), ", abs_beta = ", abs.eff.slope, ", signed_p = ", format.pval(signed.eff.pval, digits = 3), ", signed_beta = ", signed.eff.slope)
    print(abs.mdl.name)
    print(anova(lmdl.PD.abs.V3.fup.UPDRS.III))
    print(signed.mdl.name)
    print(anova(lmdl.PD.signed.V3.fup.UPDRS.III))
    abs.mdl.effects <- effect(term = "UPDRS.III.fup", mod = lmdl.PD.abs.V3.fup.UPDRS.III, xlevels = 20)
    abs.df.mdl <- as.data.frame(abs.mdl.effects)
    signed.mdl.effects <- effect(term = "UPDRS.III.fup", mod = lmdl.PD.signed.V3.fup.UPDRS.III, xlevels = 20)
    signed.df.mdl <- as.data.frame(signed.mdl.effects)
    png.name <- paste("PD.abs.signed.V3.lmdl.fup.UPDRS.III.",ROI,sep="")
    png(filename = paste(output.dir,filesep, "plot_",png.name,".png",sep=""), 
        width = 650, height = 600, units = "px")
    plot.mdl <- ggplot() + theme_minimal() +
      geom_point(data = PD.V3.abs.QSMdata, aes(UPDRS.III.fup, PD.V1.abs.QSMdata[,ROI]), color = "red", size = 3.5, stroke=0,  shape = 16, alpha=0.5) +
      geom_line(data = abs.df.mdl, aes(x=UPDRS.III.fup, y=fit), color = "red", linewidth=0.8) +
      geom_ribbon(data = abs.df.mdl, aes(x=UPDRS.III.fup, ymin=lower, ymax=upper), alpha = 0.2, fill="red") +
      geom_point(data = PD.V1.signed.QSMdata, aes(UPDRS.III.fup, PD.V1.signed.QSMdata[,ROI]), color = "blue", size = 3.5, stroke=0,  shape = 16, alpha=0.5) +
      geom_line(data = signed.df.mdl, aes(x=UPDRS.III.fup, y=fit), color = "blue", linewidth=0.8) +
      geom_ribbon(data = signed.df.mdl, aes(x=UPDRS.III.fup, ymin=lower, ymax=upper), alpha = 0.2, fill="blue") +
      labs(x="UPDRS-III (follow-up)", y=paste(ROI, "follow-up")) +
      theme(axis.title = element_text(size = 20), axis.text = element_text(size=25)) +
      ggtitle(ggtitle(heading))
    print(plot.mdl)
    dev.off()
    
  }
  assign(abs.mdl.name,lmdl.PD.abs.V3.fup.UPDRS.III)
  assign(signed.mdl.name,lmdl.PD.signed.V3.fup.UPDRS.III)
}
