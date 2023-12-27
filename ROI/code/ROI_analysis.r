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

# change as appropriate
absQSM <- TRUE

# directories
filesep <- .Platform$file.sep
d0 <- getwd()
input.dir <- paste(d0,filesep,"..",filesep,"data",sep = "")
output.dir <- paste(d0,filesep,"..",filesep,"results",sep = "")

if(absQSM ){
  data.root <- "ROIs_bilat_abs_QSM_data_PD_v1+v3"
} else if(!absQSM) {
  data.root <- "ROIs_bilat_QSM_data_PD_v1+v3"
} 

# do some housekeeping

QSMdata <- read.csv(paste(input.dir,filesep,data.root,".csv",sep = ""))
QSMdata$Sex <- as.factor(QSMdata$Sex)
ROI.names <- colnames(QSMdata)[-1:-15]

a_thr <- 0.05
mc_thr <- a_thr/length(ROI.names)

# get data subsets required
PD.QSMdata = QSMdata;
PD.V1.QSMdata <- subset(QSMdata, Timepoint == 1)
PD.V3.QSMdata <- subset(QSMdata, Timepoint == 2)

##### PD-only growth model ----

stats.PD.QSM.LMM <- matrix(rep(NaN,times=length(ROI.names)*5), ncol = 5, byrow = T)
colnames(stats.PD.QSM.LMM) <- c("marginal R2","conditional R2","beta","p-value","fdr p-value")
rownames(stats.PD.QSM.LMM) <- ROI.names

# loop through ROIs
for (ROI in ROI.names) {
  # define a linear mixed model for each ROI
  mdl.PD = lmer(PD.QSMdata[,ROI] ~ Age.baseline + Sex + fup.interval +
                  (1 | Subject), data = PD.QSMdata, REML = FALSE)
  if (absQSM) {
    mdl.name <- paste("PD.abs.QSM.growth.lmm.",ROI,sep="")
  } else {
    mdl.name <- paste("PD.signed.QSM.growth.lmm.",ROI,sep="")
  }
  res <- anova(mdl.PD)
  # marginal + conditional R squared
  mRsq <- r.squaredGLMM(mdl.PD)[1]
  cRsq <- r.squaredGLMM(mdl.PD)[2]
  beta <- coef(summary(mdl.PD))[4,1]
  pval <- res$`Pr(>F)`[3]
  stats.PD.QSM.LMM[ROI,1] <- round(mRsq,3)
  stats.PD.QSM.LMM[ROI,2] <- round(cRsq,3)
  stats.PD.QSM.LMM[ROI,3] <- beta
  stats.PD.QSM.LMM[ROI,4] <- pval
  if(pval < 1) {
    if (pval < mc_thr){
      heading <- paste0("p = ", format.pval(pval, digits = 3), "mc, Rsq = ", round(mRsq,3))
    } else if (pval < a_thr) {
      heading <- paste0("p = ", format.pval(pval, digits = 3), "*, Rsq = ", round(mRsq,3))
    } else {
      heading <- paste0("p = ", format.pval(pval, digits = 3), ", Rsq = ", round(mRsq,3))
    }
    print(mdl.name)
    print(res)
    mdl.effects <- effect(term = "fup.interval", mod = mdl.PD, xlevels = 20)
    df.mdl <- as.data.frame(mdl.effects)
    png(filename = paste(output.dir,filesep, "plot_",mdl.name,".png",sep=""),
        width = 650, height = 600, units = "px")
    plot.mdl <- ggplot() + theme_minimal() +
      geom_point(data=PD.QSMdata, aes(fup.interval, PD.QSMdata[,ROI]), color = "blue", size = 3.5, stroke=0,  shape = 16, alpha=0.5) +
      geom_line(data=df.mdl, aes(x=fup.interval, y=fit), color="blue",linewidth=0.8) +
      geom_ribbon(data= df.mdl, aes(x=fup.interval, ymin=lower, ymax=upper), alpha= 0.2, fill="blue") +
      labs(x="fup.interval (mos)", y=ROI) +
      theme(axis.title = element_text(size = 20), axis.text = element_text(size=25)) +
      ggtitle(heading)
    print(plot.mdl)
    dev.off()
    
  }
  assign(mdl.name,mdl.PD)
}

stats.PD.QSM.LMM <- as.data.frame(stats.PD.QSM.LMM)
stats.PD.QSM.LMM$`p-value` <- as.numeric(stats.PD.QSM.LMM$`p-value`)
stats.PD.QSM.LMM$`fdr p-value` <- p.adjust(stats.PD.QSM.LMM$`p-value`, method = "fdr")
stats.PD.QSM.LMM$`p-value` <- format.pval(stats.PD.QSM.LMM$`p-value`,digits = 3)
stats.PD.QSM.LMM$`fdr p-value` <- format.pval(stats.PD.QSM.LMM$`fdr p-value`,digits = 3)
csv.name <- paste("LMM.stats.", data.root , sep="")
write.csv(stats.PD.QSM.LMM, file = paste(output.dir,filesep,csv.name,".csv",sep=""))

##### PD only linear regression V1 QSM against bl CCog ----

# create a table to store overall results
stats.PD.V1.QSM.bl.CCog <- matrix(rep(NaN,times=length(ROI.names)*6), ncol = 6, byrow = T)
colnames(stats.PD.V1.QSM.bl.CCog) <- c("adjusted R2","multiple R2", "partial R2","beta","p-value","fdr p-value")
rownames(stats.PD.V1.QSM.bl.CCog) <- ROI.names

for (ROI in ROI.names) {
  # define a linear regression model for each ROI
  lmdl.PD.V1.CCog.baseline <- lm(PD.V1.QSMdata[,ROI] ~ Age.baseline + Sex + CCog.baseline, data = PD.V1.QSMdata)
  relimp <- calc.relimp.lm(lmdl.PD.V1.CCog.baseline)
  if (absQSM) {
    mdl.name <- paste("PD.V1.abs.QSM.lmdl.bl.CCog.",ROI,sep="")
  } else {
    mdl.name <- paste("PD.V1.signed.QSM.lmdl.bl.CCog.",ROI,sep="")
  }
  res <- summary(lmdl.PD.V1.CCog.baseline)
  mdl.pval <- pf(res$fstatistic[1],res$fstatistic[2],res$fstatistic[3], lower.tail = FALSE)
  eff.pval <- anova(lmdl.PD.V1.CCog.baseline)$`Pr(>F)`[3]
  eff.slope <- formatC(coef(lmdl.PD.V1.CCog.baseline)[4], format = "e", digits = 2)
  aRsq <- res$adj.r.squared
  mRsq <- res$r.squared
  rel.Rsq <- relimp@lmg[3]
  stats.PD.V1.QSM.bl.CCog[ROI,1] <- round(aRsq,3)
  stats.PD.V1.QSM.bl.CCog[ROI,2] <- round(mRsq,3)
  stats.PD.V1.QSM.bl.CCog[ROI,3] <- round(rel.Rsq,3)
  stats.PD.V1.QSM.bl.CCog[ROI,4] <- eff.slope
  stats.PD.V1.QSM.bl.CCog[ROI,5] <- eff.pval
  if (eff.pval < 1 && mdl.pval < 1){
    if (eff.pval < mc_thr){
      heading <- paste0("p = ", format.pval(eff.pval, digits = 3), "*mc, aRsq = ", round(aRsq,3), ", lmg = ", round(rel.Rsq,3), ", beta = ", eff.slope)
    } else if (eff.pval < a_thr) {
      heading <- paste0("p = ", format.pval(eff.pval, digits = 3), "*, aRsq = ", round(aRsq,3), ", lmg = ", round(rel.Rsq,3), ", beta = ", eff.slope)
    } else {
      heading <- paste0("p = ", format.pval(eff.pval, digits = 3), ", aRsq = ", round(aRsq,3), ", lmg = ", round(rel.Rsq,3), ", beta = ", eff.slope)
    }
    print(mdl.name)
    print(anova(lmdl.PD.V1.CCog.baseline))
    mdl.effects <- effect(term = "CCog.baseline", mod = lmdl.PD.V1.CCog.baseline, xlevels = 20)
    df.mdl <- as.data.frame(mdl.effects)
    png(filename = paste(output.dir,filesep, "plot_",mdl.name,".png",sep=""))
    plot.mdl <- ggplot() + theme_minimal() +
      geom_point(data = PD.V1.QSMdata, aes(CCog.baseline, PD.V1.QSMdata[,ROI]), color = "orangered") +
      geom_line(data = df.mdl, aes(x=CCog.baseline, y=fit), color = "orangered") +
      geom_ribbon(data = df.mdl, aes(x=CCog.baseline, ymin=lower, ymax=upper), alpha = 0.2, fill="orangered") +
      labs(x="combined cognitive score (baseline)", y=paste(ROI, "baseline")) +
      ggtitle(ggtitle(heading))
    print(plot.mdl)
    dev.off()
    
  }
  assign(mdl.name,lmdl.PD.V1.CCog.baseline)
}

stats.PD.V1.QSM.bl.CCog <- as.data.frame(stats.PD.V1.QSM.bl.CCog)
stats.PD.V1.QSM.bl.CCog$`p-value` <- as.numeric(stats.PD.V1.QSM.bl.CCog$`p-value`)
stats.PD.V1.QSM.bl.CCog$`fdr p-value` <- p.adjust(stats.PD.V1.QSM.bl.CCog$`p-value`, method = "fdr")
stats.PD.V1.QSM.bl.CCog$`p-value` <- format.pval(stats.PD.V1.QSM.bl.CCog$`p-value`,digits = 3)
stats.PD.V1.QSM.bl.CCog$`fdr p-value` <- format.pval(stats.PD.V1.QSM.bl.CCog$`fdr p-value`,digits = 3)
csv.name <- paste("Lmdl.stats.V1.bl.CCog.", data.root , sep = "")
write.csv(stats.PD.V1.QSM.bl.CCog, file = paste(output.dir,filesep,csv.name,".csv",sep=""))
# 


##### PD only linear regression V1 QSM against bl UPDRS-III ----

stats.PD.V1.QSM.bl.UPDRS.III <- matrix(rep(NaN,times=length(ROI.names)*6), ncol = 6, byrow = T)
colnames(stats.PD.V1.QSM.bl.UPDRS.III) <- c("adjusted R2","multiple R2", "partial R2","beta","p-value","fdr p-value")
rownames(stats.PD.V1.QSM.bl.UPDRS.III) <- ROI.names

for (ROI in ROI.names) {
  # define a linear regression model for each ROI
  lmdl.PD.V1.UPDRS.III.baseline <- lm(PD.V1.QSMdata[,ROI] ~ Age.baseline + Sex + UPDRS.III.baseline, data = PD.V1.QSMdata)
  if (absQSM) {
    mdl.name <- paste("PD.V1.abs.QSM.lmdl.bl.UPDRS.III.",ROI,sep="")
  } else {
    mdl.name <- paste("PD.V1.signed.QSM.lmdl.bl.UPDRS.III.",ROI,sep="")
  }
  relimp <- calc.relimp.lm(lmdl.PD.V1.UPDRS.III.baseline)
  res <- summary(lmdl.PD.V1.UPDRS.III.baseline)
  mdl.pval <- pf(res$fstatistic[1],res$fstatistic[2],res$fstatistic[3], lower.tail = FALSE)
  eff.pval <- anova(lmdl.PD.V1.UPDRS.III.baseline)$`Pr(>F)`[3]
  eff.slope <- formatC(coef(lmdl.PD.V1.UPDRS.III.baseline)[4], format = "e", digits = 2)
  aRsq = res$adj.r.squared
  mRsq <- res$r.squared
  rel.Rsq <- relimp@lmg[3]
  stats.PD.V1.QSM.bl.UPDRS.III[ROI,1] <- round(aRsq,3)
  stats.PD.V1.QSM.bl.UPDRS.III[ROI,2] <- round(mRsq,3)
  stats.PD.V1.QSM.bl.UPDRS.III[ROI,3] <- round(rel.Rsq,3)
  stats.PD.V1.QSM.bl.UPDRS.III[ROI,4] <- eff.slope
  stats.PD.V1.QSM.bl.UPDRS.III[ROI,5] <- eff.pval
  if (eff.pval < 1 && mdl.pval < 1){
    if (eff.pval < mc_thr){
      heading <- paste0("p = ", format.pval(eff.pval, digits = 3), "*mc, aRsq = ", round(aRsq,3), ", lmg = ", round(rel.Rsq,3), ", beta = ", eff.slope)
    } else if (eff.pval < a_thr) {
      heading <- paste0("p = ", format.pval(eff.pval, digits = 3), "*, aRsq = ", round(aRsq,3), ", lmg = ", round(rel.Rsq,3), ", beta = ", eff.slope)
    } else {
      heading <- paste0("p = ", format.pval(eff.pval, digits = 3), ", aRsq = ", round(aRsq,3), ", lmg = ", round(rel.Rsq,3), ", beta = ", eff.slope)
    }
    print(mdl.name)
    print(anova(lmdl.PD.V1.UPDRS.III.baseline))
    mdl.effects <- effect(term = "UPDRS.III.baseline", mod = lmdl.PD.V1.UPDRS.III.baseline, xlevels = 20)
    df.mdl <- as.data.frame(mdl.effects)
    png(filename = paste(output.dir,filesep, "plot_",mdl.name,".png",sep=""))
    plot.mdl <- ggplot() + theme_minimal() +
      geom_point(data = PD.V1.QSMdata, aes(UPDRS.III.baseline, PD.V1.QSMdata[,ROI]), color = "darkcyan") +
      geom_line(data = df.mdl, aes(x=UPDRS.III.baseline, y=fit), color = "darkcyan") +
      geom_ribbon(data = df.mdl, aes(x=UPDRS.III.baseline, ymin=lower, ymax=upper), alpha = 0.2, fill="darkcyan") +
      labs(x="UPDRS-III (follow-up)", y=paste(ROI, "baseline")) +
      ggtitle(heading)
    print(plot.mdl)
    dev.off()

  }
  assign(mdl.name,lmdl.PD.V1.UPDRS.III.baseline)
}

stats.PD.V1.QSM.bl.UPDRS.III <- as.data.frame(stats.PD.V1.QSM.bl.UPDRS.III)
stats.PD.V1.QSM.bl.UPDRS.III$`p-value` <- as.numeric(stats.PD.V1.QSM.bl.UPDRS.III$`p-value`)
stats.PD.V1.QSM.bl.UPDRS.III$`fdr p-value` <- p.adjust(stats.PD.V1.QSM.bl.UPDRS.III$`p-value`, method = "fdr")
stats.PD.V1.QSM.bl.UPDRS.III$`p-value` <- format.pval(stats.PD.V1.QSM.bl.UPDRS.III$`p-value`,digits = 3)
stats.PD.V1.QSM.bl.UPDRS.III$`fdr p-value` <- format.pval(stats.PD.V1.QSM.bl.UPDRS.III$`fdr p-value`,digits = 3)
csv.name <- paste("Lmdl.stats.V1.bl.UPDRS.III.", data.root , sep = "")
write.csv(stats.PD.V1.QSM.bl.UPDRS.III, file = paste(output.dir,filesep,csv.name,".csv",sep=""))
# 

##### PD only linear regression V1 QSM against fup CCog ----

# create a table to store overall results
stats.PD.V1.QSM.fup.CCog <- matrix(rep(NaN,times=length(ROI.names)*6), ncol = 6, byrow = T)
colnames(stats.PD.V1.QSM.fup.CCog) <- c("adjusted R2","multiple R2", "partial R2","beta","p-value","fdr p-value")
rownames(stats.PD.V1.QSM.fup.CCog) <- ROI.names

for (ROI in ROI.names) {
  # define a linear regression model for each ROI
  # add in CCog_baseline to the below
  lmdl.PD.V1.CCog <- lm(PD.V1.QSMdata[,ROI] ~ Age.baseline + Sex + fup.interval.const + CCog.fup, data = PD.V1.QSMdata)
  relimp <- calc.relimp.lm(lmdl.PD.V1.CCog)
  if (absQSM) {
    mdl.name <- paste("PD.V1.abs.QSM.lmdl.fup.CCog.",ROI,sep="")
  } else {
    mdl.name <- paste("PD.V1.signed.QSM.lmdl.fup.CCog.",ROI,sep="")
  }
  res <- summary(lmdl.PD.V1.CCog)
  mdl.pval <- pf(res$fstatistic[1],res$fstatistic[2],res$fstatistic[3], lower.tail = FALSE)
  eff.pval <- anova(lmdl.PD.V1.CCog)$`Pr(>F)`[4]
  eff.slope <- formatC(coef(lmdl.PD.V1.CCog)[5], format = "e", digits = 2)
  aRsq = res$adj.r.squared
  mRsq <- res$r.squared
  rel.Rsq <- relimp@lmg[5]
  stats.PD.V1.QSM.fup.CCog[ROI,1] <- round(aRsq,3)
  stats.PD.V1.QSM.fup.CCog[ROI,2] <- round(mRsq,3)
  stats.PD.V1.QSM.fup.CCog[ROI,3] <- round(rel.Rsq,3)
  stats.PD.V1.QSM.fup.CCog[ROI,4] <- eff.slope
  stats.PD.V1.QSM.fup.CCog[ROI,5] <- eff.pval
  if (eff.pval < 1 && mdl.pval < 1){
    if (eff.pval < mc_thr){
      heading <- paste0("p = ", format.pval(eff.pval, digits = 3), "*mc, aRsq = ", round(aRsq,3), ", lmg = ", round(rel.Rsq,3), ", beta = ", eff.slope)
    } else if (eff.pval < a_thr) {
      heading <- paste0("p = ", format.pval(eff.pval, digits = 3), "*, aRsq = ", round(aRsq,3), ", lmg = ", round(rel.Rsq,3), ", beta = ", eff.slope)
    } else {
      heading <- paste0("p = ", format.pval(eff.pval, digits = 3), ", aRsq = ", round(aRsq,3), ", lmg = ", round(rel.Rsq,3), ", beta = ", eff.slope)
    }
    print(mdl.name)
    print(anova(lmdl.PD.V1.CCog))
    mdl.effects <- effect(term = "CCog.fup", mod = lmdl.PD.V1.CCog, xlevels = 20)
    df.mdl <- as.data.frame(mdl.effects)
    png(filename = paste(output.dir,filesep, "plot_",mdl.name,".png",sep=""))
    plot.mdl <- ggplot() + theme_minimal() +
      geom_point(data = PD.V1.QSMdata, aes(CCog.fup, PD.V1.QSMdata[,ROI]), color = "orangered") +
      geom_line(data = df.mdl, aes(x=CCog.fup, y=fit), color = "orangered") +
      geom_ribbon(data = df.mdl, aes(x=CCog.fup, ymin=lower, ymax=upper), alpha = 0.2, fill="orangered") +
      labs(x="combined cognitive score (follow-up)", y=paste(ROI, "baseline")) +
      ggtitle(ggtitle(heading))
    print(plot.mdl)
    dev.off()

  }
  assign(mdl.name,lmdl.PD.V1.CCog)
}

stats.PD.V1.QSM.fup.CCog <- as.data.frame(stats.PD.V1.QSM.fup.CCog)
stats.PD.V1.QSM.fup.CCog$`p-value` <- as.numeric(stats.PD.V1.QSM.fup.CCog$`p-value`)
stats.PD.V1.QSM.fup.CCog$`fdr p-value` <- p.adjust(stats.PD.V1.QSM.fup.CCog$`p-value`, method = "fdr")
stats.PD.V1.QSM.fup.CCog$`p-value` <- format.pval(stats.PD.V1.QSM.fup.CCog$`p-value`,digits = 3)
stats.PD.V1.QSM.fup.CCog$`fdr p-value` <- format.pval(stats.PD.V1.QSM.fup.CCog$`fdr p-value`,digits = 3)
csv.name <- paste("Lmdl.stats.V1.fup.CCog.", data.root , sep = "")
write.csv(stats.PD.V1.QSM.fup.CCog, file = paste(output.dir,filesep,csv.name,".csv",sep=""))
# 

##### PD only linear regression V1 QSM against fup UPDRS-III ----

# create a table to store overall results
stats.PD.V1.QSM.fup.UPDRS.III <- matrix(rep(NaN,times=length(ROI.names)*6), ncol = 6, byrow = T)
colnames(stats.PD.V1.QSM.fup.UPDRS.III) <- c("adjusted R2","multiple R2", "partial R2","beta","p-value","fdr p-value")
rownames(stats.PD.V1.QSM.fup.UPDRS.III) <- ROI.names

for (ROI in ROI.names) {
  # define a linear regression model for each ROI
  # add in UPDRS_III baseline to the below
  lmdl.PD.V1.UPDRS.III <- lm(PD.V1.QSMdata[,ROI] ~ Age.baseline + Sex + fup.interval.const + UPDRS.III.fup, data = PD.V1.QSMdata)
  relimp <- calc.relimp.lm(lmdl.PD.V1.UPDRS.III)
  if (absQSM) {
    mdl.name <- paste("PD.V1.abs.QSM.lmdl.fup.UPDRS.III.",ROI,sep="")
  } else {
    mdl.name <- paste("PD.V1.signed.QSM.lmdl.fup.UPDRS.III.",ROI,sep="")
  }
  res <- summary(lmdl.PD.V1.UPDRS.III)
  mdl.pval <- pf(res$fstatistic[1],res$fstatistic[2],res$fstatistic[3], lower.tail = FALSE)
  eff.pval <- anova(lmdl.PD.V1.UPDRS.III)$`Pr(>F)`[4]
  eff.slope <- formatC(coef(lmdl.PD.V1.UPDRS.III)[5], format = "e", digits = 2)
  aRsq = res$adj.r.squared
  mRsq <- res$r.squared
  rel.Rsq <- relimp@lmg[5]
  stats.PD.V1.QSM.fup.UPDRS.III[ROI,1] <- round(aRsq,3)
  stats.PD.V1.QSM.fup.UPDRS.III[ROI,2] <- round(mRsq,3)
  stats.PD.V1.QSM.fup.UPDRS.III[ROI,3] <- round(rel.Rsq,3)
  stats.PD.V1.QSM.fup.UPDRS.III[ROI,4] <- eff.slope
  stats.PD.V1.QSM.fup.UPDRS.III[ROI,5] <- eff.pval
  if (eff.pval < 1 && mdl.pval < 1){
    if (eff.pval < mc_thr){
      heading <- paste0("p = ", format.pval(eff.pval, digits = 3), "*mc, aRsq = ", round(aRsq,3), ", lmg = ", round(rel.Rsq,3), ", beta = ", eff.slope)
    } else if (eff.pval < a_thr) {
      heading <- paste0("p = ", format.pval(eff.pval, digits = 3), "*, aRsq = ", round(aRsq,3), ", lmg = ", round(rel.Rsq,3), ", beta = ", eff.slope)
    } else {
      heading <- paste0("p = ", format.pval(eff.pval, digits = 3), ", aRsq = ", round(aRsq,3), ", lmg = ", round(rel.Rsq,3), ", beta = ", eff.slope)
    }
    print(mdl.name)
    print(anova(lmdl.PD.V1.UPDRS.III))
    mdl.effects <- effect(term = "UPDRS.III.fup", mod = lmdl.PD.V1.UPDRS.III, xlevels = 20)
    df.mdl <- as.data.frame(mdl.effects)
    png(filename = paste(output.dir,filesep, "plot_",mdl.name,".png",sep=""))
    plot.mdl <- ggplot() + theme_minimal() +
      geom_point(data = PD.V1.QSMdata, aes(UPDRS.III.fup, PD.V1.QSMdata[,ROI]), color = "darkcyan") +
      geom_line(data = df.mdl, aes(x=UPDRS.III.fup, y=fit), color = "darkcyan") +
      geom_ribbon(data = df.mdl, aes(x=UPDRS.III.fup, ymin=lower, ymax=upper), alpha = 0.2, fill="darkcyan") +
      labs(x="UPDRS-III (follow-up)", y=paste(ROI, "baseline")) +
      ggtitle(heading)
    print(plot.mdl)
    dev.off()
  }
  assign(mdl.name,lmdl.PD.V1.UPDRS.III)
}

stats.PD.V1.QSM.fup.UPDRS.III <- as.data.frame(stats.PD.V1.QSM.fup.UPDRS.III)
stats.PD.V1.QSM.fup.UPDRS.III$`p-value` <- as.numeric(stats.PD.V1.QSM.fup.UPDRS.III$`p-value`)
stats.PD.V1.QSM.fup.UPDRS.III$`fdr p-value` <- p.adjust(stats.PD.V1.QSM.fup.UPDRS.III$`p-value`, method = "fdr")
stats.PD.V1.QSM.fup.UPDRS.III$`p-value` <- format.pval(stats.PD.V1.QSM.fup.UPDRS.III$`p-value`,digits = 3)
stats.PD.V1.QSM.fup.UPDRS.III$`fdr p-value` <- format.pval(stats.PD.V1.QSM.fup.UPDRS.III$`fdr p-value`,digits = 3)
csv.name <- paste("Lmdl.stats.V1.fup.UPDRS.III.", data.root , sep = "")
write.csv(stats.PD.V1.QSM.fup.UPDRS.III, file = paste(output.dir,filesep,csv.name,".csv",sep=""))
#  

##### PD only linear regression V3 QSM against fup CCog ----

# create a table to store overall results
stats.PD.V3.QSM.fup.CCog <- matrix(rep(NaN,times=length(ROI.names)*6), ncol = 6, byrow = T)
colnames(stats.PD.V3.QSM.fup.CCog) <- c("adjusted R2","multiple R2", "partial R2","beta","p-value", "fdr p-value")
rownames(stats.PD.V3.QSM.fup.CCog) <- ROI.names

# select just the PD V3 data
PD.V3.QSMdata <- subset(QSMdata,Group == "PD" & Timepoint == 2)

for (ROI in ROI.names) {
  # define a linear regression model for each ROI
  lmdl.PD.V3.CCog <- lm(PD.V3.QSMdata[,ROI] ~ Age.fup + Sex + CCog.fup, data = PD.V3.QSMdata)
  relimp <- calc.relimp.lm(lmdl.PD.V3.CCog)
  if (absQSM) {
    mdl.name <- paste("PD.V3.abs.QSM.lmdl.fup.CCog.",ROI,sep="")
  } else {
    mdl.name <- paste("PD.V3.signed.QSM.lmdl.fup.CCog.",ROI,sep="")
  }
  res <- summary(lmdl.PD.V3.CCog)
  mdl.pval <- pf(res$fstatistic[1],res$fstatistic[2],res$fstatistic[3], lower.tail = FALSE)
  eff.pval <- anova(lmdl.PD.V3.CCog)$`Pr(>F)`[3]
  eff.slope <- formatC(coef(lmdl.PD.V3.CCog)[4], format = "e", digits = 2)
  aRsq = res$adj.r.squared
  mRsq <- res$r.squared
  rel.Rsq <- relimp@lmg[3]
  stats.PD.V3.QSM.fup.CCog[ROI,1] <- round(aRsq,3)
  stats.PD.V3.QSM.fup.CCog[ROI,2] <- round(mRsq,3)
  stats.PD.V3.QSM.fup.CCog[ROI,3] <- round(rel.Rsq,3)
  stats.PD.V3.QSM.fup.CCog[ROI,4] <- eff.slope
  stats.PD.V3.QSM.fup.CCog[ROI,5] <- eff.pval
  if (eff.pval < 1 && mdl.pval < 1){
    if (eff.pval < mc_thr){
      heading <- paste0("p = ", format.pval(eff.pval, digits = 3), "*mc, aRsq = ", round(aRsq,3), ", lmg = ", round(rel.Rsq,3), ", beta = ", eff.slope)
    } else if (eff.pval < a_thr) {
      heading <- paste0("p = ", format.pval(eff.pval, digits = 3), "*, aRsq = ", round(aRsq,3), ", lmg = ", round(rel.Rsq,3), ", beta = ", eff.slope)
    } else {
      heading <- paste0("p = ", format.pval(eff.pval, digits = 3), ", aRsq = ", round(aRsq,3), ", lmg = ", round(rel.Rsq,3), ", beta = ", eff.slope)
    }
    print(mdl.name)
    print(anova(lmdl.PD.V3.CCog))
    mdl.effects <- effect(term = "CCog.fup", mod = lmdl.PD.V3.CCog, xlevels = 20)
    df.mdl <- as.data.frame(mdl.effects)
    png(filename = paste(output.dir,filesep, "plot_",mdl.name,".png",sep=""))
    plot.mdl <- ggplot() + theme_minimal() +
      geom_point(data = PD.V3.QSMdata, aes(CCog.fup, PD.V3.QSMdata[,ROI]), color = "orangered") +
      geom_line(data = df.mdl, aes(x=CCog.fup, y=fit), color = "orangered") +
      geom_ribbon(data = df.mdl, aes(x=CCog.fup, ymin=lower, ymax=upper), alpha = 0.2, fill="orangered") +
      labs(x="combined cognitive score (follow-up)", y=paste(ROI, "follow-up")) +
      ggtitle(ggtitle(heading))
    print(plot.mdl)
    dev.off()
    
  }
  assign(mdl.name,lmdl.PD.V3.CCog)
}

stats.PD.V3.QSM.fup.CCog <- as.data.frame(stats.PD.V3.QSM.fup.CCog)
stats.PD.V3.QSM.fup.CCog$`p-value` <- as.numeric(stats.PD.V3.QSM.fup.CCog$`p-value`)
stats.PD.V3.QSM.fup.CCog$`fdr p-value` <- p.adjust(stats.PD.V3.QSM.fup.CCog$`p-value`, method = "fdr")
stats.PD.V3.QSM.fup.CCog$`p-value` <- format.pval(stats.PD.V3.QSM.fup.CCog$`p-value`,digits = 3)
stats.PD.V3.QSM.fup.CCog$`fdr p-value` <- format.pval(stats.PD.V3.QSM.fup.CCog$`fdr p-value`,digits = 3)
csv.name <- paste("Lmdl.stats.V3.fup.CCog.", data.root , sep = "")
write.csv(stats.PD.V3.QSM.fup.CCog, file = paste(output.dir,filesep,csv.name,".csv",sep=""))
# 

##### PD only linear regression V3 QSM against fup UPDRS-III ----

# create a table to store overall results
stats.PD.V3.QSM.fup.UPDRS.III <- matrix(rep(NaN,times=length(ROI.names)*6), ncol = 6, byrow = T)
colnames(stats.PD.V3.QSM.fup.UPDRS.III) <- c("adjusted R2","multiple R2", "partial R2","beta","p-value", "fdr p-value")
rownames(stats.PD.V3.QSM.fup.UPDRS.III) <- ROI.names

for (ROI in ROI.names) {
  # define a linear regression model for each ROI
  lmdl.PD.V3.UPDRS.III <- lm(PD.V3.QSMdata[,ROI] ~ Age.fup + Sex + UPDRS.III.fup, data = PD.V3.QSMdata)
  relimp <- calc.relimp.lm(lmdl.PD.V3.UPDRS.III)
  if (absQSM) {
    mdl.name <- paste("PD.V3.abs.QSM.lmdl.fup.UPDRS.III.",ROI,sep="")
  } else {
    mdl.name <- paste("PD.V3.signed.QSM.lmdl.fup.UPDRS.III.",ROI,sep="")
  }
  res <- summary(lmdl.PD.V3.UPDRS.III)
  mdl.pval <- pf(res$fstatistic[1],res$fstatistic[2],res$fstatistic[3], lower.tail = FALSE)
  eff.pval <- anova(lmdl.PD.V3.UPDRS.III)$`Pr(>F)`[3]
  eff.slope <- formatC(coef(lmdl.PD.V3.UPDRS.III)[4], format = "e", digits = 2)
  aRsq = res$adj.r.squared
  mRsq <- res$r.squared
  rel.Rsq <- relimp@lmg[3]
  stats.PD.V3.QSM.fup.UPDRS.III[ROI,1] <- round(aRsq,3)
  stats.PD.V3.QSM.fup.UPDRS.III[ROI,2] <- round(mRsq,3)
  stats.PD.V3.QSM.fup.UPDRS.III[ROI,3] <- round(rel.Rsq,3)
  stats.PD.V3.QSM.fup.UPDRS.III[ROI,4] <- eff.slope
  stats.PD.V3.QSM.fup.UPDRS.III[ROI,5] <- eff.pval
  if (eff.pval < 1 && mdl.pval < 1){
    if (eff.pval < mc_thr){
      heading <- paste0("p = ", format.pval(eff.pval, digits = 3), "*mc, aRsq = ", round(aRsq,3), ", lmg = ", round(rel.Rsq,3))
    } else if (eff.pval < a_thr) {
      heading <- paste0("p = ", format.pval(eff.pval, digits = 3), "*, aRsq = ", round(aRsq,3), ", lmg = ", round(rel.Rsq,3))
    } else {
      heading <- paste0("p = ", format.pval(eff.pval, digits = 3), ", aRsq = ", round(aRsq,3), ", lmg = ", round(rel.Rsq,3))
    }
    print(mdl.name)
    print(anova(lmdl.PD.V3.UPDRS.III))
    mdl.effects <- effect(term = "UPDRS.III.fup", mod = lmdl.PD.V3.UPDRS.III, xlevels = 20)
    df.mdl <- as.data.frame(mdl.effects)
    png(filename = paste(output.dir,filesep, "plot_",mdl.name,".png",sep=""))
    plot.mdl <- ggplot() + theme_minimal() +
      geom_point(data = PD.V3.QSMdata, aes(UPDRS.III.fup, PD.V3.QSMdata[,ROI]), color = "darkcyan") +
      geom_line(data = df.mdl, aes(x=UPDRS.III.fup, y=fit), color = "darkcyan") +
      geom_ribbon(data = df.mdl, aes(x=UPDRS.III.fup, ymin=lower, ymax=upper), alpha = 0.2, fill="darkcyan") +
      labs(x="UPDRS-III (follow-up)", y=paste(ROI, "follow-up")) +
      ggtitle(heading)
    print(plot.mdl)
    dev.off()

  }
  assign(mdl.name,lmdl.PD.V3.UPDRS.III)
}

stats.PD.V3.QSM.fup.UPDRS.III <- as.data.frame(stats.PD.V3.QSM.fup.UPDRS.III)
stats.PD.V3.QSM.fup.UPDRS.III$`p-value` <- as.numeric(stats.PD.V3.QSM.fup.UPDRS.III$`p-value`)
stats.PD.V3.QSM.fup.UPDRS.III$`fdr p-value` <- p.adjust(stats.PD.V3.QSM.fup.UPDRS.III$`p-value`, method = "fdr")
stats.PD.V3.QSM.fup.UPDRS.III$`p-value` <- format.pval(stats.PD.V3.QSM.fup.UPDRS.III$`p-value`,digits = 3)
stats.PD.V3.QSM.fup.UPDRS.III$`fdr p-value` <- format.pval(stats.PD.V3.QSM.fup.UPDRS.III$`fdr p-value`,digits = 3)
csv.name <- paste("Lmdl.stats.V3.fup.UPDRS.III.", data.root , sep = "")
write.csv(stats.PD.V3.QSM.fup.UPDRS.III, file = paste(output.dir,filesep,csv.name,".csv",sep=""))
#  
