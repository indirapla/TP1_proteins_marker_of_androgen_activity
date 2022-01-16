#========01 TP STUDY(proteomics-infertile cohort) =============

#   Paper: Novel protein markers of androgen activity in humans: proteomic study of plasma from young chemically castrated men.
# Authors: Aleksander Giwercman1*, K Barbara Sahlin*, Indira Pla, Krzysztof Pawlowski, Carl Fehniger, 
#          Yvonne Lundberg Giwercman, Irene Leijonhufvud, Roger Appelqvist, György Marko-Varga, 
#          Aniel Sanchez???, and Johan Malm???


#====INSTALL PACKAGES=====

# List of packages to install
.packages = c("pROC","ggplot2","nlme","reshape","dplyr","FactoMineR","verification","ggpubr","gridExtra")

# Install packages if not installed
.inst <- .packages %in% installed.packages()
if(length(.packages[!.inst]) > 0) install.packages(.packages[!.inst])

# Loading packages
lapply(.packages, require, character.only=TRUE)


#===Function to performs ROC analysis and export roc curves per each of the proteins ================
     
#       matrix = dataframe (e.g DATA1) containing variables as columns and samples (observation) as rows.
#     from.col = column number of the first variable that will be plot as ROC curve. 
#       to.col = column number of the last variable that will be plot as ROC curve.
#    group.col = column number of the dichotomic variable (numeric).
#    plots: should the plots be generated? (True or False)
#    files.ID = an id for the files that will be printed (e.g "IR" for those file that come from Insulin resistance analysis)

#    The output is saved in a folder called Results within the work directory (getwd()).
     
plot_roc.fn <- function(matrix,from.col,to.col, group.col, files.ID,plots=F){
 
  dir.create("Results")
  
  matrix <- na.omit(matrix)
  
  # matrix[,group.col] <- as.character(matrix[,group.col])
  # 
  # for (j in 1:nrow(matrix)) {
  #   if (matrix[j,group.col]=="positive"){matrix[j,group.col]<-1} else {matrix[j,group.col]<-0}
  # }
  # 
  # matrix[,group.col] <- as.numeric(matrix[,group.col])
  
  colnames1 <- as.character(colnames(matrix))
  ROC.matrix <- matrix(nrow=ncol(matrix),ncol = 7)
  

  p= 0
  for(i in from.col:to.col) {
    #hist(matrix[,i], main=colnames(matrix)[i])
    
    png(filename = paste("Results/","roc_",files.ID,colnames1[i],".png",sep=""),width = 480, height = 480, units = "px", pointsize = 16,
        bg = "white", res = NA, family = "", restoreConsole = TRUE,type = c("windows", "cairo", "cairo-png"))
    
    
    roc_i <- roc(matrix[,group.col] ~ matrix[,i],data=matrix,percent=T,
                 smooth=F,#levels=c("positive", "negative"),
                 legacy.axes=F, 
                 ci=T, 
                 boot.n=2000, 
                 ci.alpha=0.7, 
                 stratified=T,
                 plot=TRUE, 
                 auc.polygon=T, 
                 max.auc.polygon=T, 
                 grid=TRUE,
                 print.auc=TRUE, 
                 print.thres='best', 
                 print.thres.best.method='y',   # "closest.topleft", "youden"
                 print.thres.adj=c(-0.05, 1.25),#,c(-0.05, 2.0),c(-0.05, 5.0)), c(-0.05, 1.25),c(-0.7, 1.5)
                 print.thres.pattern="Cut-off: %.3f    \n\n\nSpec: %.1f \nSens: %.1f",
                 print.thres.pattern.cex = 1.1,
                 print.auc.cex = 1.1,
                 print.auc.y=7,
                 print.auc.x=83,
                 cex.axis=1.5,
                 cex.lab=1.5,
                 print.thres.pch=16,
                 print.thres.cex=2.0,
                 cex.main=1.5,
                 na.rm=T)
    
    print(roc_i)
    dev.off()
    

    
    ROC.matrix[i,1]<- roc_i$auc
    
    ROC_table <- cbind(roc_i$thresholds, roc_i$specificities, roc_i$sensitivities)
    cut_off <- ROC_table[which.max(ROC_table[, 2] + ROC_table[, 3]), ] 
    
    roc.A <- verification::roc.area(roc_i$original.response,roc_i$original.predictor)
    
    CI <- paste(round(roc_i$ci[1],2),"-",round(roc_i$ci[3],2),sep = "") ## coeficient interval of AUC%
    
    ROC.matrix[i,2] <- CI             ## coeficient interval of AUC%
    ROC.matrix[i,3] <- roc.A$n.total  ## N number of samples
    ROC.matrix[i,4] <- roc.A$p.value  ## p.value
    ROC.matrix[i,5] <- cut_off[1]     ## cut_off = best threshold based on "youden" method
    ROC.matrix[i,6] <- cut_off[2]     ## Sp
    ROC.matrix[i,7] <- cut_off[3]     ## Se
    
    
    roc_i.newlist <- list(AUC = as.matrix(roc_i$auc),
                          sens.= as.matrix(roc_i$sensitivities),
                          Spec.= as.matrix(roc_i$specificities),
                          Thesholds = as.matrix(roc_i$thresholds),
                          predictor = as.matrix(roc_i$predictor),
                          response = as.data.frame(roc_i$response))
    
    write.infile(roc_i.newlist, paste("Results/","roc_",files.ID,colnames1[i],".csv",sep=''), sep = "\t")
    
    p= p+1
  }
  
  row.names(ROC.matrix) <- colnames(matrix)
  colnames(ROC.matrix) <- c("AUC(%)", "CI(%)","N","p.value","cut_off", "Sp", "Se")
  
  ROC.matrix <- ROC.matrix[-1,]
  
  write.csv(ROC.matrix, paste("Results/","ROC_matrix_",files.ID,".csv",sep=""))

}

## ================== ROC per comorbidity ==========================================

#== Low testosterone---by --Total testosterone==========
     
  #== To build the MMA variable ===
     
     glm.ht = glm(Low.testo ~ HPPD + IGFBP6, data = DATA, family="binomial")
     
     MMA <- as.data.frame(predict(glm.ht, type = "response"))   # save the predicted log-odds (of being low testosterone) for each observation
     colnames(MMA) <- "MMA"
     MMA$Kod <- row.names(MMA)
     
     DATA.1 <- plyr::join_all(list(DATA,MMA),by="Kod")
     
     Low.test <- DATA.1[, c("Low.testo", "HPPD","ALDOB","IGFBP6","MMA")]
      
     plot_roc.fn(Low.test,from.col = 2, to.col = ncol(Low.test),group.col = "Low.testo",files.ID="Low.T_",plots = F)
      
     LT.model <- Low.test$Low.testo ~ HPPD+ALDOB+IGFBP6+MMA
     roc_i.list.LT <- roc(LT.model,percent=T, smooth=F, legacy.axes=F,ci=TRUE, data = Low.test)
     roc_i.list.LT

  # Graph Multicurve
   gg.LT<- ggroc(roc_i.list.LT, linetype=1,size = 0.75)+
          theme_bw()+labs()+ggtitle("Low testosterone")+
            theme(legend.title = element_blank())+
          geom_segment(aes(x = 100, xend = 0, y = 0, yend = 100), color="black", linetype="dashed",size = 0.5)
   gg.LT

   #theme(legend.position = "none")   # to dont show the legend
   
  #comparing ROC curves
   
   pROC::roc.test(roc_i.list.LT$HPPD,roc_i.list.LT$Total_testosterone, paired=T, method="delong") # HPPD vs Total Testosterone
   pROC::roc.test(roc_i.list.LT$IGFBP6,roc_i.list.LT$Total_testosterone, paired=T, method="delong") # IGFBP6 vs Total Testosterone
   pROC::roc.test(roc_i.list.LT$ALDOB,roc_i.list.LT$Total_testosterone, paired=T, method="delong") # ALDOB vs Total Testosterone
   pROC::roc.test(roc_i.list.LT$MMA,roc_i.list.LT$Total_testosterone, paired=T, method="delong") # MMA vs Total Testosterone
   
#== Low testosterone 2---by --Total testosterone (excluding borderline low testosterone) ==========
   
   #== To build the MMA (model1) variable ===
   
   glm.ht = glm(LT_nBLin ~ HPPD + IGFBP6, data = DATA, family="binomial")
   
   MMA <- as.data.frame(predict(glm.ht, type = "response"))   # save the predicted log-odds (of being low testosterone) for each observation
   colnames(MMA) <- "MMA"
   MMA$Kod <- row.names(MMA)
   
   DATA.1 <- plyr::join_all(list(DATA,MMA),by="Kod")
   
   Low.test2 <- DATA.1[, c("LT_nBLin", "HPPD","ALDOB","IGFBP6","MMA")]
   
   plot_roc.fn(Low.test2,from.col = 2, to.col = ncol(Low.test2),group.col = "LT_nBLin",files.ID="Low.T2_LT_B_")
   
   LT.model2 <- Low.test2$LT_nBLin ~ HPPD+ALDOB+IGFBP6+MMA
   roc_i.list.LT2 <- roc(LT.model2,percent=T, smooth=F, legacy.axes=F,ci=TRUE, data = Low.test2)
   roc_i.list.LT2
   
   # Graph Multicurve
   gg.LT<- ggroc(roc_i.list.LT2, linetype=1,size = 0.75)+
     theme_bw()+labs()+ggtitle("Low testosterone (LT_B)")+
     theme(legend.title = element_blank())+
     geom_segment(aes(x = 100, xend = 0, y = 0, yend = 100), color="black", linetype="dashed",size = 0.5)
   gg.LT
   
   #theme(legend.position = "none")   # to dont show the legend
   
  #comparing ROC curves
   
   pROC::roc.test(roc_i.list.LT2$HPPD,roc_i.list.LT2$Total_testosterone, paired=T, method="delong") # HPPD vs Total Testosterone
   pROC::roc.test(roc_i.list.LT2$IGFBP6,roc_i.list.LT2$Total_testosterone, paired=T, method="delong") # IGFBP6 vs Total Testosterone
   pROC::roc.test(roc_i.list.LT2$ALDOB,roc_i.list.LT2$Total_testosterone, paired=T, method="delong") # ALDOB vs Total Testosterone
   pROC::roc.test(roc_i.list.LT2$MMA,roc_i.list.LT2$Total_testosterone, paired=T, method="delong") # MMA vs Total Testosterone

#== Low calculated Free testosterone---==========
   
   #== To build the MMA (model1) variable ===
   
   glm.ht = glm(LcFT ~ HPPD + IGFBP6, data = DATA, family="binomial")
   
   MMA <- as.data.frame(predict(glm.ht, type = "response"))   # save the predicted log-odds (of being low testosterone) for each observation
   colnames(MMA) <- "MMA"
   MMA$Kod <- row.names(MMA)
   
   DATA.1 <- plyr::join_all(list(DATA,MMA),by="Kod")
   
   LcFT.test <- DATA.1[, c("LcFT", "HPPD","ALDOB","IGFBP6","MMA")]
   
   plot_roc.fn(LcFT.test,from.col = 2, to.col = ncol(LcFT.test),group.col = "LcFT",files.ID="Low.cFT_")
   

   LcFT.model <- LcFT.test$LcFT ~ HPPD+ALDOB+IGFBP6+MMA
   roc_i.list.LcFT <- roc(LcFT.model,percent=T, smooth=F, legacy.axes=F,ci=TRUE, data = LcFT.test)
   roc_i.list.LcFT
   
   # Graph Multicurve
   gg.LT<- ggroc(roc_i.list.LcFT, linetype=1,size = 0.75)+
     theme_bw()+labs()+ggtitle("Low calc.Free testosterone")+
     theme(legend.title = element_blank())+
     geom_segment(aes(x = 100, xend = 0, y = 0, yend = 100), color="black", linetype="dashed",size = 0.5)
   gg.LT
   
   #theme(legend.position = "none")   # to dont show the legend

 #comparing ROC curves
   
   pROC::roc.test(roc_i.list.LcFT$HPPD,roc_i.list.LcFT$Total_testosterone, paired=T, method="delong") # HPPD vs Total Testosterone
   pROC::roc.test(roc_i.list.LcFT$IGFBP6,roc_i.list.LcFT$Total_testosterone, paired=T, method="delong") # IGFBP6 vs Total Testosterone
   pROC::roc.test(roc_i.list.LcFT$ALDOB,roc_i.list.LcFT$Total_testosterone, paired=T, method="delong") # ALDOB vs Total Testosterone
   pROC::roc.test(roc_i.list.LcFT$MMA,roc_i.list.LcFT$Total_testosterone, paired=T, method="delong") # MMA vs Total Testosterone

#== HOMA --by --Insuline resistance=========================

    #== To build the MMA (model1) variable ===
      
      glm.ht = glm(IR ~ HPPD + IGFBP6, data = DATA, family="binomial")
      
      MMA <- as.data.frame(predict(glm.ht, type = "response"))   # save the predicted log-odds (of being low testosterone) for each observation
      colnames(MMA) <- "MMA"
      MMA$Kod <- row.names(MMA)
      
      DATA.1 <- plyr::join_all(list(DATA,MMA),by="Kod")

      HOMA.ir <- DATA.1[, c("IR", "HPPD","ALDOB","IGFBP6","Total_testosterone","MMA")]
    
      plot_roc.fn(HOMA.ir,from.col = 2, to.col = ncol(HOMA.ir),group.col = "IR",files.ID="IR_",plots = F)
    
      HOMA.model <- HOMA.ir$IR ~ HPPD + ALDOB + IGFBP6 + Total_testosterone + MMA
      roc_i.list.HOMA <- roc(HOMA.model,percent=T, smooth=F, legacy.axes=F,ci=TRUE, data = HOMA.ir)
      roc_i.list.HOMA

  # Graph Multicurve
    gg.IR<-ggroc(roc_i.list.HOMA, linetype=1,size = 0.75)+
            theme_bw()+labs()+ggtitle("Insulin resistance")+theme(legend.position = "none")+
            geom_segment(aes(x = 100, xend = 0, y = 0, yend = 100), color="black", linetype="dashed",size = 0.5)
    gg.IR 
    
  #comparing ROC curves
    
    pROC::roc.test(roc_i.list.HOMA$HPPD,roc_i.list.HOMA$Total_testosterone, paired=T, method="delong") # HPPD vs Total Testosterone
    pROC::roc.test(roc_i.list.HOMA$IGFBP6,roc_i.list.HOMA$Total_testosterone, paired=T, method="delong") # IGFBP6 vs Total Testosterone
    pROC::roc.test(roc_i.list.HOMA$ALDOB,roc_i.list.HOMA$Total_testosterone, paired=T, method="delong") # ALDOB vs Total Testosterone
    pROC::roc.test(roc_i.list.HOMA$MMA,roc_i.list.HOMA$Total_testosterone, paired=T, method="delong") # MMA vs Total Testosterone
    
#== Diabetes --==========================================================
    
  #== To build the MMA (model1) variable ===
    DATA$DM <- as.numeric(DATA$DM)
    glm.ht = glm(DM ~ HPPD + IGFBP6, data = DATA, family="binomial")
    
    MMA.DM <- as.data.frame(predict(glm.ht, type = "response"))   # save the predicted log-odds (of being low testosterone) for each observation
    colnames(MMA.DM) <- "MMA.DM"
    MMA.DM$Kod <- row.names(MMA.DM)
    
    DATA.1 <- plyr::join_all(list(DATA,MMA.DM),by="Kod")
    
    #=

    DM <- DATA.1[, c("DM", "HPPD","ALDOB","IGFBP6","Total_testosterone","MMA")]
    
    plot_roc.fn(DM,from.col = 2, to.col = ncol(DM),group.col = "DM",files.ID="DM_",plots = F)
    
    DM.model <- DM ~ HPPD + ALDOB + IGFBP6 + Total_testosterone + MMA
    roc_i.list.DM <- roc(DM.model,percent=T, smooth=F, legacy.axes=F,ci=TRUE, data = DM)
    roc_i.list.DM
    
    # Graph Multicurve
    gg.DM<-ggroc(roc_i.list.DM, linetype=1,size = 0.75)+
      theme_bw()+labs()+ggtitle("DM")+theme(legend.position = "none")+
      geom_segment(aes(x = 100, xend = 0, y = 0, yend = 100), color="black", linetype="dashed",size = 0.5)
    gg.DM
    
 #comparing ROC curves
    
    pROC::roc.test(roc_i.list.DM$HPPD,roc_i.list.DM$Total_testosterone, paired=T, method="delong") # HPPD vs Total Testosterone
    pROC::roc.test(roc_i.list.DM$IGFBP6,roc_i.list.DM$Total_testosterone, paired=T, method="delong") # IGFBP6 vs Total Testosterone
    pROC::roc.test(roc_i.list.DM$ALDOB,roc_i.list.DM$Total_testosterone, paired=T, method="delong") # ALDOB vs Total Testosterone
    pROC::roc.test(roc_i.list.DM$MMA,roc_i.list.DM$Total_testosterone, paired=T, method="delong") # MMA vs Total Testosterone
    
#== Cardivacular disease ---ratio ApoB/ApoA1 ===============

      #== To build the MMA (model1) variable ===
      
      glm.ht = glm(CVD ~ HPPD + IGFBP6, data = DATA, family="binomial")
      
      MMA <- as.data.frame(predict(glm.ht, type = "response"))   # save the predicted log-odds (of being low testosterone) for each observation
      colnames(MMA) <- "MMA"
      MMA$Kod <- row.names(MMA)
      
      DATA.1 <- plyr::join_all(list(DATA,MMA),by="Kod")
      
      Apo_B_A1_ratio <- DATA.1[, c("CVD", "HPPD","ALDOB","IGFBP6","Total_testosterone","MMA")]
      
      plot_roc.fn(Apo_B_A1_ratio,from.col = 2, to.col = ncol(Apo_B_A1_ratio),
                  group.col = "CVD", files.ID="CVD_",plots = F)
      
      Apo_BA1ratio.model <- Apo_B_A1_ratio$CVD ~ HPPD+ALDOB+IGFBP6+Total_testosterone+MMA
      roc_i.list.apo <- roc(Apo_BA1ratio.model,percent=T, smooth=F, legacy.axes=F,ci=TRUE, data = Apo_B_A1_ratio)
      roc_i.list.apo

      # Graph Multicurve
      gg.cvd <- ggroc(roc_i.list.apo, linetype=1,size = 0.75)+
                     theme_bw()+labs()+ ggtitle("CVR")+theme(legend.position = "none")+
                     geom_segment(aes(x = 100, xend = 0, y = 0, yend = 100), color="black", linetype="dashed",size = 0.5)
      gg.cvd

  #comparing ROC curves
      
      pROC::roc.test(roc_i.list.apo$HPPD,roc_i.list.apo$Total_testosterone, paired=T, method="delong") # HPPD vs Total Testosterone
      pROC::roc.test(roc_i.list.apo$IGFBP6,roc_i.list.apo$Total_testosterone, paired=T, method="delong") # IGFBP6 vs Total Testosterone
      pROC::roc.test(roc_i.list.apo$ALDOB,roc_i.list.apo$Total_testosterone, paired=T, method="delong") # ALDOB vs Total Testosterone
      pROC::roc.test(roc_i.list.apo$MMA,roc_i.list.apo$Total_testosterone, paired=T, method="delong") # MMA vs Total Testosterone
      
#=== Metabolic Syndrome =====================================

    #== To build the MMA (model1) variable ===
      
      glm.ht = glm(MetS_AHA_NHLBI.1 ~ HPPD + IGFBP6, data = DATA, family="binomial")
            
      MMA <- as.data.frame(predict(glm.ht, type = "response"))   # save the predicted log-odds (of being low testosterone) for each observation
      colnames(MMA) <- "MMA"
      MMA$Kod <- row.names(MMA)
      
      DATA.1 <- plyr::join_all(list(DATA,MMA),by="Kod")
            
      Met.S <- DATA.1[, c("MetS_AHA_NHLBI.1", "HPPD","ALDOB","IGFBP6","Total_testosterone","MMA")]
      
      plot_roc.fn(Met.S,from.col = 2, to.col = ncol(Met.S),group.col = "MetS_AHA_NHLBI.1",files.ID="Met.S_", plots = F)
      
      Met.S.model <- Met.S$MetS_AHA_NHLBI.1 ~ HPPD+ALDOB+IGFBP6+Total_testosterone+MMA
      roc_i.list.Met.S <- roc(Met.S.model,percent=T, smooth=F, legacy.axes=F,ci=TRUE, data = Met.S)
      roc_i.list.Met.S

      # Graph Multicurve
      gg.MS<- ggroc(roc_i.list.Met.S, linetype=1,size = 0.75)+
              theme_bw()+labs()+ ggtitle("MetS")+theme(legend.position = "none")+
              geom_segment(aes(x = 100, xend = 0, y = 0, yend = 100), color="black", linetype="dashed",size = 0.5)
      gg.MS

  #comparing ROC curves
      
      pROC::roc.test(roc_i.list.Met.S$HPPD,roc_i.list.Met.S$Total_testosterone, paired=T, method="delong") # HPPD vs Total Testosterone
      pROC::roc.test(roc_i.list.Met.S$IGFBP6,roc_i.list.Met.S$Total_testosterone, paired=T, method="delong") # IGFBP6 vs Total Testosterone
      pROC::roc.test(roc_i.list.Met.S$ALDOB,roc_i.list.Met.S$Total_testosterone, paired=T, method="delong") # ALDOB vs Total Testosterone
      pROC::roc.test(roc_i.list.Met.S$MMA,roc_i.list.Met.S$Total_testosterone, paired=T, method="delong") # MMA vs Total Testosterone
      
#=== Bone density-----by-------DEXA_lumbar_Z_score============

    #== To build the MMA (model1) variable ===
      
      glm.ht = glm(LBD ~ HPPD + IGFBP6, data = DATA, family="binomial")
      
      MMA <- as.data.frame(predict(glm.ht, type = "response"))   # save the predicted log-odds (of being low testosterone) for each observation
      colnames(MMA) <- "MMA"
      MMA$Kod <- row.names(MMA)
      
      DATA.1 <- plyr::join_all(list(DATA,MMA),by="Kod")

      DEXA_lumbar_Z_score <- DATA.1[, c("LBD", "HPPD","ALDOB","IGFBP6","Total_testosterone","MMA")]
      
      plot_roc.fn(DEXA_lumbar_Z_score,from.col = 2, to.col = ncol(DEXA_lumbar_Z_score),
                  group.col = "LBD", files.ID = "LBD_",plots = F)
      
      Z_score.model <- DEXA_lumbar_Z_score$LBD ~ HPPD+ALDOB+IGFBP6+Total_testosterone+MMA
      roc_i.list.z.score <- roc(Z_score.model,percent=T, smooth=F, legacy.axes=F,ci=TRUE, data = DEXA_lumbar_Z_score)
      roc_i.list.z.score

      # Graph Multicurve
      gg.LBD<- ggroc(roc_i.list.z.score, linetype=1,size = 0.75)+
                     theme_bw()+labs()+ ggtitle("Low bone density")+theme(legend.position = "none")+
                     geom_segment(aes(x = 100, xend = 0, y = 0, yend = 100), color="black", linetype="dashed",size = 0.5)
      gg.LBD
      
  #comparing ROC curves
      
      pROC::roc.test(roc_i.list.z.score$HPPD,roc_i.list.z.score$Total_testosterone, paired=T, method="delong") # HPPD vs Total Testosterone
      pROC::roc.test(roc_i.list.z.score$IGFBP6,roc_i.list.z.score$Total_testosterone, paired=T, method="delong") # IGFBP6 vs Total Testosterone
      pROC::roc.test(roc_i.list.z.score$ALDOB,roc_i.list.z.score$Total_testosterone, paired=T, method="delong") # ALDOB vs Total Testosterone
      pROC::roc.test(roc_i.list.z.score$MMA,roc_i.list.z.score$Total_testosterone, paired=T, method="delong") # MMA vs Total Testosterone
      
