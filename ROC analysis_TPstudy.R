#========01 TP STUDY(proteomics) =============

#   Paper: Novel potentially clinically valuable protein markers of androgen activity in humans 
# Authors: Aleksander Giwercman1*, K Barbara Sahlin*, Indira Pla, Krzysztof Pawlowski, Carl Fehniger, 
#          Yvonne Lundberg Giwercman, Irene Leijonhufvud, Roger Appelqvist, Gy�rgy Marko-Varga, 
#          Aniel Sanchez???, and Johan Malm???

#====INSTALL PACKAGES=====

# List of packages to install
.packages = c("pROC","ggplot2","nlme","reshape","dplyr","FactoMineR")

# Install packages if not installed
.inst <- .packages %in% installed.packages()
if(length(.packages[!.inst]) > 0) install.packages(.packages[!.inst])

# Loading packages
lapply(.packages, require, character.only=TRUE)

#========= Receiver operating characteristic (ROC) analysis==================================

#data("TO roc plots3_proteins.xlsx")


#===Function plot_roc.fn=============================

# This function performs ROC analysis and export roc curves per each of the proteins.  

#       matrix: dataframe containing variables as columns and samples (observation) as rows.
#     from.col: number of the column that contains the first variable to plot as ROC curve. 
#       to.col: number of the column that contains the last variable to plot as ROC curve.
#    group.col: number of the column that contains the dichotomic variable.
#        plots: should the plots be generated? (True or False)

plot_roc.fn <- function(matrix,from.col,to.col, group.col, plots=F){
  
  colnames1 <- as.character(colnames(matrix))
  AUC.matrix <- matrix(nrow=ncol(matrix),ncol = 3)
  
  colnames(AUC.matrix) <- c("AUC.CI","LL (95% CI)","UL (95% CI)")
  
  for(i in from.col:to.col) {
    #hist(matrix[,i], main=colnames(matrix)[i])
    
    png(filename = paste("roc_",colnames1[i],".png",sep=""),width = 480, height = 480, units = "px", pointsize = 16,
        bg = "white", res = NA, family = "", restoreConsole = TRUE,type = c("windows", "cairo", "cairo-png"))

    
    roc_i <- roc(matrix[,group.col] ~ matrix[,i],data=matrix,percent=T,
                 smooth=F,
                 legacy.axes=F,
                 ci=F, 
                 boot.n=2000, 
                 ci.alpha=0.7, 
                 stratified=T,
                 plot=plots, 
                 auc.polygon=T, 
                 max.auc.polygon=T, 
                 grid=TRUE,
                 print.auc=TRUE, 
                 print.thres='best', 
                 print.thres.best.method='y',   # "closest.topleft", "youden"
                 print.thres.adj=c(-0.7, 1.5),#,c(-0.05, 2.0),c(-0.05, 5.0)),
                 print.thres.pattern="Cut-off: %.2f    \n\n\nSpec: %.1f \nSens: %.1f",
                 print.thres.pattern.cex = 1.1,
                 print.auc.cex = 1.1,
                 print.auc.y=7,
                 print.auc.x=83,
                 cex.axis=1.5,
                 cex.lab=1.5,
                 print.thres.pch=16,
                 print.thres.cex=2.0,
                 cex.main=1.5)

    
    print(roc_i)
    
    dev.off()
    
    
    # Compute the confidence interval of the AUC (bootstrap method)
    ci.auc1 <- ci.auc(roc_i, 
                      conf.level=0.95, 
                      boot.n = 10000,
                      method = "bootstrap")  # Default: 2000.
    

    
    AUC.matrix[i,1] <- roc_i$auc
    AUC.matrix[i,2] <- ci.auc1[1]
    AUC.matrix[i,3] <- ci.auc1[3]
    
    # AUC.vector <- c(paste("auc = ",round(roc_i$auc,3), sep = ""),
    #                 paste("LL(95% CI) = ",round(ci.auc1[1],3), sep = ""),
    #                 paste("UP(95% CI) = ",round(ci.auc1[3],3), sep = ""))
    
    roc_i.newlist <- list(AUC = as.matrix(roc_i$auc),
                          sens.= as.matrix(roc_i$sensitivities),
                          Spec.= as.matrix(roc_i$specificities),
                          Thesholds = as.matrix(roc_i$thresholds),
                          predictor = as.matrix(roc_i$predictor),
                          response = as.data.frame(roc_i$response),
                          plots=F)
    
    write.infile(roc_i.newlist, paste("roc_",colnames1[i],".txt",sep=''), sep = "\t")

  }
  row.names(AUC.matrix) <- colnames(matrix)
  
  write.csv(AUC.matrix, "AUC.matrix.csv")
}


