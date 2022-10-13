

#====INSTALL PACKAGES=====

# List of packages
.packages = c("beeswarm","ggplot2",'dplyr','tidyverse')     

# Installing packages
.inst <- .packages %in% installed.packages()
if(length(.packages[!.inst]) > 0) install.packages(.packages[!.inst])

# loading packages
lapply(.packages, require, character.only=TRUE)


#==========Beeswarm - boxplot=========

# m = protein expression matrix
# labels.beeswarm <- c("A","B","C")
# labels.boxplot <- c("","","")
# xlab = ''
# ylab = 'Intensity (normalized)'
# main = as.character(ylabels[,2])
# method = "swarm"                  #  method = c("swarm", "center", "hex", "square")
# maindata_first.colum = 3          # numero de la columna donde empiezan las variables a graficar
# class.colum = 2                   # numero de la columna que tiene los grupos



beeswarm.function <- function(m,maindata_first.colum, class.colum, labels.beeswarm = c("A","B","C"),
                              labels.boxplot=c("","",""), xlab = '', ylab = 'Intensity (normalized)', 
                              method = "swarm",format.plot=c('svg','png')){
  
  m1 <-  as.data.frame(m[,maindata_first.colum:ncol(m)])
  rownames(m1) <- rownames(m)
  
  groups2 <- as.data.frame(m[,class.colum])
  groups2$sample <- rownames(m)
  
  
  dir.create("output/plot_box")
  #colnames1 <- as.character(colnames(m1))
  #-----
  
  for (i in 1:ncol(m1)) {
    
    #vaiable.name <- colnames1[i]
    
    prot.data <- as.data.frame(m1[,i])
    colnames(prot.data) <- 'protein'
    rownames(prot.data) <- rownames(m1)
    
    prot.data$sample <- rownames(prot.data)
    
    data <- plyr::join_all(list(groups2,prot.data),by="sample")
    row.names(data) <- data$sample
    data <- data %>% dplyr::select(-sample)
    data <- na.omit(data)
    data$protein <- as.numeric(data$protein)
    
    
    
    if(format.plot=='svg'){
      
      svg(paste("output/plot_box/","beeswarm_",colnames(m1)[i],"_",".svg",sep=""))
      
      beeswarm(data[,2]~ data[,1], data = data, method = method,         
               pch = 16, cex=0.6, corral = "gutter", labels = labels.beeswarm,
               xlab = xlab, ylab = ylab, main = colnames(m1)[i], col = c("yellowgreen", "red","blue"))
      
      graph <- boxplot(data[,2]~ data[,1], data = data, add = T,
                       names = labels.boxplot, col="#0000ff22")
      
      print(graph)
      dev.off()
      
    } else {
      
      png(filename = paste("output/plot_box/","beeswarm_",colnames(m1)[i],"_",".png",sep=""),width = 480, height = 480, units = "px", pointsize = 15,
          bg = "white", res = NA, family = "", restoreConsole = TRUE,type = c("windows", "cairo", "cairo-png"))
      
      beeswarm(data[,2]~ data[,1], data = data, method = method,         
               pch = 16, cex=0.6, corral = "gutter", labels = labels.beeswarm,
               xlab = xlab, ylab = ylab, main = colnames(m1)[i], col = c("yellowgreen", "red","blue"))
      
      graph <- boxplot(data[,2]~ data[,1], data = data, add = T,
                       names = labels.boxplot, col="#0000ff22")
      
      print(graph)
      dev.off()
    }
  }
}

beeswarm.function(m=Big_Tabla,maindata_first.colum=3, class.colum=2, 
                  labels.beeswarm = c("A","B","C"),
                  labels.boxplot=c("","",""), 
                  xlab = '', 
                  ylab = 'Intensity (normalized)', 
                  method = "swarm",format.plot='svg')





#===========all together

plot_all <- function(m,class.colum, 
                     from.col,to.col,
                     method = "swarm",
                     labels.boxplot = c("","",""),
                     labels.beeswarm = c("A","B","C"),
                     xlab,ylab){
  
  m1 <-  as.data.frame(m[,from.col:to.col])
  rownames(m1) <- rownames(m)
  
  groups2 <- as.data.frame(m[,class.colum])
  groups2$sample <- rownames(m)
  
  #-----
  
  for (i in 1:ncol(m1)) {
    
    
    prot.data <- as.data.frame(m1[,i])
    colnames(prot.data) <- 'protein'
    rownames(prot.data) <- rownames(m1)
    
    prot.data$sample <- rownames(prot.data)
    
    data <- plyr::join_all(list(groups2,prot.data),by="sample")
    row.names(data) <- data$sample
    data <- data %>% dplyr::select(-sample)
    data <- na.omit(data)
    data$protein <- as.numeric(data$protein)
    
    
    beeswarm(data[,2]~ data[,1], data = data, method = method,         
             pch = 16, cex=0.6, corral = "gutter", labels = labels.beeswarm,
             xlab = xlab, ylab = ylab, main = colnames(m1)[i], col = c("yellowgreen", "red","blue"))
    
    graph <- boxplot(data[,2]~ data[,1], data = data, add = T,
                     names = labels.boxplot, col="#0000ff22")
    
    graph
  }
}

par(mfrow = c(3,5), mar=c(2,4,2,1))   #con esto dividimos el ?rea para los graficos en 3 columnas y 1 fila

plot_all(m = Big_Tabla, class.colum = 2, 
         from.col = 18,to.col = 30,
         method = "swarm",
         labels.beeswarm = c("A","B","C"),
         labels.boxplot = c("","",""),
         xlab ='',ylab = 'Intensity (normalized)')








