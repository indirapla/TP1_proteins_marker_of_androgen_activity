#========01 TP STUDY(proteomics) =============

#   Paper: Novel protein markers of androgen activity in humans: proteomic study of plasma from young chemically castrated men.
# Authors: Aleksander Giwercman1*, K Barbara Sahlin*, Indira Pla, Krzysztof Pawlowski, Carl Fehniger, 
#          Yvonne Lundberg Giwercman, Irene Leijonhufvud, Roger Appelqvist, György Marko-Varga, 
#          Aniel Sanchez†, and Johan Malm†



#====Install packages=====

# List of packages to install
.packages = c("devtools","ggplot2", "reshape2", 
              "tidyverse","dplyr","ez")

# Install packages if not installed
.inst <- .packages %in% installed.packages()
if(length(.packages[!.inst]) > 0) {install.packages(.packages[!.inst])}

# Loading packages
lapply(.packages, require, character.only=TRUE)


#== one-way repeated measurements ANOVA =======


#-------Function ANOVA_paired-------- 

# This function performe one-way repeated measurements ANOVA per protein and then 
# performs a post-hoc test based on pairwise t-test (two-tails and paired) 

#      m: dataframe contating protein expresion per sample.
#   n.gr: number of groups to compare
# groups: dataframe containing 3 columns with sample infomation (column1= Patients, column2 = cond, column3 = sample)


ANOVA_paired <- function(m, groups, n.gr){
  

  k = 2  # binary comparison
  comb <- factorial(n.gr)/(factorial(k)*factorial(n.gr-k))  # number of combinations
  
  p.value <- matrix(nrow=nrow(m),ncol = 2)
  protein.name <- vector(mode = "character")
  mean.new <- matrix(nrow=nrow(m),ncol = comb)
  Fold.change <- matrix(nrow=nrow(m),ncol = comb)
  t.test.posthoc <- matrix(nrow=nrow(m),ncol = comb+1)
  t.test.posthoc.adj <- matrix(nrow=nrow(m),ncol = comb+1)
  num.samples <- vector(mode = "numeric")
  
  
  
  for (i in 1:nrow(m)) {
    data <- cbind(as.data.frame(t(m[i,])),groups)
    prot.name <- colnames(data[1])
    
    data1 <- cbind(subset(data,cond == "A"),subset(data,cond == "B"),subset(data,cond == "C"))
    
    data1.complete <- data1[complete.cases(data1),]
    
    num.samples[i] <- nrow(data1.complete)
    
    data2 <- rbind(data1.complete[1:4],data1.complete[5:8],data1.complete[9:12])
    row.names(data2) <- data2$sample
    
    colnames(data2) <- c("protein", colnames(data2[2:length(colnames(data2))]))
   
    ez= ezANOVA(data=data2, dv=protein, wid=Patient, within = cond)      #ez= ezANOVA(data=data, dv=P63104, wid=Patient, within = cond) #Anova 
    p.value[i,1] <- ez$ANOVA$p
    p.value[i,2] <- ez$ANOVA$`p<.05`
    
    
    protein.name[i] <- prot.name
    
    
    # ANOVA post hoc
    t.test.pv<- with(data2, pairwise.t.test(protein, cond, p.adjust.method="none", paired =T))
    t.test.pv <- reshape2::melt(t.test.pv$p.value)
    t.test.posthoc[i,] <- t.test.pv[,3]
    colnames(t.test.posthoc) <- c(paste("t.test.paired","(",t.test.pv[1,1],"-",t.test.pv[1,2],")"),
                                  paste("t.test.paired","(",t.test.pv[2,1],"-",t.test.pv[2,2],")"),
                                  paste("t.test.paired","(",t.test.pv[3,1],"-",t.test.pv[3,2],")"),
                                  paste("t.test.paired","(",t.test.pv[4,1],"-",t.test.pv[4,2],")"))
    # ANOVA post hoc-adjusted
    t.test.pv<- with(data2, pairwise.t.test(protein, cond, p.adjust.method="fdr", paired =T))
    t.test.pv <- reshape2::melt(t.test.pv$p.value)
    t.test.posthoc.adj[i,] <- t.test.pv[,3]
    colnames(t.test.posthoc.adj) <- c(paste("t.test.paired.adj","(",t.test.pv[1,1],"-",t.test.pv[1,2],")"),
                                  paste("t.test.paired.adj","(",t.test.pv[2,1],"-",t.test.pv[2,2],")"),
                                  paste("t.test.paired.adj","(",t.test.pv[3,1],"-",t.test.pv[3,2],")"),
                                  paste("t.test.paired.adj","(",t.test.pv[4,1],"-",t.test.pv[4,2],")"))
    
    mean.new_A <- mean(subset(data2,cond=="A")$protein)
    mean.new_B <- mean(subset(data2,cond=="B")$protein)
    mean.new_C <- mean(subset(data2,cond=="C")$protein)
    
    mean.new[i,1] <-mean.new_A
    mean.new[i,2] <-mean.new_B
    mean.new[i,3] <-mean.new_C
    
    fold.change.AB <- mean.new_B-mean.new_A
    fold.change.BC <- mean.new_C-mean.new_B
    fold.change.AC <- mean.new_C-mean.new_A
    
    Fold.change[i,1] <- fold.change.AB 
    Fold.change[i,2] <- fold.change.BC 
    Fold.change[i,3] <- fold.change.AC 

  }
  colnames(p.value) <- c("ANOVA p.value", "ANOVA p<.05")
  colnames(Fold.change) <- c("Log2.FC(B-A)", "Log2.FC(C-B)","Log2.FC(C-A)")
  colnames(mean.new) <- c("mean.new.A", "mean.new.B","mean.new.C")
  
  final.data <- cbind(protein.name,mean.new,p.value,t.test.posthoc,t.test.posthoc.adj,Fold.change, num.samples)
  
  return(final.data)
}




