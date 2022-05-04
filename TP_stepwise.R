#========01 TP STUDY(proteomics-infertile cohort) =============

#   Paper: Novel protein markers of androgen activity in humans: proteomic study of plasma from young chemically castrated men.
# Authors: Aleksander Giwercman1*, K Barbara Sahlin*, Indira Pla, Krzysztof Pawlowski, Carl Fehniger, 
#          Yvonne Lundberg Giwercman, Irene Leijonhufvud, Roger Appelqvist, György Marko-Varga, 
#          Aniel Sanchez, and Johan Malm
#

#  Because the combination of different markers may improve the discriminative power 
#  to diagnose hypogonadism and predict its long term sequelae, proteins selected from 
#  the ROC-AUC analysis were included as predictors in a stepwise regression (method: backward) 
#  to select the best combination of markers that predict the odds of being low T. Bootstrap 
#  resampling with replacement method was applied to assess consistency of predictors selected 
#  with the stepwise regression. 


#====INSTALL PACKAGES=====

# List of packages to install
.packages = c("MASS", "bootStepAIC")

# Install packages if not installed
.inst <- .packages %in% installed.packages()
if(length(.packages[!.inst]) > 0) install.packages(.packages[!.inst])

# Loading packages
lapply(.packages, require, character.only=TRUE)

#======END INSTALL PACKAGES====

#=====Loading dataset================================================================================== 
# Download file "Healthy_model_Signif_biomarkers.xlsx" to your computer and read (upload) it from RStudio 
# to start doing the analysis.

Big.table <- as.data.frame(readxl::read_excel("Healthy_model_Signif_biomarkers.xlsx")) # Read the file data
rownames(Big.table) <- Big.table[,1]
Big.table <- dplyr::select(Big.table, -id)

data1 <- dplyr::select(Big.table, -Patient)

# Run Logistic Regression model with predictors. 

mod1 <- glm(Time_point ~ `4HPPD`+ ALDOB + IGFBP6, 
                          data = data1, 
                          family = "binomial")

summary(mod1)

library(MASS)

mod_step <- stepAIC(mod1, direction = "backward", trace = F)
mod_step


# Bootstraps the Stepwise Algorithm of stepAIC() for Choosing a Model by AIC
library(bootStepAIC)

mod_boot <- boot.stepAIC(mod1, data = data1, B=50)
mod_boot
