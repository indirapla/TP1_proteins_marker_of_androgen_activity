# TP1_proteins_marker_of_androgen_activity

This repository contains R code used to perform the data analysis descrived in the manuscript:
 
### Novel potentially clinically valuable protein markers of androgen activity in humans 

#### File: one-way repeated measurement ANOVA.R:

Differentially expressed proteins were determined by doing one-way repeated measures ANOVA (R function: ezANOVA{ez}) to reveal overall differences between conditions (three time points) and differences between individual conditions were detected by performing a post-hoc test based on pairwise t-test (two-tails and paired) (R function: pairwise.t.test{stats}). The ‘pairwise.t.test’ function utilized the proteins with significant overall changes (ANOVA p-value < 0.05) to perform pairwise comparisons between conditions while corrected for multiple pairwise testing. Proteins with adjusted p-values (‘fdr’ method) < 0.05 following the pairwise t-test were considered significant

#### File: ROC analysis_TPstudy.R

Receiver operating characteristic (ROC) analysis to select proteins capable to discriminate between normal and low Testosterone. (Healthy human model)

#### File: TP_stepwise.R

A stepwise regression (method: backward) to select the best combination of markers that predict the odds of being low Testosterone. Bootstrap resampling with replacement method was applied to assess consistency of predictors selected with the stepwise regression.

#### File: ROC analysis_TPstudy_infertile cohort.R

ROC analysis to discriminate patients with MetS, IR, CVRLP, DM or LBD within the cohort infertile men using the expression level of the candidate biomarkers and Testosterone hormone
