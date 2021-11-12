# TP1_proteins_marker_of_androgen_activity

This repository contains R code used to perform the data analysis descrived in the manuscript:
 
### Novel potentially clinically valuable protein markers of androgen activity in humans 

#### one-way repeated measurement ANOVA.R:

Differentially expressed proteins were determined by doing one-way repeated measures ANOVA (R function: ezANOVA{ez}) to reveal overall differences between conditions (three time points) and differences between individual conditions were detected by performing a post-hoc test based on pairwise t-test (two-tails and paired) (R function: pairwise.t.test{stats}). The ‘pairwise.t.test’ function utilized the proteins with significant overall changes (ANOVA p-value < 0.05) to perform pairwise comparisons between conditions while corrected for multiple pairwise testing. Proteins with adjusted p-values (‘fdr’ method) < 0.05 following the pairwise t-test were considered significant
