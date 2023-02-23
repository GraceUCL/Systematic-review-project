# UCL_project
library(meta)
all<-read.csv(file="NEW2021demenita.csv",encoding="UTF-8")
source("plot_forest_sesp20211227utf8.R")
sensitivity_all <- metaprop(all$TP, all$TP+all$FN, fixed=FALSE, random=TRUE, sm="PLOGIT", method.ci="CP", studlab=all$Author, subgroup =all$Test)
specificity_all <- metaprop(all$TN, all$TN+all$FP,fixed=FALSE, random=TRUE, sm="PLOGIT", method.ci="CP", studlab=all$Author, subgroup=all$Test)
plotsesp(sensitivity_all,specificity_all,digit=3,savepdf="dementia0117-mix.pdf",overall=FALSE,headline=TRUE,noRawdata=TRUE,gtestgray=40)
