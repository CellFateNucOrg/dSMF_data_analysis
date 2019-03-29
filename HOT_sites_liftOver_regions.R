library(rtracklayer)
library(liftOver)

setwd("C:/Users/pmeister/AppData/Local/Packages/CanonicalGroupLimited.UbuntuonWindows_79rhkp1fndgsc/LocalState/rootfs/home/pmeister/dSMF_data_analysis")


HOT_sites <- import.bed("worm_ce10_HOT_core_midpoints.bed")
chain <- import ("liftOver_chains/ce10ToCe11.over.chain")
HOT_sites_ce11 <- unlist(liftOver(HOT_sites,chain))
HOT_sites_ce11 <- resize (HOT_sites_ce11, fix="center", width=2000)
saveRDS(HOT_sites_ce11,"HOT_sites_ce11_2kb.RDS")
