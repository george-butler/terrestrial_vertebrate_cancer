library(adephylo)
library(pbmcapply)

#this script is designed to work with a posterior distribution of rate scaled tree from Cooney & Thomas "Heterogeneous relationships between rates of speciation and body size evolution across vertebrate clades"

path_measure<-function(t,tip_ids){
  opt<-lapply(c("patristic","nNodes"), function(e) distRoot(t,tips = tip_ids, method=e))
  names(opt)<- c("path_length","no_nodes")
  opt<-as.data.frame(opt)
  opt$species<-row.names(opt)
  return(opt)
}

#set number of cores
no_cores=50

cancer_df<-read.csv("~/data/species_data.csv")

v<-c("amphibians","birds","mammals","squamates")
for (i in v){
  #point to the location of the Cooney and Thomas posterior trees
  data<-readRDS(paste0("~/",i,"/traits.txt.Output.trees.rds"))
  data<-data[c(2:length(data))]
  tip_labels<-data[[1]]$tip.label
  target_tips<-tip_labels[tip_labels %in% cancer_df$Species]
  summary_stats<-pbmclapply(data,path_measure,tip_ids=target_tips, mc.cores=no_cores,ignore.interactive = T)
  summary_stats<-data.table::rbindlist(summary_stats, idcol = TRUE)
  write.table(summary_stats,paste0("~/data/",i,"/full_pl_data.txt"),row.names = FALSE,sep="\t",quote = FALSE)
  
  summary_stats<-subset(summary_stats, select = -c(.id))
  median_summary_stats<-aggregate(summary_stats[, 1:2], list(summary_stats$species), median)
  write.table(median_summary_stats,paste0("~/data/",i,"/median_pl_data.txt"),row.names = FALSE,sep="\t",quote = FALSE)
}

