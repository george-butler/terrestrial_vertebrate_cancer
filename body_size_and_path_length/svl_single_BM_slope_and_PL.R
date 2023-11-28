library(MCMCglmm)
library(ggplot2)
library(ape)
library(pbmcapply)
##########################################################################
#read in the tree 
tree <- read.nexus(file="~/data/terrestrial_vertebrate_tree.nexus.trees")
tree <- ladderize(tree, right=TRUE)
tree<- di2multi(tree,0.000001)
ultra <- is.ultrametric(tree)
ultra
treeAinv <- inverseA(tree)$Ainv

#read in the data 
data<-read.csv('~/data/species_data.csv')
data<-data[data$Class %in% c("Amphibia","Reptilia"),]
#########################################################################
data$Class<-as.factor(data$Class)
class_levels<-levels(data$Class)
number_iterations<-1e7
sample_interval<-1e3
burnin_iterations<-9e6

number_of_samples<-((number_iterations - burnin_iterations) / sample_interval)
#########################################################################
data$reptile = ifelse(data$Class == "Reptilia", 1, 0)
data$log_bm_rep = data$reptile*data$log_bm
data$log_pl_rep = data$reptile*data$log_pl
#########################################################################
#set up the model priors 
pr1<- list(R=list(V=diag(2), nu=2, fix=2),B=list(mu=(rep(0,12)),V=diag(12)*1e8),
           G=list(G1=list(V=diag(2), nu=2, alpha.mu=cbind(0,0), alpha.V=diag(2)*25^2)))
#########################################################################
species_list<-as.list(unique(data$Species))

#cross validation function
cross_validate_seperate_slopes<-function(taxa_remove){
  data<-data[data$Species != taxa_remove,]
  mod1<-MCMCglmm(cbind(neo,mal) ~ trait+trait:reptile+trait:log_rec + trait:log_bm + trait:log_pl + trait:log_pl_rep, family=cbind("poisson","poisson"), pr=TRUE,prior = pr1, data=data, pl=TRUE, 
                 rcov = ~ us(trait):units, random = ~us(trait):Species,nitt=number_iterations,thin=sample_interval,burnin=burnin_iterations, verbose=FALSE, DIC=TRUE, ginverse=list(Species=treeAinv))
  saveRDS(mod1,file=paste0("~/results/n_1_single_svl_seperate_pl_",taxa_remove,".rds"))
}

#set the number of cores
no_cores<-50
cross_validation_results<-data.table::rbindlist(pbmclapply(species_list,cross_validate_seperate_slopes, mc.cores=no_cores,ignore.interactive = T))

neoplasia_outlier_results<-array(data=NA,c(0,2))
malignant_outlier_results<-array(data=NA,c(0,2))
for (i in 1:length(species_list)){
  mod1<-readRDS(paste0("~/results/n_1_single_svl_seperate_pl_",species_list[i],".rds"))
  pred<-as.data.frame(predict.MCMCglmm(mod1,data,interval="prediction",type="terms",marginal = NULL))
  
  pred_df<-cbind(pred[1:nrow(data),],pred[(nrow(data)+1):(2*nrow(data)),])
  pred_df<-pred_df[,c(1,4)]
  colnames(pred_df)<-c("neo_fit","mal_fit")
  
  output<-cbind(data,pred_df)
  output$neo_res<-output$neo - output$neo_fit
  neo_res_val<-output[output$Species == species_list[i],]$neo_res/sd(output[output$Species != species_list[i],]$neo_res)
  tmp<-cbind(species_list[i],neo_res_val)
  neoplasia_outlier_results<-rbind(neoplasia_outlier_results,tmp)
  
  output$mal_res<-output$mal - output$mal_fit
  mal_res_val<-output[output$Species == species_list[i],]$mal_res/sd(output[output$Species != species_list[i],]$mal_res)
  tmp<-cbind(species_list[i],mal_res_val)
  malignant_outlier_results<-rbind(malignant_outlier_results,tmp)
}
neoplasia_outlier_results<-as.data.frame(neoplasia_outlier_results)
colnames(neoplasia_outlier_results)[1]<-c("species")
neoplasia_outlier_results$species<-as.character(neoplasia_outlier_results$species)
neoplasia_outlier_results$neo_res_val<-as.numeric(neoplasia_outlier_results$neo_res_val)
neoplasia_outlier_results$neo_abs_res_val<-abs(neoplasia_outlier_results$neo_res_val)

malignant_outlier_results<-as.data.frame(malignant_outlier_results)
colnames(malignant_outlier_results)[1]<-c("species")
malignant_outlier_results$species<-as.character(malignant_outlier_results$species)
malignant_outlier_results$mal_res_val<-as.numeric(malignant_outlier_results$mal_res_val)
malignant_outlier_results$mal_abs_res_val<-abs(malignant_outlier_results$mal_res_val)


write.csv(neoplasia_outlier_results, '~/results/svl_neoplasia_single_bm_seperate_pl_studentized_residuals.csv',row.names=FALSE)
write.csv(malignant_outlier_results, '~/results/svl_malignant_single_bm_seperate_pl_studentized_residuals.csv',row.names=FALSE)
#########################################################################
neoplasia_sig_outlier_results<-neoplasia_outlier_results[neoplasia_outlier_results$neo_abs_res_val > 3,]
malignant_sig_outlier_results<-malignant_outlier_results[malignant_outlier_results$mal_abs_res_val > 3,]

sig_species<-as.vector(unique(c(neoplasia_sig_outlier_results$species,malignant_sig_outlier_results$species)))

filtered_data<-data[!(data$Species %in% sig_species),]

mod1<-MCMCglmm(cbind(neo,mal) ~ trait+trait:reptile+trait:log_rec + trait:log_bm  + trait:log_pl + trait:log_pl_rep, family=cbind("poisson","poisson"), pr=TRUE,prior = pr1, data=filtered_data, pl=TRUE, 
               rcov = ~ us(trait):units, random = ~us(trait):Species,nitt=number_iterations,thin=sample_interval,burnin=burnin_iterations, verbose=TRUE, DIC=TRUE, ginverse=list(Species=treeAinv))

saveRDS(mod1,'~/results/single_svl_seperate_pl.rds')

tmp<-as.data.frame(mod1$Sol[,c(1:12)])
class_vec<-rep(c("Amphibian","Reptilia"),each=nrow(tmp)*2)
growth<-rep(rep(c("neo","mal"),each=nrow(tmp)),2)
neo_amphib<-cbind(tmp[,1],tmp[,5],tmp[,7],tmp[,9])
mal_amphib<-cbind(tmp[,1]+tmp[,2],tmp[,6],tmp[,8],tmp[,10])
neo_reptil<-cbind(tmp[,1]+tmp[,3],tmp[,5],tmp[,7],tmp[,9]+tmp[,11])
mal_reptil<-cbind(tmp[,1]+tmp[,2]+tmp[,4],tmp[,6],tmp[,8],tmp[,10]+tmp[,12])

output<-as.data.frame(cbind(class_vec,growth,rbind(neo_amphib,mal_amphib,neo_reptil,mal_reptil)))
colnames(output)<-c("class","growth","intercept","records","bm_slope","pl_slope")
output[3:6]<-lapply(output[3:6],as.numeric)
write.csv(output,'~/results/single_svl_seperate_pl_parameter_estimates.csv',row.names=FALSE)
###################################################################################
#Single PL slope
pr1<- list(R=list(V=diag(2), nu=2, fix=2),B=list(mu=(rep(0,10)),V=diag(10)*1e8),
           G=list(G1=list(V=diag(2), nu=2, alpha.mu=cbind(0,0), alpha.V=diag(2)*25^2)))

species_list<-as.list(unique(data$Species))

cross_validate_seperate_slopes<-function(taxa_remove){
  data<-data[data$Species != taxa_remove,]
  mod1<-MCMCglmm(cbind(neo,mal) ~ trait+trait:reptile+trait:log_rec + trait:log_bm+ trait:log_pl, family=cbind("poisson","poisson"), pr=TRUE,prior = pr1, data=data, pl=TRUE, 
                 rcov = ~ us(trait):units, random = ~us(trait):Species,nitt=number_iterations,thin=sample_interval,burnin=burnin_iterations, verbose=FALSE, DIC=TRUE, ginverse=list(Species=treeAinv))
  saveRDS(mod1,file=paste0("~/results/n_1_single_svl_single_pl_",taxa_remove,".rds"))
}


no_cores<-50
cross_validation_results<-data.table::rbindlist(pbmclapply(species_list,cross_validate_seperate_slopes, mc.cores=no_cores,ignore.interactive = T))

neoplasia_outlier_results<-array(data=NA,c(0,2))
malignant_outlier_results<-array(data=NA,c(0,2))
for (i in 1:length(species_list)){
  mod1<-readRDS(paste0("~/results/n_1_single_svl_single_pl_",species_list[i],".rds"))
  pred<-as.data.frame(predict.MCMCglmm(mod1,data,interval="prediction",type="terms",marginal = NULL))
  pred_df<-cbind(pred[1:nrow(data),],pred[(nrow(data)+1):(2*nrow(data)),])
  pred_df<-pred_df[,c(1,4)]
  colnames(pred_df)<-c("neo_fit","mal_fit")
  
  output<-cbind(data,pred_df)
  output$neo_res<-output$neo - output$neo_fit
  neo_res_val<-output[output$Species == species_list[i],]$neo_res/sd(output[output$Species != species_list[i],]$neo_res)
  tmp<-cbind(species_list[i],neo_res_val)
  neoplasia_outlier_results<-rbind(neoplasia_outlier_results,tmp)
  
  output$mal_res<-output$mal - output$mal_fit
  mal_res_val<-output[output$Species == species_list[i],]$mal_res/sd(output[output$Species != species_list[i],]$mal_res)
  tmp<-cbind(species_list[i],mal_res_val)
  malignant_outlier_results<-rbind(malignant_outlier_results,tmp)
}
neoplasia_outlier_results<-as.data.frame(neoplasia_outlier_results)
colnames(neoplasia_outlier_results)[1]<-c("species")
neoplasia_outlier_results$species<-as.character(neoplasia_outlier_results$species)
neoplasia_outlier_results$neo_res_val<-as.numeric(neoplasia_outlier_results$neo_res_val)
neoplasia_outlier_results$neo_abs_res_val<-abs(neoplasia_outlier_results$neo_res_val)

malignant_outlier_results<-as.data.frame(malignant_outlier_results)
colnames(malignant_outlier_results)[1]<-c("species")
malignant_outlier_results$species<-as.character(malignant_outlier_results$species)
malignant_outlier_results$mal_res_val<-as.numeric(malignant_outlier_results$mal_res_val)
malignant_outlier_results$mal_abs_res_val<-abs(malignant_outlier_results$mal_res_val)

write.csv(neoplasia_outlier_results, '~/results/svl_neoplasia_single_bm_single_pl_studentized_residuals.csv',row.names=FALSE)
write.csv(malignant_outlier_results, '~/results/svl_malignant_single_bm_single_pl_studentized_residuals.csv',row.names=FALSE)
#########################################################################
neoplasia_sig_outlier_results<-neoplasia_outlier_results[neoplasia_outlier_results$neo_abs_res_val > 3,]
malignant_sig_outlier_results<-malignant_outlier_results[malignant_outlier_results$mal_abs_res_val > 3,]

sig_species<-as.vector(unique(c(neoplasia_sig_outlier_results$species,malignant_sig_outlier_results$species)))

filtered_data<-data[!(data$Species %in% sig_species),]

mod1<-MCMCglmm(cbind(neo,mal) ~ trait+trait:reptile+trait:log_rec + trait:log_bm + trait:log_pl, family=cbind("poisson","poisson"), pr=TRUE,prior = pr1, data=filtered_data, pl=TRUE, 
               rcov = ~ us(trait):units, random = ~us(trait):Species,nitt=number_iterations,thin=sample_interval,burnin=burnin_iterations, verbose=TRUE, DIC=TRUE, ginverse=list(Species=treeAinv))

saveRDS(mod1,'~/results/single_svl_single_pl.rds')

tmp<-as.data.frame(mod1$Sol[,c(1:10)])
class_vec<-rep(c("Amphibian","Reptilia"),each=nrow(tmp)*2)
growth<-rep(rep(c("neo","mal"),each=nrow(tmp)),2)
neo_amphib<-cbind(tmp[,1],tmp[,5],tmp[,7],tmp[,9])
mal_amphib<-cbind(tmp[,1]+tmp[,2],tmp[,6],tmp[,8],tmp[,10])
neo_reptil<-cbind(tmp[,1]+tmp[,3],tmp[,5],tmp[,7],tmp[,9])
mal_reptil<-cbind(tmp[,1]+tmp[,2]+tmp[,4],tmp[,6],tmp[,8],tmp[,10])

output<-as.data.frame(cbind(class_vec,growth,rbind(neo_amphib,mal_amphib,neo_reptil,mal_reptil)))
colnames(output)<-c("class","growth","intercept","records","bm_slope","pl_slope")
output[3:6]<-lapply(output[3:6],as.numeric)
write.csv(output,'~/results/single_svl_single_pl_parameter_estimates.csv',row.names=FALSE)

