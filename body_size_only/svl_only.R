library(MCMCglmm)
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
data<-data[!(data$Species %in% c("Heloderma_suspectum", "Lampropeltis_getula", "Lampropeltis_triangulum", "Leptodactylus_fallax", "Pantherophis_guttatus","Peltophryne_lemur")),]
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
#########################################################################
#set up the model priors
pr1<- list(R=list(V=diag(2), nu=2, fix=2),B=list(mu=(rep(0,10)),V=diag(10)*1e8),
           G=list(G1=list(V=diag(2), nu=2, alpha.mu=cbind(0,0), alpha.V=diag(2)*25^2)))

mod1<-MCMCglmm(cbind(neo,mal) ~ trait+trait:reptile+trait:log_rec + trait:log_bm + trait:log_bm_rep, family=cbind("poisson","poisson"), pr=TRUE,prior = pr1, data=data, pl=TRUE,
               rcov = ~ us(trait):units, random = ~us(trait):Species,nitt=number_iterations,thin=sample_interval,burnin=burnin_iterations, verbose=TRUE, DIC=TRUE, ginverse=list(Species=treeAinv))

saveRDS(mod1,'~/results/svl_seperate_slopes_svl_only.rds')

#########################################################################
#repeat for SINGLE SLOPE MODEL
#set up the model priors
pr1<- list(R=list(V=diag(2), nu=2, fix=2),B=list(mu=(rep(0,8)),V=diag(8)*1e8),
           G=list(G1=list(V=diag(2), nu=2, alpha.mu=cbind(0,0), alpha.V=diag(2)*25^2)))

mod1<-MCMCglmm(cbind(neo,mal) ~ trait+trait:reptile+trait:log_rec + trait:log_bm, family=cbind("poisson","poisson"), pr=TRUE,prior = pr1, data=data, pl=TRUE, 
               rcov = ~ us(trait):units, random = ~us(trait):Species,nitt=number_iterations,thin=sample_interval,burnin=burnin_iterations, verbose=TRUE, DIC=TRUE, ginverse=list(Species=treeAinv))

saveRDS(mod1,'~/results/svl_cross_val/svl_single_slope_svl_only.rds')

######################################################################
pr1<- list(R=list(V=diag(2), nu=2, fix=2),B=list(mu=(rep(0,4)),V=diag(4)*1e8),
           G=list(G1=list(V=diag(2), nu=2, alpha.mu=cbind(0,0), alpha.V=diag(2)*25^2)))

null_mod<-MCMCglmm(cbind(neo,mal) ~ trait+trait:reptile, family=cbind("poisson","poisson"), pr=TRUE,prior = pr1, data=data, pl=TRUE, 
                   rcov = ~ us(trait):units, random = ~us(trait):Species,nitt=number_iterations,thin=sample_interval,burnin=burnin_iterations, verbose=TRUE, DIC=TRUE, ginverse=list(Species=treeAinv))

saveRDS(null_mod,'~/results/svl_intercept.rds')

r_squared_calc<-function(null_model,alt_model){
  #this function is based off of the code from Nakagawa at https://github.com/itchyshin/R2/blob/master/R/R_code_lmer_MCMCglmm.R
  #null model is intercept only model
  #alt model is model of interest
  mVdis <- log(1 + 1/exp(null_model$Sol[,1] + 0.5*(rowSums(alt_model$VCV))))
  vmVarF1<-numeric(nrow(alt_model$Sol))
  for(i in 1:nrow(alt_model$Sol)){
    Var<-var(as.vector(alt_model$Sol[i,c(1:alt_model$Fixed$nfl)] %*% t(alt_model$X)))
    vmVarF1[i]<-Var
  }
  R2m<-vmVarF1/(vmVarF1+alt_model$VCV[,1]+alt_model$VCV[,2]+alt_model$VCV[,3] + mean(mVdis))
  R2c<-(vmVarF1+alt_model$VCV[,1]+alt_model$VCV[,2])/(vmVarF1+alt_model$VCV[,1]+alt_model$VCV[,2]+alt_model$VCV[,3]+ mVdis)
  return(c(mean(R2m),mean(R2c)))#returns marginal rsquared first and then conditional rsquared
}

r_squared_calc(null_mod,mod1)



