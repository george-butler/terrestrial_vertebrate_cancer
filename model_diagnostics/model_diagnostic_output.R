library(MCMCglmm)

model<-readRDS("~/models/birds_and_mammals/single_bm_single_pl.rds")

pdf(file="model_diagnostics.pdf", width=20, height=15)
plot(model$Sol[,c(1:model$Fixed$nfl)])
plot(model$VCV)
dev.off()

