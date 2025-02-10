library(MCMCglmm)
library(ape)
library(pbmcapply)
library(reshape2)

#this function is lifted from MCMCglmm.predict with the only modification being that it returns the entire posterior rather than just the upper and lower bound
full_posterior_predicit_mcmc_glmm<-function (object, newdata = NULL, marginal = object$Random$formula, 
                                             type = "response", interval = "none", level = 0.95, it = NULL, 
                                             posterior = "all", verbose = FALSE, approx = "numerical", 
                                             ...){
  rm.obs <- c()
  if (interval == "prediction") {
    post.pred <- t(simulate(object = object, nsim = nrow(object$Sol), 
                            newdata = newdata, marginal = marginal, type = type, 
                            it = it, posterior = posterior, verbose = verbose))
  }
  else {
    if (type %in% c("response", "terms") == FALSE) {
      stop("type must be response or terms")
    }
    if (interval %in% c("none", "confidence", "prediction") == 
        FALSE) {
      stop("interval must be none, confidence or prediction")
    }
    if (!is.null(posterior)) {
      if (!posterior %in% c("distribution", "mean", "mode", 
                            "all")) {
        stop("posterior argument must be either distribution, mean, mode or all")
      }
    }
    if (!is.null(marginal)) {
      if (!is(marginal, "formula")) {
        stop("marginal should be NULL or a formula")
      }
    }
    if (!is.null(it)) {
      if (length(it) > 1) {
        stop("it should be an integer")
      }
      if (it > nrow(object$Sol) | it < 1) {
        stop("it should be less than or equal to the number of iterations")
      }
    }
    rcomponents <- split.direct.sum(as.character(object$Random$formula)[2])
    mcomponents <- split.direct.sum(as.character(marginal)[2])
    if (is.null(object$meta)) {
      object$meta <- FALSE
    }
    if ((length(rcomponents) + object$meta) != length(object$Random$nrt)) {
      stop("if mev was used, add my_model$meta=TRUE and rerun. If not, the model is a covu model and predictions are not yet implemented for this type of model")
    }
    if (any(mcomponents %in% rcomponents == FALSE)) {
      stop("marginal formula does not correspond to model formula")
    }
    marginalise <- rep(as.numeric(rcomponents %in% mcomponents), 
                       object$Random$nrt)
    if (!is.null(newdata)) {
      suppressWarnings(object2 <- MCMCglmm(fixed = object$Fixed$formula, 
                                           random = object$Random$formula, rcov = object$Residual$formula, 
                                           family = object$Residual$original.family, data = newdata, 
                                           nitt = 1, thin = 1, burnin = 0, ginverse = object$ginverse, 
                                           verbose = FALSE, pr = any(marginalise == 0)))
      find.fixed <- match(colnames(object2$Sol)[1:object2$Fixed$nfl], 
                          colnames(object$Sol))
      find.random <- match(colnames(object2$Sol)[-c(1:object2$Fixed$nfl)], 
                           colnames(object$Sol))
      if (any(is.na(find.fixed))) {
        stop("model for newdata has fixed effects not present in original model")
      }
      if (any(!colnames(object$Sol)[1:object$Fixed$nfl] %in% 
              colnames(object2$Sol)[1:object2$Fixed$nfl])) {
        stop("model for newdata has fixed effects absent from the original model")
      }
      if (verbose) {
        if (any(is.na(find.random))) {
          missing.random <- colnames(object2$Sol)[which(is.na(find.random)) + 
                                                    object2$Fixed$nfl]
          warning(paste("model for newdata has random effects not present in original model:", 
                        paste(missing.random, collapse = ", ")))
        }
        missing.fixed <- which(!colnames(object$Sol)[1:object$Fixed$nfl] %in% 
                                 colnames(object2$Sol)[1:object2$Fixed$nfl])
        if (length(missing.fixed) > 0) {
          missing.fixed <- colnames(object$Sol)[1:object$Fixed$nfl][missing.fixed]
          warning(paste("original model has fixed effects not present in newdata:", 
                        paste(missing.fixed, collapse = ", ")))
        }
      }
      object2$Sol <- object$Sol[, c(find.fixed, find.random), 
                                drop = FALSE]
      find.vcv <- match(colnames(object2$VCV), colnames(object$VCV))
      if (any(is.na(find.vcv))) {
        stop("model for newdata has (co)variance terms not in original model")
      }
      object2$VCV <- object$VCV[, find.vcv, drop = FALSE]
      if (!is.null(object2$CP)) {
        find.cp <- match(colnames(object2$CP), colnames(object$CP))
        if (any(is.na(find.cp))) {
          stop("model for newdata has cutpoints not in original model")
        }
        object2$CP <- object$CP[, find.cp, drop = FALSE]
      }
      object <- object2
      rm(object2)
    }
    if (posterior == "mean") {
      object$VCV <- matrix(colMeans(object$VCV), 1, ncol(object$VCV))
      object$Sol <- matrix(colMeans(object$Sol), 1, ncol(object$Sol))
      it <- 1
    }
    if (posterior == "mode") {
      object$VCV <- matrix(posterior.mode(object$VCV, ...), 
                           1, ncol(object$VCV))
      object$Sol <- matrix(posterior.mode(object$Sol, ...), 
                           1, ncol(object$Sol))
      it <- 1
    }
    if (is.null(it)) {
      if (posterior == "distribution") {
        it <- sample(nrow(object$Sol), 1)
      }
      if (posterior == "all") {
        it <- 1:nrow(object$Sol)
      }
    }
    object$Sol <- object$Sol[it, , drop = FALSE]
    object$VCV <- object$VCV[it, , drop = FALSE]
    if (!is.null(object$Lambda)) {
      object$Lambda <- object$Lambda[it, , drop = FALSE]
    }
    if (!is.null(object$CP)) {
      object$CP <- object$CP[it, , drop = FALSE]
    }
    if (is.null(object$X)) {
      stop("fixed effect design matrix not saved: pass saveX=TRUE to MCMCglmm")
    }
    if (any(marginalise == 0) & dim(object$Sol)[2] == dim(object$X)[2]) {
      stop("posterior distribution of random effects not saved: pass pr=TRUE to MCMCglmm")
    }
    if (any(marginalise == 0) & is.null(object$Z)) {
      stop("random effect design matrix not saved: pass saveZ=TRUE to MCMCglmm")
    }
    if (is.null(object$Random$nfl) == FALSE) {
      st <- c(1, cumsum(object$Random$nrl * object$Random$nfl) + 
                1)
      st <- st[-length(st)]
      end <- cumsum(object$Random$nrl * object$Random$nfl)
      keep <- unlist(mapply(st[which(marginalise == 0)], 
                            end[which(marginalise == 0)], FUN = ":"))
    }
    else {
      keep <- NULL
    }
    if (!is.null(newdata) & !is.null(object$Random$nfl)) {
      missing.random <- which(is.na(object$Sol[1, ]))
      if (posterior == "mean" | posterior == "mode") {
        if (any(is.na(object$Sol))) {
          object$Sol[missing.random] <- 0
        }
      }
      else {
        dv <- 1:ncol(object$Sol)
        cnt <- 0
        for (j in 1:length(object$Random$nfl)) {
          nfl <- object$Random$nfl[j]
          nrl <- object$Random$nrl[j]
          dv[object$Fixed$nfl + cnt + 1:(nfl * nrl)] <- rep(diag(matrix(1:(nfl^2), 
                                                                        nfl, nfl)), each = nrl)
          cnt <- cnt + (nfl * nrl)
        }
        if (any(is.na(object$Sol))) {
          object$Sol[, missing.random, drop = FALSE] <- rnorm(nrow(object$Sol) * 
                                                                length(missing.random), 0, sqrt(object$VCV[, 
                                                                                                           dv[missing.random]]))
        }
      }
    }
    object$Sol <- object$Sol[, c(1:object$Fixed$nfl, object$Fixed$nfl + 
                                   keep), drop = FALSE]
    W <- cbind(object$X, object$Z)
    W <- W[, c(1:object$Fixed$nfl, object$Fixed$nfl + keep), 
           drop = FALSE]
    post.pred <- t(apply(object$Sol, 1, function(x) {
      (W %*% x)@x
    }))
    if (type == "response") {
      if (any(object$family != "gaussian" & object$family != 
              "cengaussian" & object$family != "ncst")) {
        post.var <- buildV(object, marginal = marginal, 
                           diag = TRUE, it = NULL, posterior = "all", 
                           ...)
      }
      super.trait <- 1:length(object$Residual$family)
      cnt <- 1
      i <- 1
      while (i <= length(super.trait)) {
        if (!grepl("hu|zi|za|multinomial", object$Residual$family[i])) {
          super.trait[i] <- cnt
          i <- i + 1
          cnt <- cnt + 1
        }
        else {
          if (grepl("multinomial", object$Residual$family[i])) {
            nm <- as.numeric(substr(object$Residual$family[i], 
                                    12, nchar(object$Residual$family[i]))) - 
              1
            super.trait[i + 1:nm - 1] <- cnt
            i <- i + nm
            cnt <- cnt + 1
          }
          else {
            super.trait[i] <- cnt
            super.trait[i + 1] <- cnt
            i <- i + 2
            cnt <- cnt + 1
          }
        }
      }
      normal.evd <- function(mu, v, approx) {
        if (approx == "numerical" || approx == "taylor2") {
          if (approx == "numerical") {
            int.foo <- function(x, mu, v) {
              exp(-exp(x)) * dnorm(x, mu, sqrt(v))
            }
            return(integrate(int.foo, qnorm(1e-06, mu, 
                                            sqrt(v)), qnorm(1 - 1e-06, mu, sqrt(v)), 
                             mu, v)[[1]])
          }
          else {
            return(exp(-exp(mu)) + 0.5 * v * exp(-exp(mu) + 
                                                   mu) * (exp(mu) - 1))
          }
        }
        else {
          stop(paste(approx, "approximation not implemented for this response"))
        }
      }
      normal.ztp <- function(mu, v, approx) {
        if (approx == "numerical" || approx == "taylor2") {
          if (approx == "numerical") {
            int.foo <- function(x, mu, v) {
              exp(x)/(1 - exp(-exp(x))) * dnorm(x, mu, 
                                                sqrt(v))
            }
            return(integrate(int.foo, qnorm(1e-04, mu, 
                                            sqrt(v)), qnorm(0.9999, mu, sqrt(v)), mu, 
                             v)[[1]])
          }
          else {
            return(exp(mu)/(1 - exp(-exp(mu))) + 0.5 * 
                     v * exp(exp(mu) + mu) * (1 - 2 * exp(exp(mu)) + 
                                                exp(2 * exp(mu)) + 3 * exp(mu) + 2 * exp(2 * 
                                                                                           mu) - 3 * exp(exp(mu) + mu) + exp(exp(mu) + 
                                                                                                                               2 * mu))/(exp(exp(mu)) - 1)^3)
          }
        }
        else {
          stop(paste(approx, "approximation not implemented for this response"))
        }
      }
      normal.ztb <- function(mu, v, size, approx) {
        if (approx == "numerical" || approx == "taylor2") {
          if (approx == "numerical") {
            int.foo <- function(x, mu, v) {
              (size * plogis(x)/(1 - (1 - plogis(x))^size)) * 
                dnorm(x, mu, sqrt(v))
            }
            return(integrate(int.foo, qnorm(1e-04, mu, 
                                            sqrt(v)), qnorm(0.9999, mu, sqrt(v)), mu, 
                             v)[[1]])
          }
          else {
            stop(paste(approx, "taylor approximation not implemented for this response"))
          }
        }
        else {
          stop(paste(approx, "approximation not implemented for this response"))
        }
      }
      normal.logistic <- function(mu, v, approx) {
        if (approx == "numerical") {
          int.foo <- function(x, mu, v) {
            plogis(x) * dnorm(x, mu, sqrt(v))
          }
          return(integrate(int.foo, qnorm(1e-04, mu, 
                                          sqrt(v)), qnorm(0.9999, mu, sqrt(v)), mu, 
                           v)[[1]])
        }
        if (approx == "diggle") {
          return(plogis(mu/sqrt(1 + v * (16 * sqrt(3)/(15 * 
                                                         pi))^2)))
        }
        if (approx == "mcculloch") {
          return(plogis(mu - 0.5 * v * tanh(mu * (1 + 
                                                    2 * exp(-0.5 * v))/6)))
        }
        if (approx == "taylor2") {
          return(plogis(mu) - 0.5 * v * exp(mu) * (exp(mu) - 
                                                     1)/(1 + exp(mu))^3)
        }
      }
      normal.clogistic <- function(mu, v, size, approx) {
        if (approx == "numerical" || approx == "taylor2") {
          if (approx == "numerical") {
            int.foo <- function(x, mu, v, size) {
              (1 - ((1 - plogis(mu))^size)) * dnorm(x, 
                                                    mu, sqrt(v))
            }
            return(integrate(int.foo, qnorm(1e-04, mu, 
                                            sqrt(v)), qnorm(0.9999, mu, sqrt(v)), mu, 
                             v, size)[[1]])
          }
          else {
            return((1 - ((1 - plogis(mu))^size)) - 0.5 * 
                     v * exp(mu) * size * (size * exp(mu) - 
                                             1)/(1 + exp(mu))^(2 + size))
          }
        }
        else {
          stop(paste(approx, "approximation not implemented for this response"))
        }
      }
      normal.multilogistic <- function(mu, v, approx) {
        if (approx == "numerical") {
          int.foo <- function(x, mu, v, i) {
            (exp(x[i])/(1 + sum(exp(x)))) * prod(dnorm(x, 
                                                       mu, sqrt(v)))
          }
          res <- 1:length(mu)
          for (i in 1:length(mu)) {
            res[i] <- cubature::adaptIntegrate(int.foo, 
                                               qnorm(1e-04, mu, sqrt(v)), qnorm(0.9999, 
                                                                                mu, sqrt(v)), mu = mu, v = v, i = i)[[1]]
          }
          return(res)
        }
        else {
          stop(paste(approx, "approximation not implemented for this response"))
        }
      }
      if (!is.null(object$Lambda)) {
        if (any(object$family != "gaussian")) {
          stop("sorry - prediction not yet available for sir models applied to non-Gaussian data")
        }
        I <- Diagonal(nrow(object$XL))
        post.pred <- t(sapply(1:nrow(object$Sol), function(i) {
          as.vector(solve(I - object$XL %*% kronecker(object$Lambda[i, 
          ], I), post.pred[i, ]))
        }))
      }
      if (!is.null(object$ThetaS)) {
        stop("sorry - prediction not yet available for theta_scale models")
      }
      if (any(object$family %in% c("poisson", "cenpoisson"))) {
        keep <- which(object$family %in% c("poisson", 
                                           "cenpoisson"))
        post.pred[, keep] <- exp(post.pred[, keep] + 
                                   0.5 * post.var[, keep])
      }
      if (any(object$family %in% c("ztpoisson"))) {
        keep <- which(object$family %in% c("ztpoisson"))
        post.pred[, keep] <- mapply(post.pred[, keep], 
                                    post.var[, keep], FUN = function(mu, v) {
                                      normal.ztp(mu, v, approx)
                                    })
      }
      if (any(object$family %in% c("exponential", "cenexponential", 
                                   "geometric", "cengeometric"))) {
        keep <- which(object$family %in% c("exponential", 
                                           "cenexponential", "geometric", "cengeometric"))
        post.pred[, keep] <- exp(-post.pred[, keep] + 
                                   0.5 * post.var[, keep])
      }
      if (any(object$family %in% c("ordinal"))) {
        nord <- unique(object$error.term[which(object$family == 
                                                 "ordinal")])
        cp.names <- substr(colnames(object$CP), 10, regexpr("\\.[1-9]$|\\.[1-9][0-9]$", 
                                                            colnames(object$CP)) - 1)
        cp.names <- match(cp.names, unique(cp.names))
        for (k in 1:length(nord)) {
          keep <- which(object$family %in% c("ordinal") & 
                          object$error.term == nord[k])
          CP <- cbind(-Inf, 0, object$CP[, which(cp.names == 
                                                   k), drop = FALSE], Inf)
          q <- matrix(0, dim(post.pred)[1], length(keep))
          for (i in 2:(dim(CP)[2] - 1)) {
            q <- q + (pnorm(CP[, i + 1] - post.pred[, 
                                                    keep], 0, sqrt(post.var[, keep] + 1)) - 
                        pnorm(CP[, i] - post.pred[, keep], 0, sqrt(post.var[, 
                                                                            keep] + 1))) * (i - 1)
          }
          post.pred[, keep] <- q
          rm(q)
        }
      }
      if (any(object$family %in% c("threshold"))) {
        nord <- unique(object$error.term[which(object$family == 
                                                 "threshold")])
        cp.names <- substr(colnames(object$CP), 10, regexpr("\\.[1-9]$|\\.[1-9][0-9]$", 
                                                            colnames(object$CP)) - 1)
        cp.names <- match(cp.names, unique(cp.names))
        for (k in 1:length(nord)) {
          keep <- which(object$family %in% c("threshold") & 
                          object$error.term == nord[k])
          CP <- cbind(-Inf, 0, object$CP[, which(cp.names == 
                                                   k), drop = FALSE], Inf)
          q <- matrix(0, dim(post.pred)[1], length(keep))
          for (i in 2:(dim(CP)[2] - 1)) {
            q <- q + (pnorm(CP[, i + 1] - post.pred[, 
                                                    keep], 0, sqrt(post.var[, keep])) - pnorm(CP[, 
                                                                                                 i] - post.pred[, keep], 0, sqrt(post.var[, 
                                                                                                                                          keep]))) * (i - 1)
          }
          post.pred[, keep] <- q
          rm(q)
        }
      }
      if (any(grepl("nzbinom", object$family))) {
        keep <- grep("nzbinom", object$family)
        size <- t(matrix(object$y.additional[keep, 1], 
                         length(keep), nrow(post.pred)))
        post.pred[, keep] <- mapply(post.pred[, keep], 
                                    post.var[, keep], size, FUN = function(mu, 
                                                                           v, size) {
                                      normal.clogistic(mu, v, size, approx)
                                    })
      }
      if (any(grepl("ncst", object$family))) {
        keep <- grep("ncst", object$family)
        scale <- t(matrix(object$y.additional[keep, 1], 
                          length(keep), nrow(post.pred)))
        df <- t(matrix(object$y.additional[keep, 2], 
                       length(keep), nrow(post.pred)))
        post.pred[, keep] <- exp(0.5 * log(df/2) + lgamma(df/2 - 
                                                            1/2) - lgamma(df/2)) * post.pred[, keep]
      }
      for (k in unique(super.trait)) {
        if (any(grepl("multinomial", object$Residual$family))) {
          keep <- which(object$error.term %in% which(super.trait == 
                                                       k))
          size <- object$y.additional[keep, 1]
          ncat <- sum(super.trait == k)
          for (j in 1:nrow(post.pred)) {
            prob <- matrix(post.pred[j, keep], length(keep)/sum(super.trait == 
                                                                  k), ncat)
            pvar <- matrix(post.var[j, keep], length(keep)/sum(super.trait == 
                                                                 k), ncat)
            if (ncat == 1) {
              post.pred[j, keep] <- t(sapply(1:nrow(prob), 
                                             function(x) {
                                               normal.logistic(prob[x, ], pvar[x, 
                                               ], approx)
                                             })) * size
            }
            else {
              post.pred[j, keep] <- t(sapply(1:nrow(prob), 
                                             function(x) {
                                               normal.multilogistic(prob[x, ], pvar[x, 
                                               ], approx)
                                             })) * size
            }
          }
        }
        if (any(grepl("hupoisson", object$Residual$family))) {
          keep <- which(object$error.term %in% which(super.trait == 
                                                       k))
          keep <- keep[-c(1:(length(keep)/2))]
          rm.obs <- c(rm.obs, keep)
          post.pred[, keep] <- mapply(post.pred[, keep], 
                                      post.var[, keep], FUN = function(mu, v) {
                                        1 - normal.logistic(mu, v, approx)
                                      })
          post.pred[, keep - length(keep)] <- post.pred[, 
                                                        keep] * mapply(post.pred[, keep - length(keep)], 
                                                                       post.var[, keep - length(keep)], FUN = function(mu, 
                                                                                                                       v) {
                                                                         normal.ztp(mu, v, approx)
                                                                       })
        }
        if (any(grepl("zapoisson", object$Residual$family))) {
          keep <- which(object$error.term %in% which(super.trait == 
                                                       k))
          keep <- keep[-c(1:(length(keep)/2))]
          rm.obs <- c(rm.obs, keep)
          post.pred[, keep] <- 1 - mapply(post.pred[, 
                                                    keep], post.var[, keep], FUN = function(mu, 
                                                                                            v) {
                                                      normal.evd(mu, v, approx)
                                                    })
          post.pred[, keep - length(keep)] <- post.pred[, 
                                                        keep] * mapply(post.pred[, keep - length(keep)], 
                                                                       post.var[, keep - length(keep)], FUN = function(mu, 
                                                                                                                       v) {
                                                                         normal.ztp(mu, v, approx)
                                                                       })
        }
        if (any(grepl("zipoisson", object$Residual$family[which(super.trait == 
                                                                k)]))) {
          keep <- which(object$error.term %in% which(super.trait == 
                                                       k))
          keep <- keep[-c(1:(length(keep)/2))]
          rm.obs <- c(rm.obs, keep)
          post.pred[, keep] <- mapply(post.pred[, keep], 
                                      post.var[, keep], FUN = function(mu, v) {
                                        1 - normal.logistic(mu, v, approx)
                                      })
          post.pred[, keep - length(keep)] <- post.pred[, 
                                                        keep] * exp(post.pred[, keep - length(keep)] + 
                                                                      post.var[, keep - length(keep)]/2)
        }
        if (any(grepl("zibinomial", object$Residual$family))) {
          keep <- which(object$error.term %in% which(super.trait == 
                                                       k))
          size <- object$y.additional[keep[1:(length(keep)/2)], 
                                      1]
          keep <- keep[-c(1:(length(keep)/2))]
          rm.obs <- c(rm.obs, keep)
          post.pred[, keep] <- mapply(post.pred[, keep], 
                                      post.var[, keep], FUN = function(mu, v) {
                                        1 - normal.logistic(mu, v, approx)
                                      })
          post.pred[, keep - length(keep)] <- post.pred[, 
                                                        keep] * mapply(post.pred[, keep - length(keep)], 
                                                                       post.var[, keep - length(keep)], FUN = function(mu, 
                                                                                                                       v) {
                                                                         normal.logistic(mu, v, approx)
                                                                       }) * size
        }
        if (any(grepl("hubinomial", object$Residual$family))) {
          keep <- which(object$error.term %in% which(super.trait == 
                                                       k))
          size <- object$y.additional[keep[1:(length(keep)/2)], 
                                      1]
          keep <- keep[-c(1:(length(keep)/2))]
          rm.obs <- c(rm.obs, keep)
          post.pred[, keep] <- mapply(post.pred[, keep], 
                                      post.var[, keep], FUN = function(mu, v) {
                                        1 - normal.logistic(mu, v, approx)
                                      })
          post.pred[, keep - length(keep)] <- post.pred[, 
                                                        keep] * mapply(post.pred[, keep - length(keep)], 
                                                                       post.var[, keep - length(keep)], rep(size, 
                                                                                                            each = nrow(post.pred)), FUN = function(mu, 
                                                                                                                                                    v, size) {
                                                                                                              normal.ztb(mu, v, size, approx)
                                                                                                            })
        }
      }
    }
  }
  return(t(post.pred))
}
prediction_posterior_formatting<-function(data,sl,response){
  if (response == "neo"){
    df<-data[1:(nrow(data)/2),]
  }
  if (response == "mal"){
    df<-data[(1+(nrow(data)/2)):nrow(data),]
  }
  colnames(df)<-sub('.', '', colnames(df))
  df$species<-unlist(sl)
  tmp<-melt(df, id.vars=c("species"))
  colnames(tmp)[2]<-c("iteration")
  tmp$iteration<-as.numeric(tmp$iteration)
  t<-tmp[with(tmp, order(species, iteration)), ]
  return(t)
}

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
data<-data[data$Class %in% c("Aves","Mammalia"),]
#########################################################################
#set up the model priors for neoplasia analysis
data$Class<-as.factor(data$Class)
class_levels<-levels(data$Class)
number_iterations<-1e7
sample_interval<-1e3
burnin_iterations<-9e6

number_of_samples<-((number_iterations - burnin_iterations) / sample_interval)
data$mammal = ifelse(data$Class == "Mammalia", 1, 0)
data$Class<-as.character(data$Class)

df_full<-data
data<-data[!(data$Species %in% c("Atelerix_albiventris","Gallus_gallus","Melopsittacus_undulatus","Mustela_putorius","Rattus_norvegicus")),]
#########################################################################
#set up the model priors for malignant analysis
pr1<- list(R=list(V=diag(2), nu=2, fix=2),B=list(mu=(rep(0,10)),V=diag(10)*1e8),
           G=list(G1=list(V=diag(2), nu=2, alpha.mu=cbind(0,0), alpha.V=diag(2)*25^2)))
#########################################################################
species_list<-as.list(unique(df_full$Species))

mod1<-MCMCglmm(cbind(neo,mal) ~ trait+trait:mammal+trait:log_rec + trait:log_bm + trait:log_pl, family=cbind("poisson","poisson"), pr=TRUE,prior = pr1, data=data, pl=TRUE, 
               rcov = ~ us(trait):units, random = ~us(trait):Species,nitt=number_iterations,thin=sample_interval,burnin=burnin_iterations, verbose=TRUE, DIC=TRUE, ginverse=list(Species=treeAinv))

saveRDS(mod1,file=paste0("~/results/exceptional_analysis_birds_and_mammals_single_bm_single_pl.rds"))

mod1<-readRDS("~/results/exceptional_analysis_birds_and_mammals_single_bm_single_pl.rds")

prediction_data<-as.data.frame(full_posterior_predicit_mcmc_glmm(mod1,df_full,interval="prediction",type="response",marginal = NULL))
neoplasia_prediction_data<-prediction_posterior_formatting(prediction_data,species_list,"neo")
malignant_prediction_data<-prediction_posterior_formatting(prediction_data,species_list,"mal")

prediction_results = array(data=NA,c(0,5))
for (i in species_list){
  neo_results<-c(i,df_full[df_full$Species == i,]$Class,"neo",median(neoplasia_prediction_data[neoplasia_prediction_data$species == i,]$value),df_full[df_full$Species == i,]$neo)
  mal_results<-c(i,df_full[df_full$Species == i,]$Class,"mal",median(malignant_prediction_data[malignant_prediction_data$species == i,]$value),df_full[df_full$Species == i,]$mal)
  prediction_results = rbind(prediction_results,neo_results)
  prediction_results = rbind(prediction_results,mal_results)
}
prediction_results<-as.data.frame(prediction_results)
colnames(prediction_results)<-c("species","class","response","median_prediction","actual_val")
prediction_results[,c(4:5)]<-lapply(prediction_results[,c(4:5)],as.numeric)
prediction_results$diff<-prediction_results$actual_val - prediction_results$median_prediction
prediction_results$outcome<-"M"
prediction_results[prediction_results$diff < 0,]$outcome<-"L"
prediction_results[prediction_results$diff == 0,]$outcome<-"S"
prediction_results$abs_diff<-abs(prediction_results$diff)

write.csv(prediction_results,"~/results/exceptional_mammals_and_birds_results_bm_and_pl.csv",row.names=FALSE)
