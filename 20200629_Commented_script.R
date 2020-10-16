##this code allows to perform a simulation to assess performance in terms of
## sensitivity and specificity of prespecified statistical methods
##used to find the true predictors of an health among the exposome under diverse
##causal structures involving a causal link from the exposome to an outcome.
##Some of them use an intermediary layer, some not. It contains 5 ##parts:

##1. defining the functions allowing to generate a realistic dataset of exposome
##intermediate layer and outcome. The three layers(E, M and Y) can be linearly 
##related to simulate various causal structures. 

#It needs real datasets (exposome/intermediate layer) as inputs, as well 
##as parameters allowing to define the association within the three layers (number of
#predictors, variability explained, correlation..)

##2. defining the methods assessed

##3. defining some functions used to assess methods performance

##4. defining the simulation function, which, for a given scenario, generates
##the datasets, applies the methods and assess their performance. This function
##allows to parallelize the simulation.

##5. runnning the simulation itself with parallelization, repeating X times the
##function defined in 4. for each scenario and saving the results.
###############################################################################

##load packages
library(mvtnorm)
library(boot)
library(parallel)
library(reshape)
library(glmnet)
library(DSA)

#####################################
##1. define the generating functions
####################################
simulator <-
  function(E_true, ##real exposome
           M_true, ## real intermediate layer, eg methylome
           R2_tot = 0.5, ##total varibility of the outcome explained by
           #the predictors of E and M
           propmE = 0, ##proportion of variables of M affected by E 
           ##without affecting Y
           propmEY = 0.1, ##proportion of variables of M acting as mediators 
           ##between E and Y
           propmY = 0, ##proportion of variables of M not affected by E 
           ##but affecting Y.
           BetamEY = 0.1, ##coefficient of the effect of M on Y for variables 
           ##of M acting as mediators between E and Y. It can be an unique value
           ##or a vector of length n_mEY
           BetamY = 0.1, ##coefficient of the effect of M on Y for variables 
           ##of M acting as mediators between E and Y. It can be an unique value
           ##or a vector of length n_mY
           n_mE = NULL, ##alternative way to specify the three previous set 
           ##of predictors: #directly giving the number of predictors 
           n_mEY = NULL,
           n_mY = NULL,
           n_EmE = 0, ##number of exposures affecting variables of M without 
           ##acting through M on the outcome
           n_EmEY = 3,##number of exposures having an effect on Y through M
           n_Ey = 0, ##number of exposures acting directly on M
           ##the 4 next variables specify the intersection between the different
           ## sets of predictors in E
           n_EmE_U_n_EmEY = 0,
           n_Ey_U_n_EmE = 0,
           n_Ey_U_n_EmEY = 0,
           n_Ey_U_n_EmE_U_n_EmEY = 0,
           ##the three next variables specify the coefficients of effects for 
           ##the different sets of predictors in E. They can be an unique value
           ##or a vector of size respectively n_EmE, n_EmEY, n_Ey...
           BetaEmE = 0.1, 
           BetaEmEY = 0.1,
           BetaEy = 0.1,
           test_and_training = TRUE ##generating only a training set or
           ##alternatively also a test set of same size
           ) {
    ##check of consistency within the input
    if ((BetaEmEY == 0 &
         n_EmEY != 0)) {
      stop("error: BetaEmEY and n_EmEY not compatible")
    }
    if ((n_mEY == 0 &
         n_EmEY != 0)) {
      stop("error: n_mEY and n_mEY not compatible")
    }
    if ((n_mEY != 0 &
         n_EmEY == 0)) {
      stop("error: n_mEY and n_mEY not compatible")
    }
    ##sampling with replacemnt the real data for exposome
    data.X <- as.data.frame(dataExp_true)
    names_row <- rownames(data.X)
    data.X <-
      data.X[sample(1:nrow(data.X), 2 * nrow(data.X), replace = TRUE),]
    rownames(data.X) <- c(names_row, sprintf('boot%s', names_row))
    dataExp <- data.X
    remove(data.X)
    ##sampling with replacement the real intermediate data 
    data.X <- as.data.frame(M1_true)
    names_row <- rownames(data.X)
    data.X <-
      data.X[sample(1:nrow(data.X), 2 * nrow(data.X), replace = TRUE),]
    rownames(data.X) <- c(names_row, sprintf('boot%s', names_row))
    M1 <- data.X
    
    ##setting linear relationship between E and M
    
    ##E on variables on M having no effect on Y
    
    ##converting if necessary proportion of predictors to number of predictors
    if (is.null(n_mE)) {
      n_mE <- floor(propmE * nrow(M1))
    }
    ##multiplicity constraint: check if number of variables n_mE is a multiple of
    ##n_EmE
    if (n_mE != 0){
    n_pat_mE <- (n_mE / n_EmE)
    if (trunc(n_pat_mE) != n_pat_mE) {
      stop("error: number of cpgs explained by E is not a multiple of
           the number of predictors from E")
    }
    }
    ###creating the vector of predictors
    mE_predictors <- list()
    ##generating vector of effect coefficients
    Npred_mE_E <- rep(list(n_EmE), n_mE)
    Betapred_mE_E <- list()
    for (t in  seq(length.out = n_EmE)) {
      list_temp <- as.list(rep(0, n_EmE))
      if (length(BetaEmE) == 1) {
        list_temp[[t]] <- BetaEmE
      } else{
        if (length(BetaEmE) != n_EmE) {
          stop("error: Betas for E explaining M are not consistent
               with the number of predictors")
        }
        list_temp[[t]] <- BetaEmE[[t]]
      }
      list_Beta_pred_mE_E_by_pat <- rep(list(list_temp), n_pat_mE)
      Betapred_mE_E <- c(Betapred_mE_E, list_Beta_pred_mE_E_by_pat)
    }
    ##random sampling of causal exposures
    ind_mE_E <- sample(ncol(dataExp), n_EmE)
    ##random sampling of variables of M affected
    ind_mE_M <- sample(ncol(M1), n_mE)
    ##adding a linear effect from E on each variable of M
    for (k in seq(length.out = n_mE)) {
      list_temp <-
        list(ind_mE_M[k],
             colnames(M1)[ind_mE_M[k]],
             Npred_mE_E[[k]],
             Betapred_mE_E[[k]],
             ind_mE_E)
      list_temp <- c(list_temp, list(colnames(dataExp)[list_temp[[5]]]))
      names(list_temp) <-
        c(
          "indice_cpg",
          "cpg",
          "nb_exp_predictors",
          "Beta_exp_predictors",
          "ind_exp_predictors",
          "name_exp_predictors"
        )
      cpg_temp <-
        simResponseSimple(
          met = dataExp,
          Nmet = list_temp[[3]],
          beta = unlist(list_temp[[4]]),
          list_temp[[5]]
        )
      M1[, list_temp[[1]]] <-
        as.numeric(M1[, list_temp[[1]]] + cpg_temp$resp)
      list_temp[[7]] <-
        estimatedR2(dataExp, list_temp[[6]], M1[, list_temp[[1]], drop = FALSE])
      names(list_temp)[[7]] <- "R2"
      mE_predictors <- c(mE_predictors, list(list_temp))
      remove(list_temp)
      remove(cpg_temp)
    }
    ##empirical estimation of mean R2 (mean variability of M affected by E 
    ##explained by E)
    if (n_mE != 0) {
      R2_mean_mE <-
        mean(unlist(lapply(mE_predictors, function(X)
          X$R2$r.squared)))
    } else{
      R2_mean_mE = 0
    }
    
    
    
    
    ##E on variables on M mediating effect on Y
    ##converting if necessary proportion of predictors to number of predictors
    if (is.null(n_mEY)) {
      n_mEY <- floor(propmEY * nrow(M1))
    }
    ##multiplicity constraint: check if number of variables n_mEY is a multiple of
    ##n_EmEY
    if (n_mEY != 0){
    n_pat_mEY <- (n_mEY / n_EmEY)
    if (trunc(n_pat_mEY) != n_pat_mEY) {
      stop("error: number of cpgs explained by E explaining Y is not a
           multiple of the number of their predictors from E")
    }
    }
    ###creating the vector of predictors
    mEY_predictors <- list()
    ##generating vector of effect coefficients
    Npred_mEY_E <- rep(list(n_EmEY), n_mEY)
    Betapred_mEY_E <- list()
    for (t in  seq(length.out = n_EmEY)) {
      list_temp <- as.list(rep(0, n_EmEY))
      if (length(BetaEmEY) == 1) {
        list_temp[[t]] <- BetaEmEY
      } else{
        if (length(BetaEmEY) != n_EmEY) {
          stop("error: Betas for E explaining M are not consistent
               with the number of predictors")
        }
        list_temp[[t]] <- BetaEmEY[[t]]
      }
      list_Beta_pred_mEY_E_by_pat <- rep(list(list_temp), n_pat_mEY)
      Betapred_mEY_E <- c(Betapred_mEY_E, list_Beta_pred_mEY_E_by_pat)
    }
    ##random sampling of causal exposures 
    if ((n_EmE) != 0) {
      ind_mEY_E_in_pred_mE = sample(ind_mE_E, n_EmE_U_n_EmEY)
      ind_mEY_E <-
        c(ind_mEY_E_in_pred_mE, sample((1:ncol(dataExp))[-ind_mE_E], n_EmEY - n_EmE_U_n_EmEY))
    } else{
      ind_mEY_E <- sample(ncol(dataExp), n_EmEY)
    }
   ##random sampling of variables of M
    if (n_mE != 0) {
      ind_mEY_M <-
        sample((1:ncol(M1))[-ind_mE_M], n_mEY)
    } else{
      ind_mEY_M <- sample(ncol(M1), n_mEY)
    }
    
    ##adding a linear effect from E on each variable of M
        for (k in seq(length.out = n_mEY)) {
      list_temp <-
        list(ind_mEY_M[k],
             colnames(M1)[ind_mEY_M[k]],
             Npred_mEY_E[[k]],
             Betapred_mEY_E[[k]],
             ind_mEY_E)
      list_temp <- c(list_temp, list(colnames(dataExp)[list_temp[[5]]]))
      names(list_temp) <-
        c(
          "indice_cpg",
          "cpg",
          "nb_exp_predictors",
          "Beta_exp_predictors",
          "ind_exp_predictors",
          "name_exp_predictors"
        )
      cpg_temp <-
        simResponseSimple(
          met = dataExp,
          Nmet = list_temp[[3]],
          beta = unlist(list_temp[[4]]),
          list_temp[[5]]
        )
      M1[, list_temp[[1]]] <-
        as.numeric(M1[, list_temp[[1]]] + cpg_temp$resp)
      list_temp[[7]] <-
        estimatedR2(dataExp, list_temp[[6]], M1[, list_temp[[1]], drop = FALSE])
      names(list_temp)[[7]] <- "R2"
      mEY_predictors <- c(mEY_predictors, list(list_temp))
      remove(list_temp)
      remove(cpg_temp)
    }
    ##empirical estimation of mean R2 (mean variability of M affected by E 
    ##explained by E)  
    if (n_mEY != 0) {
      R2_mean_mEY <-
        mean(unlist(lapply(mEY_predictors, function(X)
          X$R2$r.squared)))
    } else{
      R2_mean_mEY = 0
    }
    
    ##variables of M not affected by E having an effect on Y
    if (is.null(n_mY)) {
      n_mY <- floor(propmY * nrow(M1))
    }
    ##generating vector of effect coefficients 
  if (length(BetamY) == 1) {
      Betapred_yM_M <- rep(BetamY, n_mY)
    } else{
      if (length(BetamY) != n_mY) {
        stop(
          "error: Betas for M explaining Y not explained by E are not consistent with the number of predictors"
        )
      }
      Betapred_yM_M <- BetamY
    }
    ##if there is not effect of mY, generating an empty yM
    if (n_mY == 0) {
      yM = list(as.matrix(rep(0, nrow(M1)), ncol = 1), NULL, NULL)
      names(yM) <- c("resp", "beta")
    } else{
      ##if there is an effect, generating yM: part of the outcome which is a 
      ##linear combination of variables of M not affected by E
      if (n_mEY == 0 & n_mE == 0) {
        ind_yM_M <- sample(ncol(M1), n_mY)
      } else{
        ind_yM_M <- sample((1:ncol(M1))[-c(ind_mEY_M, ind_mE_M)], n_mY)
      }
      yM <-
        simResponseSimple(
          met = M1,
          Nmet = length(ind_yM_M),
          beta = Betapred_yM_M,
          cpg = ind_yM_M
        )
    }
    
    
  ##variables of M mediating an effect of E on Y
  ##generating vector of effect coefficients 
    if (length(BetamEY) == 1) {
      Betapred_yME_M <- rep(BetamEY, n_mEY)
    } else{
      if (length(BetamEY) != n_mEY) {
        stop(
          "error: Betas for M explaining Y explained by E are not consistent with the number of predictors"
        )
      }
      Betapred_yM_M <- BetamY
    }
#if there is an effect, generating yME: part of the outcome which is a 
##linear combination of variables of M  affected by E
    ind_yME_M <- ind_mEY_M
    if (n_mEY != 0) {
      yME <-
        simResponseSimple(
          met = M1,
          Nmet = length(ind_yME_M),
          beta = Betapred_yME_M,
          cpg = ind_yME_M
        )
    } else{
      yME = list(as.matrix(rep(0, nrow(M1)), ncol = 1), NULL, NULL)
      names(yME) <- c("resp", "beta")
    }
    
    ##direct effect of E on Y
    ##random sampling of exposures with respect to the specification of 
    #intersections between the different groups of exposures having different 
    ##effects
    if (n_EmE == 0 & n_EmEY == 0) {
      ind_yE_E <- sample(ncol(dataExp), n_Ey)
      ind_yE_E_shared_mEY <- integer(0)
      ind_yE_E_shared_mE <- integer(0)
      ind_yE_E_shared_mE_mEY <- integer(0)
    } else{
      ind_yE_E <-
        sample((1:ncol(dataExp))[-unique(c(ind_mE_E, ind_mEY_E))],
               n_Ey - n_Ey_U_n_EmE - n_Ey_U_n_EmEY + n_Ey_U_n_EmE_U_n_EmEY)
      if (n_Ey_U_n_EmE_U_n_EmEY != 0) {
        ind_yE_E_shared_mE_mEY <-
          sample(intersect(ind_mE_E, ind_mEY_E), n_Ey_U_n_EmE_U_n_EmEY)
      } else{
        ind_yE_E_shared_mE_mEY <- integer(0)
      }
      if (n_Ey_U_n_EmE != 0) {
        ind_yE_E_shared_mE <-
          sample(ind_mE_E[!ind_mE_E %in% intersect(ind_mE_E, ind_mEY_E)], n_Ey_U_n_EmE -
                   n_Ey_U_n_EmE_U_n_EmEY)
      } else{
        ind_yE_E_shared_mE <- integer(0)
      }
      if (n_Ey_U_n_EmEY != 0) {
        ind_yE_E_shared_mEY <-
          sample(ind_mEY_E[!ind_mEY_E %in% intersect(ind_mE_E, ind_mEY_E)], n_Ey_U_n_EmEY -
                   n_Ey_U_n_EmE_U_n_EmEY)
        
      } else{
        ind_yE_E_shared_mEY <- integer(0)
      }
      ind_yE_E <-
        c(ind_yE_E,
          ind_yE_E_shared_mE,
          ind_yE_E_shared_mEY,
          ind_yE_E_shared_mE_mEY)
    }
    ##generating vector of effects coefficients
    if (length(BetaEy) == 1) {
      Betapred_yE_E <- rep(BetaEy, n_Ey)
    } else{
      if (length(BetaEy) != n_Ey) {
        stop(
          "error: Betas for M explaining Y not explained by E are not 
          consistent with the number of predictors"
        )
      }
      Betapred_yE_E <- BetaEy
    }
    ##generating yE: part of the outcome which is a 
    ##linear combination of exposures
    yE <-
      simResponseSimple(
        met = dataExp,
        Nmet = length(ind_yE_E),
        beta = Betapred_yE_E,
        cpg = ind_yE_E
      )
    
    ##Creating the final Y by adding a gaussian according to the variability 
    ##wanted to the differents parts of Y already created
    Y <- yE$resp + yME$resp + yM$resp
    if (!is.na(R2_tot)) {
      if (((R2_tot) != 0)) {
        sigma <- var(Y) * (1 / R2_tot - 1)
      } else{
        R2 = 0.00000001
        sigma <- var(Y) * (1 / R2_tot - 1)
      }
      Y <- as.matrix(Y + rnorm(length(Y), mean(Y), sqrt(sigma)), ncol = 1)
    }
    
    ##extracting all indirect predictors of Y from E 
    ##and computing the corresponding betas
    datapred <- data.frame(exp = character(0), beta = numeric(0))
    for (k in seq(length.out = n_mEY)) {
      datapred_temp <-
        cbind(
          exp = unlist(mEY_predictors[[k]]$name_exp_predictors),
          beta = unlist(mEY_predictors[[k]]$Beta_exp_predictors) * yME$beta[k]
        )
      datapred <- rbind(datapred, datapred_temp)
    }
    class(datapred$beta) <- "numeric"
    yME_E <-
      list(
        resp = yME$resp,
        beta = sapply(unique(datapred$exp), function(X)
          sum(datapred[datapred$exp == X, ]$beta)),
        predictors = unique(datapred$exp)
      )
    ##extracting all predictors of Y from E (direct and indirect) 
    ##and computing the corresponding betas
    datapred <-
      data.frame(cbind(
        exp = c(
          as.character(yE$predictors),
          as.character(yME_E$predictors)
        ),
        beta = c(as.numeric(yE$beta), as.numeric(yME_E$beta))
      ))
    datapred$beta <- as.numeric(as.character(datapred$beta))
    yE_ME <-
      list(
        resp = Y,
        beta = sapply(unique(datapred$exp), function(X)
          sum(datapred[datapred$exp == X, ]$beta)),
        predictors = unique(datapred$exp)
      )
    remove(datapred)
    ##extracting all predictors of Y from M
    ##and computing the corresponding betas
    datapred <-
      data.frame(cbind(
        cpg = c(as.character(yM$predictors), as.character(yME$predictors)),
        beta = c(as.numeric(yM$beta), as.numeric(yME$beta))
      ))
    datapred$beta <- as.numeric(as.character(datapred$beta))
    yM_ME <-
      list(
        resp = Y,
        beta = sapply(unique(datapred$cpg), function(X)
          sum(datapred[datapred$cpg == X, ]$beta)),
        predictors = unique(datapred$cpg)
      )
    remove(datapred)
    
    ##computing the effective R2
    if ((n_mEY + n_mE) != 0) {
      R2_mean_M_E <-
        (n_mEY * R2_mean_mEY + n_mE * R2_mean_mE) / ((n_mEY + n_mE))
    } else{
      R2_mean_M_E = 0
    }
    R2 <-
      list(
        BMI_all_exp = estimatedR2(dataExp, yE_ME$predictors, Y)$r.squared,
        BMI_all_M = estimatedR2(M1, yM_ME$predictors, Y)$r.squared,
        mean_M_E = R2_mean_M_E
      )
    
    ##creating a list to return
    
    resultats <-
      list(
        Y_train = Y[1:(nrow(M1) / 2), , drop = FALSE],
        E_train = dataExp[1:(nrow(M1) / 2), , drop = FALSE],
        M_train = M1[1:(nrow(M1) / 2), , drop = FALSE],
        Y_test = Y[(nrow(M1) / 2):nrow(M1), , drop = FALSE],
        E_test = dataExp[(nrow(M1) / 2):nrow(M1), , drop = FALSE],
        M_test = M1[(nrow(M1) / 2):nrow(M1), , drop = FALSE],
        yM_E = yME_E,
        y_E = yE_ME,
        yE_E = yE,
        yM_M = yM_ME,
        R2 = R2,
        list_mY_predictor = as.character(yM$predictors),
        list_mE_Y_predictor = as.character(yME$predictors)
      )
    return(resultats)
  }

##function used to create a linear response
simResponseSimple <- function(met, ##matrix of potential predictors
                              Nmet = NA, ##number of predictors
                              beta = NULL, ##vector of effects
                              cpg = NULL) { ##optionnal: directly specifying 
                            ##some of the indexes of predictors
  if (all(c(is.na(Nmet), is.na(cpg))) == TRUE) {
    return (list(
      resp = as.matrix(rep(0, nrow(met)), ncol = 1),
      beta = NA,
      predictors = NA
    ))
  }
  temp <- Nmet - length(cpg)
  if (temp != 0) {
    wh <- sample((1:ncol(met)[-cpg]), temp)
     wh <- c(cpg, wh)
  } else{
    wh <- cpg
  }
  CovMat <- as.matrix(met[, wh])
  colnames(CovMat) <- colnames(met)[wh]
  # computing the response
  mean <- CovMat %*% matrix(beta, ncol = 1)
  rownames(mean)<-rownames(met)
  names(beta) <- colnames(CovMat)
  return (list(
    resp = mean,
    beta = beta,
    predictors = colnames(met)[wh]
  ))
}

##fonction to estimate R2 from a dataframe of potential predictors, a vector of 
##predictors names and the outcome
estimatedR2 <- function(X, truepred, Y) {
  if ("y" %in% truepred) {
    stop("error: one of the true predictors is named y")
  }
  if (ncol(Y) != 1) {
    stop("error:Y is multidimensionnal")
  }
  if (nrow(X) != nrow(Y)) {
    stop("error not the same number of rows")
  }
  if (isTRUE(all.equal(rownames(X), rownames(Y))) == FALSE) {
    stop("error individuals are not ordered similarly in X and Y")
  }
  if (all(truepred %in% colnames(X))) {
    data <- X[, colnames(X) %in% truepred, drop = FALSE]
    data <- cbind(Y, data)
    colnames(data)[1] <- "y"
    mod <- lm(y ~ ., as.data.frame(data))
    toselect.x <- summary(mod)$coeff[-1, 4]
    r <-
      list(summary(mod)$r.squared,
           summary(mod)$adj.r.squared,
           names(toselect.x)[toselect.x == TRUE])
    names(r) <- c("r.squared", "adj.r.squared", "pred")
    return(r)
  } else{
    stop("error: X does not countain all true predictors")
    
  }
}


###########################################
##2. defining the methods to test
###########################################


###agnostic methods

##function to compute residuals of a linear model if covariates are specified
getresiduals_2df<-function(data_Y_in,data_covar_in,name_Y,covar){
  data_covar<-data_covar_in[,colnames(data_covar_in)%in%covar,drop=FALSE]
  data_Y<-data_Y_in[rownames(data_Y_in)%in%rownames(data_covar),colnames(data_Y_in)==name_Y,drop=FALSE]
  data_covar<-data_covar[rownames(data_covar)%in%rownames(data_Y),,drop=FALSE]
  data_covar<-data_covar[rownames(data_Y),,drop=FALSE]
  data_output<-data_Y
  data<-cbind(data_Y,data_covar)
  mod<-lm(data=data)
  data_output[,1]<-as.data.frame(residuals(mod))
  return(data_output)
}

###ExWAS
ewas<-function(data_Xs_in=NULL,data_Y_in=NULL,name_Y,data_covar_in=NULL,covar=character(0),corr="BH"){
  require(parallel)
  if (length(covar)>0){
    data_covar<-data_covar_in[rownames(data_covar_in)%in%rownames(data_Y_in)&rownames(data_covar_in)%in%rownames(data_Xs_in),colnames(data_covar_in)%in%covar,drop=FALSE]
    data_Y<-data_Y_in[rownames(data_Y_in)%in%rownames(data_covar)&rownames(data_Y_in)%in%rownames(data_Xs_in),colnames(data_Y_in)==name_Y,drop=FALSE]
    data_Xs<-data_Xs_in[rownames(data_Xs_in)%in%rownames(data_covar)&rownames(data_Xs_in)%in%rownames(data_Y),,drop=FALSE]
    data_covar<-data_covar[rownames(data_Y),,drop=FALSE]
    data_Xs<-data_Xs[rownames(data_Y),,drop=FALSE]
    data_Y<-getresiduals_2df(data_Y,data_covar,name_Y,covar)
  }else{
    data_Y<-data_Y_in[rownames(data_Y_in)%in%rownames(data_Xs_in),colnames(data_Y_in)==name_Y,drop=FALSE]
    data_Xs<-data_Xs_in[rownames(data_Xs_in)%in%rownames(data_Y_in),,drop=FALSE]
    data_Xs<-data_Xs[rownames(data_Y),,drop=FALSE]
  }
  if (is.null(data_Y)==TRUE |is.null(data_Xs)==TRUE|!(name_Y%in%colnames(data_Y))) {
    stop("Inconsistent data")
  }
  
  ##computing p.values
  p.values <- lapply(1:ncol(data_Xs), function(x,data_Xs){
    c(colnames(data_Xs)[x],confint(lm(Y~var1,
                                      data=data.frame(cbind(var1=data_Xs[,x],Y=data_Y[,1]))))[2,],
      summary(lm(Y~var1, data=data.frame(cbind(var1=data_Xs[,x],Y=data_Y[,1]))))$coefficients[2,])},
    data_Xs)
  if (length(p.values)>1){
    p.values <- cbind(matrix(unlist(p.values), ncol = 7, byrow = TRUE)[,-6])
    
  }else{
    p.values <-  as.data.frame(t(as.data.frame(unlist(p.values)))[,-6,drop=FALSE])
  }
  p.values<-as.data.frame(p.values)
  colnames(p.values) <- c("var","conf - 2.5%","conf - 97.5%", "Est","Sd","pVal")
  p.values <- p.values[p.values$var!="Intercept",]
  p.values$pVal<-as.numeric(as.character(p.values$pVal))
  p.values.adj<-p.values
  pVal <- as.numeric(as.character(p.values$pVal))
  ##add correction for multiple testing
  if(corr=="None"){
    wh <- which(pVal<=0.05)
    p.values.adj$pVal_adj<-pVal}
  if(corr=="Bon"){ wh <- which(pVal<=0.05/nrow(p.values))
  p.values.adj$pVal_adj<-pVal*nrow(p.values)}
  if(corr=="BH") {wh <- which(p.adjust(pVal,"BH")<=0.05)
  p.values.adj$pVal_adj<-p.adjust(pVal,"BH")}
  if(corr=="BY") {wh <- which(p.adjust(pVal,"BY")<=0.05)
  p.values.adj$pVal_adj<-p.adjust(pVal,"BY")}
  if(!corr%in%c("Bon","BH","BY","","None")) stop("Please specify a known correction method for
                                                 multiple testing")
  wh <- p.values$var[wh]
  a<-list(wh,p.values.adj)
  names(a)<-c("selected","pval")
  return(a)
}

###LASSO
lasso <-
  function(data_Xs_in,
           data_Y_in,
           name_Y,
           data_covar_in = NULL,
           covar = character(0)) {
    if (length(covar) > 0) {
      data_covar <-
        data_covar_in[rownames(data_covar_in) %in% rownames(data_Y_in) &
                        rownames(data_covar_in) %in% rownames(data_Xs_in), 
                      colnames(data_covar_in) %in%
                        covar, drop = FALSE]
      data_Y <-
        data_Y_in[rownames(data_Y_in) %in% rownames(data_covar) &
                    rownames(data_Y_in) %in% rownames(data_Xs_in), 
                  colnames(data_Y_in) == name_Y, drop =
                    FALSE]
      data_Xs <-
        data_Xs_in[rownames(data_Xs_in) %in% rownames(data_covar) &
                     rownames(data_Xs_in) %in% rownames(data_Y_in), ,
                   drop = FALSE]
      data_covar <- data_covar[rownames(data_Y), ]
      data_Xs <- data_Xs[rownames(data_Y), ]
      data_Y <- getresiduals_2df(data_Y, data_covar, name_Y, covar)
    } else{
      data_Y <-
        data_Y_in[rownames(data_Y_in) %in% rownames(data_Xs_in), 
                  colnames(data_Y_in) ==
                    name_Y, drop = FALSE]
      data_Xs <-
        data_Xs_in[rownames(data_Xs_in) %in% rownames(data_Y_in), , drop = FALSE]
      data_Xs <- data_Xs[rownames(data_Y), ]
    }
    data_Y <- data.matrix(data_Y)
    data_Xs <- data.matrix(data_Xs)
    model.enet <- cv.glmnet(data_Xs, data_Y, family = "gaussian",
                            alpha = 1)
    cvfit <- model.enet
    
  ##Calcul Y_predit
  Y_predit<-predict(cvfit,newx=data_Xs, s = "lambda.min")
  Y_predit<-Y_predit[rownames(Y_predit),]
  
  ##liste des CPG selectionnés
  tmp_coeffs <- coef(cvfit, s ="lambda.min")
  cg_select<-data.frame(name = tmp_coeffs@Dimnames[[1]][tmp_coeffs@i + 1], coefficient = tmp_coeffs@x)
  
  cg_select<-cg_select$name[cg_select$name!="(Intercept)"]
  a<-list()
  if (length(cg_select)!=0){
    a<-list("selected"=cg_select,"prediction"=Y_predit)  
  }else{
    a<-list("selected"=character(),"prediction"="no_prediction")
  }
  return(a)
  
}


##DSA
DSAreg <- function(Exp,resp, family = gaussian,maxsize = 15, maxsumofpow = 2,
                   maxorderint = 2){
  Exp <- data.frame(cbind(data.frame(Exp), resp=resp))
  res <- DSA(resp ~ 1, data = Exp, family = family, maxsize = maxsize, maxsumofpow
             = maxsumofpow, maxorderint = maxorderint ,nsplits=1,usersplits = NULL)
  form <- gsub("I[(]","",colnames(coefficients(res)))
  form <- gsub("[*]",":",gsub("[)]","",gsub("[:^:]1","",form)))
  if(length(grep(":",form))>0){
    nam <- strsplit(form[grep("[:]",form)],":")
    for(j in 1:length(nam)){
      nam[[j]] <- gsub("[[:space:]]","",nam[[j]])
      name <- nam[[j]][1]
      for(k in 2:length(nam[[j]]))
        name <- paste(name,":",nam[[j]][k],sep="")
      Exp <- cbind(Exp,name=apply(Exp[,nam[[j]]],1,prod))
    }}
  form2 <- "resp~1"
  if(length(form)>1)for(i in 2:length(form)) form2 <- paste(form2,"+",form[i])
  res2 <- lm(form2, data=data.frame(Exp))
  ##decomment next line and change "prediction" to pred in the return line
  ##if outcome predicted by DSA is needed (not used presently)
  #pred <- predict(res2,Exp)
  coef <- summary(res2)$coefficients
  coef <- as.character(rownames(coef)[rownames(coef)!="Intercept"])
  
  return(list(selected=coef[coef!="(Intercept)"], pred="prediction"))
}


##################################################################
####3. defining some functions used to assess methods performance
#################################################################

sensitivity<-function(truepred, predfound){
  return(length(truepred[truepred%in%predfound])/length(truepred))
}
fdp<-function(truepred, predfound){ ##false discovery proportion
  if (length(predfound)==0) {return(0)
  }else{
    return(length(predfound[!predfound%in%truepred])/length(predfound))}
}
specificity<-function(truepred, predfound,n_base){
  return(
    (n_base-length(truepred)-length(predfound[!predfound%in%truepred]))/(n_base-length(truepred)))
}


#################################################################
##4. defining the simulation function which will be parallelized
#################################################################
##it first generates datasets, then applies methods and then assessed 
##their performance
f0<-function(x){
  ##important: the parallelization is done on the seed 
  set.seed(x)
  
  ##generating datasets
  simu <-
    simulator(
      E_true = dataExp_true,
      M_true = M1_true,
      R2_tot = R2_fixed,
      propmE = 0,
      propmEY = 0.1,
      propmY = 0,
      BetamEY = BetamEY,
      BetamY = 0,
      n_mE = 0,
      n_mEY = n_mEY,
      n_mY = 0,
      n_EmE = 0,
      n_EmEY = n_EmEY,
      n_Ey = n_Ey,
      n_EmE_U_n_EmEY = 0,
      n_Ey_U_n_EmE = 0,
      n_Ey_U_n_EmEY = n_Ey_U_n_EmEY,
      n_Ey_U_n_EmE_U_n_EmEY = 0,
      BetaEmE = 0,
      BetaEmEY = BetaEmEY,
      BetaEy = BetaEy,
      test_and_training = TRUE
    )
  ###################
  ###applying methods
  ####################
  ##ExWas on the intermediate layer
  predBMI_M <-
    list(ewas_BH = ewas(
      as.data.frame(simu$M_train),
      as.data.frame(simu$Y_train),
      colnames(as.data.frame(simu$Y_train)),
      corr = "BH"
    ))
  ##ExWAS on E
  predBMI_E <-
    list(
      ewas_BH = ewas(
        as.data.frame(simu$E_train),
        as.data.frame(simu$Y_train),
        colnames(as.data.frame(simu$Y_train)),
        corr = "BH"
      ),
      ewas_Bon = ewas(
        as.data.frame(simu$E_train),
        as.data.frame(simu$Y_train),
        colnames(as.data.frame(simu$Y_train)),
        corr = "Bon"
      )
    )
    print("ewas done")
  
  #########################
    #oMITM
    if (length(predBMI_M$ewas_BH$selected) != 0) {
      select_M <-
        as.data.frame(simu$M_train[,
                                   colnames(simu$M_train) %in% predBMI_M$ewas_BH$selected,
                                   drop =  FALSE])
      print(ncol(select_M))
      rownames(select_M) <- rownames(simu$M_train)
      colnames(select_M) <- predBMI_M$ewas_BH$selected
      list <- list()
      list_exp <- list()
      list_nom <- list()
      list_ewas_signif <- list()
      ##step b using the exwas performed on M as step a
      for (i in (1:ncol(simu$E_train))) {
        predE_select_M <-
          ewas(
            as.data.frame(select_M),
            (simu$E_train[, i, drop = FALSE]),
            colnames(simu$E_train[, i, drop = FALSE]),
            corr = "None",
            data_covar_in = as.data.frame(simu$Y_train),
            covar = colnames(as.data.frame(simu$Y_train)[1])
          )
        
        list <- c(list, list(predE_select_M))
        list_nom <- c(list_nom, list(colnames(simu$E_train)[i]))
        
        
        list_exp <- c(list_exp, list(colnames(simu$E_train)[i]))
        temp_ewas <-
          cbind(predE_select_M$pval, rep(colnames(simu$E_train)[i], 
                                         nrow(predE_select_M$pval)))
        list_ewas_signif <- c(list_ewas_signif, list(temp_ewas))
        remove(temp_ewas)
        remove(predE_select_M)
      }
      df_all_ewas <- do.call("rbind", list_ewas_signif)
      if (!is.null(df_all_ewas)) {
        df_all_ewas$pVal_adj <- p.adjust(df_all_ewas$pVal, "BH")
        colnames(df_all_ewas)[8] <- "exposures"
        names(list) <- as.vector(unlist(list_nom))
      }
      exp <- df_all_ewas$exposures[df_all_ewas$pVal_adj <= 0.05]
      n_exp_select <- length(unique(exp))
      ##step c
      if (length(exp) != 0) {
        select_E <- simu$E_train[, colnames(simu$E_train) %in% exp, drop = FALSE]
        ##ExWAS implementation for step c
        predBMI_E_MITM <-
          ewas(
            as.data.frame(select_E),
            as.data.frame(simu$Y_train),
            colnames(as.data.frame(simu$Y_train)),
            corr = "BH"
          )
        ##DSA implementation for step c
        predBMI_E_MITMdsa <-
          DSAreg(
            Exp = as.data.frame(select_E),
            resp = simu$Y_train,
            maxsize = floor(ncol(simu$E_train) / 10),
            maxsumofpow = 1,
            maxorderint = 1
          )
        predReducedExp <- list(selected = unique(exp), pred = "NULL")
      } else{
        predReducedExp <- list(vector(), vector())
        names(predReducedExp) <- c("selected", "pred")
      }
      if (exists("predBMI_E_MITM")) {
        
      } else{
        predBMI_E_MITM <- list(vector(), vector())
        names(predBMI_E_MITM) <- c("selected", "pval")
      }
      if (exists("predBMI_E_MITMdsa")) {
        
      } else{
        predBMI_E_MITMdsa <- list(vector(), vector())
        names(predBMI_E_MITM) <- c("selected", "pred")
      }
    } else{
      predBMI_E_MITM <- list(vector(), vector())
      names(predBMI_E_MITM) <- c("selected", "pval")
      n_exp_select = 0
      predReducedExp <- list(vector(), vector())
      names(predReducedExp) <- c("selected", "pred")
      predBMI_E_MITMdsa <- list(vector(), vector())
      names(predBMI_E_MITMdsa) <- c("selected", "pred"
      )
    
  }
    ##storing results in a list
    predBMI_E <-
      c(
        predBMI_E,
        MITM = list(predBMI_E_MITM),
        MITMdsa = list(predBMI_E_MITMdsa),
        ReducedExp = list(predReducedExp)
      )
    print("oMITM done")
    
    ###Control method : random sampling on a random set of exposures of same
    ##dimension as the reduced exposome of oMITM
    if (n_exp_select > 0) {
      tirage <-
        ewas(
          as.data.frame(simu$E_train)[, sample(colnames(as.data.frame(simu$E_train)),
                                               n_exp_select), drop =
                                        FALSE],
          as.data.frame(simu$Y_train),
          colnames(as.data.frame(simu$Y_train)),
          corr = "BH"
        )
    } else{
      tirage <- list(selected = character(0), null = "null")
    }
    ##storing results in the same list
    predBMI_E <- c(predBMI_E, random_sampling = list(tirage))
    print(n_exp_select)
    
    ###################################
    ##mediation
    if (length(predBMI_E$ewas_BH$selected) != 0) {
      select_E <-
        as.data.frame(simu$E_train[,
                                   colnames(simu$E_train) %in% predBMI_E$ewas_BH$selected, 
                                   drop = FALSE])
      rownames(select_E) <- rownames(simu$E_train)
      colnames(select_E) <- predBMI_E$ewas_BH$selected
      #step a
      list_temp_ewas_med <- list()
      for (i in 1:ncol(simu$M_train)) {
        exp_affecting_M_all <-
          ewas(as.data.frame(select_E),
               (simu$M_train[, i, drop = FALSE]),
               colnames(simu$M_train[, i, drop = FALSE]),
               corr = "None")
        temp_ewas_med <-
          cbind(exp_affecting_M_all$pval, rep(colnames(simu$M_train)[i], 
                                              nrow(exp_affecting_M_all$pval)))
        list_temp_ewas_med <- c(list_temp_ewas_med, list(temp_ewas_med))
      }
      ewas_med <- do.call("rbind", list_temp_ewas_med)
      if (!is.null(ewas_med)) {
        ewas_med$pVal_adj_1 <- p.adjust(ewas_med$pVal, "BH")
        colnames(ewas_med)[8] <- "cpg"
        colnames(ewas_med)[1] <- "exp"
      }
      #step b
      list_temp_ewas_med_2 <- list()
      for (i in 1:ncol(select_E)) {
        M_affecting_Y_all <-
          ewas(
            as.data.frame(simu$M_train),
            as.data.frame(simu$Y_train),
            colnames(as.data.frame(simu$Y_train)),
            corr = "None",
            data_covar_in = (select_E[, i, drop = FALSE]),
            covar = colnames(select_E[, i, drop = FALSE])
          )
        temp_ewas_med_2 <-
          cbind(M_affecting_Y_all$pval, rep(colnames(select_E)[i], 
                                            nrow(M_affecting_Y_all$pval)))
        list_temp_ewas_med_2 <-
          c(list_temp_ewas_med_2, list(temp_ewas_med_2))
      }
      ewas_med_2 <- do.call("rbind", list_temp_ewas_med_2)
      if (!is.null(ewas_med_2)) {
        ewas_med_2$pVal_adj_2 <- p.adjust(ewas_med_2$pVal, "BH")
        colnames(ewas_med_2)[8] <- "exp"
        colnames(ewas_med_2)[1] <- "cpg"
      }
      
      ewas_med_tot <-
        merge(
          ewas_med,
          ewas_med_2,
          by.x = c("exp", "cpg"),
          by.y = c("exp", "cpg")
        )
      exp_med <-
        unique(ewas_med_tot$exp[ewas_med_tot$pVal_adj_2 <= 0.05 &
                                  ewas_med_tot$pVal_adj_1 <= 0.05])
      exp_med <- exp_med[!is.na(exp_med)]
      if (length(exp_med) != 0) {
        predMediation <- list(selected = exp_med, pred = NULL)
      } else{
        predMediation <- list(selected = vector(), pred = vector())
      }
      
    } else{
      predMediation <- list(selected = vector(), pred = vector())
    }
    predBMI_E <- c(predBMI_E, mediation = list(predMediation))
  
 #######################################
  ##applying agnostic methods
  ## lasso
  predlasso <-
    lasso(
      data_Xs_in = as.data.frame(simu$E_train),
      data_Y_in = as.data.frame(simu$Y_train),
      colnames(as.data.frame(simu$Y_train))
    )
  predBMI_E <- c(predBMI_E, lasso_CV = list(predlasso))
  print("lasso_ done")
  ##DSA
  predDSA <-
    DSAreg(
      Exp = simu$E_train,
      resp = simu$Y_train,
      maxsize = floor(ncol(simu$E_train) / 10),
      maxsumofpow = 1,
      maxorderint = 1
    )
  predBMI_E <- c(predBMI_E, DSA = list(predDSA))
  print("DSA done")
  
  ##########################################
  ##assessing performance for ExWAS on M
  truepred<-simu$yM_M$predictors
  for (k1 in 1:length(predBMI_M)){
    truepred<-simu$yM_M$predictors
    predfound<-predBMI_M[[k1]]$selected
    if (exists("predfound")&exists("truepred")){
      if (length(predfound)==0) {print("no predictors found")}
      a<-sensitivity(truepred,predfound)
      b<-specificity(truepred,predfound,ncol(simu$M_train))
      c<-fdp(truepred,predfound)
      d<-estimatedR2(simu$M_test,predfound,simu$Y_test)$r.squared
      # print(a)
      # print(b)
      # print(c)
      # print(d)
      remove(predfound,truepred)
    }else{
      if(!exists("predfound")&exists("truepred")){
        a<-0
        b<-specificity(truepred,numeric(0),ncol(simu$M_train))
        c<-0
        d<-0
        remove(truepred)
      }else{
      if(!exists("truepred")&exists("predfound")){
        a<-1
        b<-specificity(numeric(0),predfound,ncol(simu$M_train))
        c<-fdp(numeric(0),predfound)
        d<-estimatedR2(simu$M_test,predfound,simu$Y_test)$r.squared
        remove(predfound)
      }else{
      if(!exists("truepred")&!exists("predfound")){
        a<-1
        b<-1
        c<-0
        d<-0
      }}}
    }
      predBMI_M[[k1]]<-c(predBMI_M[[k1]],sens=a,spec=b,fdp=c,R2_test=d)
      remove(a)
      remove(b)
      remove(c)
      remove(d)

    }
  
  ##assessing performance for all methods on E
  truepred<-simu$y_E$predictors
  for (k1 in 1:length(predBMI_E)){
    truepred<-simu$y_E$predictors
    truepred<-simu$y_E$predictors
    predfound<-predBMI_E[[k1]]$selected
    if (exists("predfound")&exists("truepred")){
      if (length(predfound)==0) {print("no predictors found")}
      a<-sensitivity(truepred,predfound)
      b<-specificity(truepred,predfound,ncol(simu$E_train))
      c<-fdp(truepred,predfound)
      d<-estimatedR2(simu$E_test,predfound,simu$Y_test)$r.squared
      # print(a)
      # print(b)
      # print(c)
      # print(d)
      remove(predfound,truepred)
    }else{
      if(!exists("predfound")&exists("truepred")){
        a<-0
        b<-specificity(truepred,numeric(0),ncol(simu$E_train))
        c<-0
        d<-0
        remove(truepred)
      }else{
      if(!exists("truepred")&exists("predfound")){
        a<-1
        b<-specificity(numeric(0),predfound,ncol(simu$E_train))
        c<-fdp(numeric(0),predfound)
        d<-estimatedR2(simu$E_test,predfound,simu$Y_test)$r.squared
        remove(predfound)
      }else{
      if(!exists("truepred")&!exists("predfound")){
        a<-1
        b<-1
        c<-0
        d<-0
      }}}
    }
    predBMI_E[[k1]]<-c(predBMI_E[[k1]],sens=a,spec=b,fdp=c,R2_test=d)
    remove(a)
    remove(b)
    remove
    remove(d)
  }
  
  ##assessing the performance considering only the predictors having an indirect
  ##effect
  truepred<-simu$yM_E$predictors
  for (k1 in 1:length(predBMI_E)){
    truepred<-simu$yM_E$predictors
    predfound<-predBMI_E[[k1]]$selected
    if (exists("predfound")&exists("truepred")){
      if (length(predfound)==0) {print("no predictors found")}
      a<-sensitivity(truepred,predfound)
      b<-specificity(truepred,predfound,ncol(simu$E_train))
      c<-fdp(truepred,predfound)
      d<-estimatedR2(simu$E_test,predfound,simu$Y_test)$r.squared
      # print(a)
      # print(b)
      # print(c)
      # print(d)

      remove(predfound,truepred)
    }else{
      if(!exists("predfound")&exists("truepred")){
        a<-0
        b<-specificity(truepred,numeric(0),ncol(simu$E_train))
        c<-0
        d<-0
        remove(truepred)
      }else{
      if(!exists("truepred")&exists("predfound")){
        a<-1
        b<-specificity(numeric(0),predfound,ncol(simu$E_train))
        c<-fdp(numeric(0),predfound)
        d<-estimatedR2(simu$E_test,predfound,simu$Y_test)$r.squared
        remove(predfound)
      }else{
      if(!exists("truepred")&!exists("predfound")){
        a<-1
        b<-1
        c<-0
        d<-0
      }}}
    }
 predBMI_E[[k1]]<-c(predBMI_E[[k1]],sens=a,spec=b,fdp=c,R2_test=d)
  remove(a)
  remove(b)
  remove
  remove(d)
}
  print("performance characterized")
  ##building the list with datasets generated + results of methods + 
  ##performance to return
  A<-list(simu=simu,predBMI_E=predBMI_E,predBMI_M=predBMI_M,
          nl_exp_select=n_exp_select)
  remove(simu)
  remove(select_M)
  remove(predBMI_M)
  remove(predBMI_E)
  remove(predBMI_E_MITM_WM)
  remove(predBMI_E_MITM)
  remove(predlasso)
  remove(n_exp_select)
  gc()
  return(A)
  
}

#################################
##5. Running simulations
#################################


##loading real datasets
dataExp_true <- readRDS("20190205 Exposome simu borne.rds")
M1_true <- readRDS("20191129 Methylome simu.Rds")
M1_true <- scale(M1_true)

##initialization
list_simulated_data <- list()
list_list_predBMI_E <- list()
list_list_predBMI_M <- list()
list_list_nl_exp_select <- list()

##setting simulations parameters
n_iter <- 100 ##number of iterations for one scenarios
##parameters for generating datasets
##all combinations will be tested (each combination allows to build a scenario)
##(adapt the code of the loop if multiple values instead of single values for 
##some parameters)
c_n_my <- c(10, 18, 25, 100)
c_n_R2_fixed <- c(0.01, 0.05, 0.1, 0.4)
BetamEY = 0.01
c_BetaEy <- c(0.0001, 0.001, 0.01, 0.1, 0.5)
c_n_Ey <- c(1, 3, 10, 25)
c_BetaEmEY <- c(0.0001, 0.001, 0.01, 0.1, 0.5)
n_mE <- 0
n_mEY <- 0

##initialization of table of results
comp_method <-
  data.frame(
    Methods = vector(),
    Association_tested = vector(),
    Nb_true_predictors_of_BMI_in_M = numeric(0),
    Nb_true_predictors_of_BMI_in_E = numeric(0),
    Total_variability_of_BMI_explained_by_EandM = numeric(0),
    Total_variability_of_BMI_explained_by_E = numeric(0),
    Total_variability_of_BMI_explained_by_M = numeric(0),
    Mean_variability_of_M_explained_by_E_for_Mey = numeric(0),
    Number_iterations = numeric(0),
    Mean_number_exp_selected_to_be_randomly_tested =
      numeric(0),
    Mean_number_predictors_found = numeric(0),
    Mean_sensitivity = numeric(0),
    Mean_specificity = numeric(0),
    Mean_fdp = numeric(0),
    Mean_R2_test = numeric(0),
    Mean_mediated_sensitivity = numeric(0),
    Mean_mediated_specificity = numeric(0),
    Mean_mediated_fdp = numeric(0),
    Mean_R2_test = numeric(0),
    SD_number_predictors_found = numeric(0),
    SD_sensitivity = numeric(0),
    SD_specificity = numeric(0),
    SD_fdp = numeric(0),
    SD_R2_test = numeric(0),
    SD_mediated_sensitivity = numeric(0),
    SD_mediated_specificity = numeric(0),
    SD_mediated_fdp = numeric(0),
    SD_R2_test = numeric(0),
    Which_iteration = numeric(0)
  )


##looping on the different vectors of parameters
n = 1
for (i2 in 1:length(c_n_R2_fixed)) {
  R2_fixed <- c_n_R2_fixed[i2]
  for (i3 in 1:length(c_BetaEy)) {
    BetaEy <- c_BetaEy[i3]
    for (i4 in 1:length(c_n_Ey)) {
      n_Ey <- c_n_Ey[i4]
      n_Ey_U_n_EmEY <- n_Ey
      n_EmEY <- n_Ey
      for (i1 in 1:length(c_n_my)) {
        n_mEY <- c_n_my[i1]
        for (i5 in 1:length(c_BetaEmEY)) {
          BetaEmEY <- c_BetaEmEY[i5]
          print(n)
          iteration_OK <- TRUE
          
          if (n < 1501) {
            n <- n + 1
          } else{
            if (n_mEY != 0) {
              n_pat_mEY <- (n_mEY / n_EmEY)
              if (trunc(n_pat_mEY) != n_pat_mEY) {
                iteration_OK <- FALSE
              }
            }
            if (n_mE != 0) {
              n_pat_mE <- (n_mE / n_EmE)
              if (trunc(n_pat_mE) != n_pat_mE) {
                iteration_OK <- FALSE
              }
            }
            if (iteration_OK == TRUE) {
              n_row = nrow(comp_method)
              simulated_data <- list()
              list_predBMI_E <- list()
              list_predBMI_M <- list()
              list_nl_exp_select <- list()
              start_time <- Sys.time()
              ##parallelization of f0
              cl <- makeCluster(getOption("cl.cores", round(detectCores())))
              clusterExport(
                cl,
                list(
                  "simulator",
                  "simResponseSimple",
                  "estimatedR2",
                  "getresiduals_2df",
                  "ewas",
                  "lasso",
                  "lasso_stab",
                  "DSAreg",
                  "sensitivity",
                  "fdp",
                  "specificity",
                  "f0",
                  "dataExp_true",
                  "M1_true",
                  "R2_fixed",
                  "n_Ey",
                  "n_mEY",
                  "BetaEy",
                  "BetamEY",
                  "n_Ey_U_n_EmEY",
                  "n_EmEY",
                  "BetaEmEY"
                )
              )
              clusterEvalQ(cl, list(library("boot"), library("reshape"),
                                    library("glmnet"), library("DSA")))
              results_1_jeu <- clusterApply(cl, 1:n_iter, f0)
              stopCluster(cl)
              simulated_data <- lapply(results_1_jeu, function(x)
                x$simu)
              
              ##structure of results priorized by methods and not anymore 
              ##priorized by datasets
              list_predBMI_E <- lapply(results_1_jeu, function(x)
                x$predBMI_E)
              list_predBMI_M <- lapply(results_1_jeu, function(x)
                x$predBMI_M)
              list_nl_exp_select <-
                lapply(results_1_jeu, function(x)
                  x$nl_exp_select)
              remove(results_1_jeu)
              
              
              ###compilation of results for this scenario
              
              ##table describing the empirical characteristics of 
              ##the simulated datasets
              param_simu <-
                data.frame(
                  Parameters = vector(),
                  Fixed_or_measured = vector(),
                  Value = numeric(0)
                )
              param_simu[1, ] <-
                c(
                  as.character("Nb_true predictors of BMI in M"),
                  as.character("Fixed"),
                  mean(unlist(
                    lapply(simulated_data, function(X)
                      length(X[[10]]$beta))
                  ))
                )
              param_simu[2, ] <-
                c("Nb_true predictors of BMI in E",
                  "Measured",
                  mean(unlist(
                    lapply(simulated_data, function(X)
                      length(X[[8]]$beta))
                  )))
              param_simu[3, ] <-
                c("Total variability of BMI explained by (E+M)",
                  "Fixed",
                  R2_fixed)
              param_simu[4, ] <-
                c("Mean variability of BMI explained by E",
                  "Measured",
                  mean(unlist(
                    lapply(simulated_data, function(X)
                      (X[[11]]$BMI_all_exp))
                  )))
              param_simu[5, ] <-
                c("Mean variability of BMI explained by M",
                  "Measured",
                  mean(unlist(
                    lapply(simulated_data, function(X)
                      (X[[11]]$BMI_all_M))
                  )))
              param_simu[6, ] <-
                c(
                  "Mean variability of part of M explained by E explained by E",
                  "Measured",
                  mean(unlist(
                    lapply(simulated_data, function(X)
                      (X[[11]]$mean_M_E))
                  ))
                )
              param_simu[7, ] <- c("Number_iterations", "Fixed", n_iter)
              param_simu[8, ] <-
                c(
                  "Mean_number_exp_selected_to_be_randomly_tested",
                  "Measured",
                  mean(unlist(list_nl_exp_select))
                )
              
              
              ##summarizing each method performance by a line in comp_method datadrame within
              ##this scenario
              
              for (k1 in (1:length(list_predBMI_M[[1]]))) {
                comp_method[n_row + k1, ] <-
                  c(
                    names(list_predBMI_M[[1]][k1]),
                    "BMI - M",
                    param_simu[1, 3],
                    param_simu[2, 3],
                    R2_fixed,
                    param_simu[4, 3],
                    param_simu[5, 3],
                    param_simu[6, 3],
                    n_iter,
                    NA,
                    mean(unlist(
                      lapply(list_predBMI_M, function(X)
                        length(X[[k1]][[1]]))
                    )),
                    mean(unlist(
                      lapply(list_predBMI_M, function(X)
                        (X[[k1]][[3]]))
                    ), na.rm = TRUE),
                    mean(unlist(
                      lapply(list_predBMI_M, function(X)
                        (X[[k1]][[4]]))
                    ), na.rm = TRUE),
                    mean(unlist(
                      lapply(list_predBMI_M, function(X)
                        (X[[k1]][[5]]))
                    ), na.rm = TRUE),
                    mean(unlist(
                      lapply(list_predBMI_M, function(X)
                        (X[[k1]][[6]]))
                    ), na.rm = TRUE),
                    NA,
                    NA,
                    NA,
                    NA,
                    sd(unlist(
                      lapply(list_predBMI_M, function(X)
                        length(X[[k1]][[1]]))
                    )),
                    sd(unlist(
                      lapply(list_predBMI_M, function(X)
                        (X[[k1]][[3]]))
                    ), na.rm = TRUE),
                    sd(unlist(
                      lapply(list_predBMI_M, function(X)
                        (X[[k1]][[4]]))
                    ), na.rm = TRUE),
                    sd(unlist(
                      lapply(list_predBMI_M, function(X)
                        (X[[k1]][[5]]))
                    ), na.rm = TRUE),
                    sd(unlist(
                      lapply(list_predBMI_M, function(X)
                        (X[[k1]][[6]]))
                    ), na.rm = TRUE),
                    NA,
                    NA,
                    NA,
                    NA,
                    n
                  )
              }
              
              for (k2 in (1:length(list_predBMI_E[[1]]))) {
                comp_method[n_row + k2 + k1, ] <-
                  c(
                    names(list_predBMI_E[[1]][k2]),
                    "BMI - E",
                    param_simu[1, 3],
                    param_simu[2, 3],
                    R2_fixed,
                    param_simu[4, 3],
                    param_simu[5, 3],
                    param_simu[6, 3],
                    n_iter,
                    param_simu[8, 3],
                    mean(unlist(
                      lapply(list_predBMI_E, function(X)
                        length(X[[k2]][[1]]))
                    )),
                    mean(unlist(
                      lapply(list_predBMI_E, function(X)
                        (X[[k2]][[3]]))
                    )),
                    mean(unlist(
                      lapply(list_predBMI_E, function(X)
                        (X[[k2]][[4]]))
                    )),
                    mean(unlist(
                      lapply(list_predBMI_E, function(X)
                        (X[[k2]][[5]]))
                    )),
                    mean(unlist(
                      lapply(list_predBMI_E, function(X)
                        (X[[k2]][[6]]))
                    )),
                    
                    mean(unlist(
                      lapply(list_predBMI_E, function(X)
                        (X[[k2]][[7]]))
                    )),
                    mean(unlist(
                      lapply(list_predBMI_E, function(X)
                        (X[[k2]][[8]]))
                    )),
                    mean(unlist(
                      lapply(list_predBMI_E, function(X)
                        (X[[k2]][[9]]))
                    )),
                    mean(unlist(
                      lapply(list_predBMI_E, function(X)
                        (X[[k2]][[10]]))
                    )),
                    
                    
                    sd(unlist(
                      lapply(list_predBMI_E, function(X)
                        length(X[[k2]][[1]]))
                    )),
                    sd(unlist(
                      lapply(list_predBMI_E, function(X)
                        (X[[k2]][[3]]))
                    )),
                    sd(unlist(
                      lapply(list_predBMI_E, function(X)
                        (X[[k2]][[4]]))
                    )),
                    sd(unlist(
                      lapply(list_predBMI_E, function(X)
                        (X[[k2]][[5]]))
                    )),
                    sd(unlist(
                      lapply(list_predBMI_E, function(X)
                        (X[[k2]][[6]]))
                    )),
                    sd(unlist(
                      lapply(list_predBMI_E, function(X)
                        (X[[k2]][[7]]))
                    )),
                    sd(unlist(
                      lapply(list_predBMI_E, function(X)
                        (X[[k2]][[8]]))
                    )),
                    sd(unlist(
                      lapply(list_predBMI_E, function(X)
                        (X[[k2]][[9]]))
                    )),
                    sd(unlist(
                      lapply(list_predBMI_E, function(X)
                        (X[[k2]][[10]]))
                    )),
                    n
                  )
              }
              
              
              ##storing generated datasets and methods results for this scenario
              list_list_predBMI_E <- c(list_list_predBMI_E, list(list_predBMI_E))
              list_list_predBMI_M <- c(list_list_predBMI_M, list(list_predBMI_M))
              list_list_nl_exp_select <-
                c(list_list_nl_exp_select,
                  list(list_nl_exp_select))
              end_time <- Sys.time()
              end_time - start_time
              ##saving
              saveRDS(comp_method,
                      "comp_method_mediation_and_direct.Rds")
              saveRDS(
                simulated_data,
                file = paste(
                  n,
                  '_simulated_data_scenario_iteration_mediation_and_direct.Rds'
                )
              )
              saveRDS(
                list_list_predBMI_E,
                "list_list_predBMI_E_mediation_and_direct.rds"
              )
              saveRDS(
                list_list_predBMI_M,
                "list_list_predBMI_M_mediation_and_direct.rds"
              )
              saveRDS(
                list_list_nl_exp_select,
                "list_list_nl_exp_select_mediation_and_direct.rds"
              )
              
              remove(simulated_data)
              remove(list_predBMI_E)
              remove(list_predBMI_M)
              remove(list_nl_exp_select)
            }
            n = n + 1
          }
        }
      }
    }
  }
}


