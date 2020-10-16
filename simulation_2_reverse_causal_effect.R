##this code allows to perform a simulation to assess performance in terms of
## sensitivity and specificity of prespecified statistical methods
##used to find the true predictors of an health among the exposome. The causal
##structures considered involve a reverse causative  link from the outcome on
##the exposome Some of them use an intermediary layer, some not.
##it contains 5 ##parts:

##1. defining the functions allowing to generate a realistic dataset of exposome
##intermediate layer and outcome. The three layers(E, M and Y) can be linearly
##related to simulate various causal structures.

#It needs real datasets (exposome/intermediate layer/outcome) as inputs, as well
##as parameters allowing to define the association within the three layers (number of
#predictors, #variability explained, correlation..)

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
  function(E_true,
           ##real exposome
           M_true,
           ## real intermediate layer, eg methylome
           Y_true,
           ##real outcome
           n_mY = 0,
           #variables of M not affected by E
           ##but affecting Y.
           R2_mY = 0,
           ##variability of M not affected by E affecting Y
           BetamY = 0.1,
           ##corresponding effect coefficient
           n_yM = 50,
           ##number of variables of M affected by Y
           Beta_yM = 0.001,
           ##corresponding effect coefficient
           n_yE = 5,
           ##number of variables of E affected by Y
           Beta_yE = 0.01,
           ##corresponding effect coefficient
           test_and_training = TRUE)#generating only a training set or
           ##alternatively also a test set of same size) 
           {
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
           ##sampling with replacement the real outcome data
           data.X <- as.data.frame(Y_true)
           names_row <- rownames(data.X)
           data.X <-
             data.X[sample(1:nrow(data.X), 2 * nrow(data.X), replace = TRUE),
                    , drop =
                      FALSE]
           rownames(data.X) <- c(names_row, sprintf('boot%s', names_row))
           Y_boot <- data.X
           remove(data.X)
           
           ##setting if necessary linear relationship from M to Y
           if (n_mY != 0) {
             if (R2_mY == 0) {
               stop("error: n_MY and R2_mY non consistent")
             }
             ##generating vector of effect coefficients
             if (length(BetamY) == 1) {
               Betapred_mY <- rep(BetamY, n_mY)
             } else{
               if (length(BetamY) != n_mY) {
                 stop("error: Betas for M explaining Y are not consistent
               with the number of predictors")
               }
               Betapred_mY <- BetamY
             }
             ##random sampling of variables of M affecting Y
             ind_mY <- sample(ncol(M1), n_mY)
             ##generating yM: part of the outcome which is a
             ##linear combination of variables of M not affected by E
             yM <-
               simResponseSimple(
                 met = M1,
                 Nmet = length(ind_mY),
                 beta = Betapred_mY,
                 cpg = ind_mY
               )
             if (!is.na(R2_mY)) {
               if (((R2_mY) != 0)) {
                 sigma <- var(yM$resp) * (1 / R2_mY - 1)
               } else{
                 warning("R2 was not specified, automatic value")
                 R2_mY = 0.00000001
                 sigma <- var(Y$resp) * (1 / R2_mY - 1)
               }
               Y <-
                 as.matrix(Y$resp + rnorm(length(Y$resp), mean(Y$resp),
                                          sqrt(sigma)), ncol = 1)
               #standardiZation
               Y <- as.data.frame(scale(Y))
               
             }
             ##empirical estimation of R2 ( variability of Y explained by M)
             R2_mY_measured <- estimatedR2(M1, ind_mY, Y)$r.squared
           } else{
             ##if no effect, creating an empty yM
             yM <- list(resp = Y_boot,
                        beta = NULL,
                        predictors = NULL)
             Y <- as.data.frame(Y_boot)
             R2_mY_measured <- 0
             ind_mY <- integer()
             Betapred_mY = NULL
           }
           
           ##setting if necessary linear relationship from Y to M
           if (n_yM != 0) {
             ##random sampling of variables of M affected by Y
             if (n_mY == 0) {
               ind_yM_M <- sample(ncol(M1), n_yM)
             } else{
               ind_yM_M <- sample((1:ncol(M1))[ind_mY], n_yM)
             }
             ##generating vector of effect coefficients
             if (length(BetayM) == 1) {
               Betapred_yM_M <- rep(BetayM, n_yM)
             } else{
               if (length(BetayM) != n_yM) {
                 stop(
                   "error: Betas for Y explaining M  are not consistent with
                   the number of predicted cpgs"
                 )
               }
               Betapred_yM_M <- BetayM
             }
             ##adding a linear effect of Y on selected variables of M
             list_R2_M <- list()
             for (i in 1:n_yM) {
               M1[, ind_yM_M[i]] <- as.numeric(M1[, ind_yM_M[i]] +
                                                 simResponseSimple(
                                                   met = Y,
                                                   Nmet = 1,
                                                   beta = Betapred_yM_M[i],
                                                   cpg = 1
                                                 )$resp)
               
               list_R2_M <-
                 c(list_R2_M, list(estimatedR2(Y, colnames(Y)[1], M1[, 
                                                  ind_yM_M[i], drop = FALSE])))
             }
             ##empirical estimation of mean R2 (mean variability of M affected by Y)
             mean_R2_M <-
               mean(unlist(lapply(list_R2_M, function(x)
                 x$r.squared)), na.rm = T)
             SD_R2_M <-
               sd(unlist(lapply(list_R2_M, function(x)
                 x$r.squared)), na.rm = T)
           } else{
             ind_yM_M <- integer()
             Betapred_yM_M = NULL
             list_R2_M = NULL
             mean_R2_M = 0
           }
           #setting if necessary linear relationship from Y to E
           list_R2_E <- list()
           if (n_yE != 0) {
             ##random sampling of variables of e affected by Y
             ind_yE_E <- sample(ncol(dataExp), n_yE)
             
             ##generating vector of effect coefficients
             if (length(BetayE) == 1) {
               Betapred_yE_E <- rep(BetayE, n_yE)
             } else{
               if (length(BetayE) != n_yE) {
                 stop(
                   "error: Betas for Y explaining E  are not consistent with the
                   number of predicted exposures"
                 )
               }
               Betapred_yE_E <- BetayE
             }
             ##adding a linear effect of Y on selected variables of e
             for (i in 1:n_yE) {
               dataExp[, ind_yE_E] <- as.numeric(dataExp[, ind_yE_E[i]] +
                                                   simResponseSimple(
                                                     met = Y,
                                                     Nmet = 1,
                                                     beta = Betapred_yE_E[i],
                                                     cpg = 1
                                                   )$resp)
               
               list_R2_E <-
                 c(list_R2_E, list(estimatedR2(Y, colnames(Y)[1], dataExp[,
                                                ind_yE_E[i], drop = FALSE])))
             }
             ##empirical estimation of mean R2 (mean variability of e affected by Y)
             mean_R2_E <-
               mean(unlist(lapply(list_R2_E, function(x)
                 x$r.squared)), na.rm = T)
             SD_R2_E <-
               sd(unlist(lapply(list_R2_E, function(x)
                 x$r.squared)), na.rm = T)
           } else{
             ind_yE_E <- integer()
             Betapred_yE_E = NULL
             list_R2_E = NULL
             mean_R2_E = 0
           }
           
           
           ##Building a result object with generated datesets;
           ##vector of predictors and effects
           results <-
             
             list(
               Y_train = Y[1:(nrow(M1) / 2), , drop = FALSE],
               E_train = dataExp[1:(nrow(M1) / 2), , drop = FALSE],
               M_train = M1[1:(nrow(M1) / 2), , drop = FALSE],
               Y_test = Y[(nrow(M1) / 2):nrow(M1), , drop = FALSE],
               E_test = dataExp[(nrow(M1) / 2):nrow(M1), , drop = FALSE],
               M_test = M1[(nrow(M1) / 2):nrow(M1), , drop = FALSE],
               cpg_predictors = list(
                 name = colnames(M1)[ind_mY],
                 indices = ind_mY,
                 betas = Betapred_mY,
                 R2 = NULL
               ),
               cpg_predicted = list(
                 name = colnames(M1)[ind_yM_M],
                 indices = ind_yM_M,
                 betas = Betapred_yM_M,
                 R2 = list_R2_M
               ),
               exp_predicted = list(
                 name = colnames(dataExp)[ind_yE_E],
                 indices = ind_yE_E,
                 betas = Betapred_yE_E,
                 R2 = list_R2_E
               ),
               R2_mY_true = R2_mY,
               R2_mY_measured = R2_mY_measured,
               R2_yM_mean = mean_R2_M,
               R2_yM_SD = SD_R2_M,
               R2_yE_mean = mean_R2_E,
               R2_yE_SD = SD_R2_E
             )
           
           return(results)
           }


##function used to create a linear response
simResponseSimple <- function(met,
                              ##matrix of potential predictors
                              Nmet = NA,
                              ##number of predictors
                              beta = NULL,
                              ##vector of effects
                              cpg = NULL) {
  ##optionnal: directly specifying
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
  rownames(mean) <- rownames(met)
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
getresiduals_2df <- function(data_Y_in, data_covar_in, name_Y, covar) {
  data_covar <-
    data_covar_in[, colnames(data_covar_in) %in% covar, drop = FALSE]
  data_Y <-
    data_Y_in[rownames(data_Y_in) %in% rownames(data_covar), colnames(data_Y_in) ==
                name_Y, drop = FALSE]
  data_covar <-
    data_covar[rownames(data_covar) %in% rownames(data_Y), , drop = FALSE]
  data_covar <- data_covar[rownames(data_Y), , drop = FALSE]
  data_output <- data_Y
  data <- cbind(data_Y, data_covar)
  mod <- lm(data = data)
  data_output[, 1] <- as.data.frame(residuals(mod))
  return(data_output)
}

###ExWAS
ewas <-
  function(data_Xs_in = NULL,
           data_Y_in = NULL,
           name_Y,
           data_covar_in = NULL,
           covar = character(0),
           corr = "BH") {
    require(parallel)
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
                     rownames(data_Xs_in) %in% rownames(data_Y), , drop = FALSE]
      data_covar <- data_covar[rownames(data_Y), , drop = FALSE]
      data_Xs <- data_Xs[rownames(data_Y), , drop = FALSE]
      data_Y <- getresiduals_2df(data_Y, data_covar, name_Y, covar)
    } else{
      data_Y <-
        data_Y_in[rownames(data_Y_in) %in% rownames(data_Xs_in), colnames(data_Y_in) ==
                    name_Y, drop = FALSE]
      data_Xs <-
        data_Xs_in[rownames(data_Xs_in) %in% rownames(data_Y_in), , drop = FALSE]
      data_Xs <- data_Xs[rownames(data_Y), , drop = FALSE]
    }
    if (is.null(data_Y) == TRUE |
        is.null(data_Xs) == TRUE | !(name_Y %in% colnames(data_Y))) {
      stop("Inconsistent data")
    }
    
    ##computing p.values
    p.values <- lapply(1:ncol(data_Xs), function(x, data_Xs) {
      c(colnames(data_Xs)[x],
        confint(lm(Y ~ var1,
                   data = data.frame(
                     cbind(var1 = data_Xs[, x], Y = data_Y[, 1])
                   )))[2, ],
        summary(lm(Y ~ var1, data = data.frame(
          cbind(var1 = data_Xs[, x], Y = data_Y[, 1])
        )))$coefficients[2, ])
    },
    data_Xs)
    if (length(p.values) > 1) {
      p.values <-
        cbind(matrix(unlist(p.values), ncol = 7, byrow = TRUE)[, -6])
      
    } else{
      p.values <-
        as.data.frame(t(as.data.frame(unlist(p.values)))[, -6, drop = FALSE])
    }
    p.values <- as.data.frame(p.values)
    colnames(p.values) <-
      c("var", "conf - 2.5%", "conf - 97.5%", "Est", "Sd", "pVal")
    p.values <- p.values[p.values$var != "Intercept", ]
    p.values$pVal <- as.numeric(as.character(p.values$pVal))
    p.values.adj <- p.values
    pVal <- as.numeric(as.character(p.values$pVal))
    ##add correction for multiple testing
    if (corr == "None") {
      wh <- which(pVal <= 0.05)
      p.values.adj$pVal_adj <- pVal
    }
    if (corr == "Bon") {
      wh <- which(pVal <= 0.05 / nrow(p.values))
      p.values.adj$pVal_adj <- pVal * nrow(p.values)
    }
    if (corr == "BH") {
      wh <- which(p.adjust(pVal, "BH") <= 0.05)
      p.values.adj$pVal_adj <- p.adjust(pVal, "BH")
    }
    if (corr == "BY") {
      wh <- which(p.adjust(pVal, "BY") <= 0.05)
      p.values.adj$pVal_adj <- p.adjust(pVal, "BY")
    }
    if (!corr %in% c("Bon", "BH", "BY", "", "None"))
      stop("Please specify a known correction method for
                                                 multiple testing")
    wh <- p.values$var[wh]
    a <- list(wh, p.values.adj)
    names(a) <- c("selected", "pval")
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
      data_covar <- data_covar[rownames(data_Y),]
      data_Xs <- data_Xs[rownames(data_Y),]
      data_Y <- getresiduals_2df(data_Y, data_covar, name_Y, covar)
    } else{
      data_Y <-
        data_Y_in[rownames(data_Y_in) %in% rownames(data_Xs_in),
                  colnames(data_Y_in) ==
                    name_Y, drop = FALSE]
      data_Xs <-
        data_Xs_in[rownames(data_Xs_in) %in% rownames(data_Y_in), , drop = FALSE]
      data_Xs <- data_Xs[rownames(data_Y),]
    }
    data_Y <- data.matrix(data_Y)
    data_Xs <- data.matrix(data_Xs)
    model.enet <- cv.glmnet(data_Xs, data_Y, family = "gaussian",
                            alpha = 1)
    cvfit <- model.enet
    
    ##Calcul Y_predit
    Y_predit <- predict(cvfit, newx = data_Xs, s = "lambda.min")
    Y_predit <- Y_predit[rownames(Y_predit), ]
    
    ##liste des CPG selectionnés
    tmp_coeffs <- coef(cvfit, s = "lambda.min")
    cg_select <-
      data.frame(name = tmp_coeffs@Dimnames[[1]][tmp_coeffs@i + 1], coefficient = tmp_coeffs@x)
    
    cg_select <- cg_select$name[cg_select$name != "(Intercept)"]
    a <- list()
    if (length(cg_select) != 0) {
      a <- list("selected" = cg_select, "prediction" = Y_predit)
    } else{
      a <- list("selected" = character(), "prediction" = "no_prediction")
    }
    return(a)
    
  }


##DSA
DSAreg <-
  function(Exp,
           resp,
           family = gaussian,
           maxsize = 15,
           maxsumofpow = 2,
           maxorderint = 2) {
    Exp <- data.frame(cbind(data.frame(Exp), resp = resp))
    res <-
      DSA(
        resp ~ 1,
        data = Exp,
        family = family,
        maxsize = maxsize,
        maxsumofpow
        = maxsumofpow,
        maxorderint = maxorderint ,
        nsplits = 1,
        usersplits = NULL
      )
    form <- gsub("I[(]", "", colnames(coefficients(res)))
    form <- gsub("[*]", ":", gsub("[)]", "", gsub("[:^:]1", "", form)))
    if (length(grep(":", form)) > 0) {
      nam <- strsplit(form[grep("[:]", form)], ":")
      for (j in 1:length(nam)) {
        nam[[j]] <- gsub("[[:space:]]", "", nam[[j]])
        name <- nam[[j]][1]
        for (k in 2:length(nam[[j]]))
          name <- paste(name, ":", nam[[j]][k], sep = "")
        Exp <- cbind(Exp, name = apply(Exp[, nam[[j]]], 1, prod))
      }
    }
    form2 <- "resp~1"
    if (length(form) > 1)
      for(i in 2:length(form))
        form2 <- paste(form2, "+", form[i])
    res2 <- lm(form2, data = data.frame(Exp))
    ##decomment next line and change "prediction" to pred in the return line
    ##if outcome predicted by DSA is needed (not used presently)
    #pred <- predict(res2,Exp)
    coef <- summary(res2)$coefficients
    coef <- as.character(rownames(coef)[rownames(coef) != "Intercept"])
    
    return(list(selected = coef[coef != "(Intercept)"], pred = "prediction"))
  }


##################################################################
####3. defining some functions used to assess methods performance
#################################################################

sensitivity <- function(truepred, predfound) {
  return(length(truepred[truepred %in% predfound]) / length(truepred))
}
fdp <- function(truepred, predfound) {
  ##false discovery proportion
  if (length(predfound) == 0) {
    return(0)
  } else{
    return(length(predfound[!predfound %in% truepred]) / length(predfound))
  }
}
specificity <- function(truepred, predfound, n_base) {
  return((n_base - length(truepred) - length(predfound[!predfound %in% truepred])) /
           (n_base - length(truepred)))
}


#################################################################
##4. defining the simulation function which will be parallelized
#################################################################
##it first generates datasets, then applies methods and then assessed
##their performance
##simulation d'un jeu, application des méthodes, évaluation des méthodes
f0 <- function(x) {
  ##important: the parallelization is made on the seed
  set.seed(x)
  ##generating datasets
  simu <- simulator(
    E_true = dataExp_true,
    M_true = M1_true,
    Y_true = Y_true,
    n_mY = n_mY,
    R2_mY = R2_mY,
    BetamY = BetamY,
    n_yM = n_yM,
    Beta_yM = Beta_yM,
    n_yE = n_yE,
    Beta_yE = Beta_yE,
    test_and_training = TRUE
  )
  
  #simulated_data<-c(simulated_data,list(simu))
  
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
      as.data.frame(simu$M_train[, colnames(simu$M_train) %in% predBMI_M$ewas_BH$selected,
                                 drop =
                                   FALSE])
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
        cbind(predE_select_M$pval, rep(colnames(simu$E_train)[i], nrow(predE_select_M$pval)))
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
    ##step c
    n_exp_select <-
      length(unique(exp)) ##nb of exposures in reduced exposome
    
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
    names(predBMI_E_MITMdsa) <- c("selected", "pred")
    
  }
  ##storing results in a list
  predBMI_E <-
    c(
      predBMI_E,
      MITM = list(predBMI_E_MITM),
      MITMdsa = list(predBMI_E_MITMdsa),
      ReducedExp = list(predReducedExp)
    )
  print("oMITM")
  ###Control method : random sampling on a random set of exposures of same
  ##dimension as the reduced exposome of oMITM
  
  if (n_exp_select > 0) {
    tirage <-
      ewas(
        as.data.frame(simu$E_train)[, sample(colnames(as.data.frame(simu$E_train)), n_exp_select), drop =
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
      as.data.frame(simu$E_train[, colnames(simu$E_train) %in% predBMI_E$ewas_BH$selected, drop =
                                   FALSE])
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
        cbind(exp_affecting_M_all$pval, rep(
          colnames(simu$M_train)[i],
          nrow(exp_affecting_M_all$pval)
        ))
      list_temp_ewas_med <- c(list_temp_ewas_med, list(temp_ewas_med))
    }
    ewas_med <- do.call("rbind", list_temp_ewas_med)
    if (!is.null(ewas_med)) {
      ewas_med$pVal_adj_1 <- p.adjust(ewas_med$pVal, "BH")
      colnames(ewas_med)[8] <- "cpg"
      colnames(ewas_med)[1] <- "exp"
    }
    #step  b
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
        cbind(M_affecting_Y_all$pval, rep(colnames(select_E)[i], nrow(M_affecting_Y_all$pval)))
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
      ewas_med_tot$exp[ewas_med_tot$pVal_adj_2 <= 0.05 &
                         ewas_med_tot$pVal_adj_1 <= 0.05]
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
  print("lasso done")
  ##DSA
  predDSA <-
    DSAreg(
      Exp = simu$E_train,
      resp = as.data.frame(simu$Y_train),
      maxsize = floor(ncol(simu$E_train) / 10),
      maxsumofpow = 1,
      maxorderint = 1
    )
  predBMI_E <- c(predBMI_E, DSA = list(predDSA))
  print("DSA done")
  
  ##########################################
  ##assessing performance linked to reverse causality between M and Y
  ##(sensibility to predicted exposures) for all methods
  truepred <- simu$cpg_predicted$name
  for (k1 in 1:length(predBMI_M)) {
    truepred <- as.character(simu$cpg_predicted$name)
    predfound <- as.character(predBMI_M[[k1]]$selected)
    if (exists("predfound") & exists("truepred")) {
      if (length(predfound) == 0) {
        print("no predictors found")
      }
      a <- sensitivity(truepred, predfound)
      b <- specificity(truepred, predfound, ncol(simu$M_train))
      c <- fdp(truepred, predfound)
      d <- estimatedR2(simu$M_test, predfound, simu$Y_test)$r.squared
      # print(a)
      # print(b)
      # print(c)
      # print(d)
      remove(predfound, truepred)
    } else{
      if (!exists("predfound") & exists("truepred")) {
        a <- 0
        b <- specificity(truepred, numeric(0), ncol(simu$M_train))
        c <- 0
        d <- 0
        remove(truepred)
      } else{
        if (!exists("truepred") & exists("predfound")) {
          a <- 1
          b <- specificity(numeric(0), predfound, ncol(simu$M_train))
          c <- fdp(numeric(0), predfound)
          d <-
            estimatedR2(simu$M_test, predfound, simu$Y_test)$r.squared
          remove(predfound)
        } else{
          if (!exists("truepred") & !exists("predfound")) {
            a <- 1
            b <- 1
            c <- 0
            d <- 0
          }
        }
      }
    }
    predBMI_M[[k1]] <-
      c(
        predBMI_M[[k1]],
        sens_rev_caus = a,
        spec_rev_caus = b,
        fdp_rev_caus = c,
        R2_test_rev_caus = d
      )
    remove(a)
    remove(b)
    remove(c)
    remove(d)
    
  }
  
  ##########################################
  ##assessing performance for ExWAS on M
  truepred <- simu$cpg_predictors$name
  for (k1 in 1:length(predBMI_M)) {
    truepred <- as.character(simu$cpg_predictors$name)
    predfound <- as.character(predBMI_M[[k1]]$selected)
    if (exists("predfound") & exists("truepred")) {
      if (length(predfound) == 0) {
        print("no predictors found")
      }
      a <- sensitivity(truepred, predfound)
      b <- specificity(truepred, predfound, ncol(simu$M_train))
      c <- fdp(truepred, predfound)
      d <- estimatedR2(simu$M_test, predfound, simu$Y_test)$r.squared
      # print(a)
      # print(b)
      # print(c)
      # print(d)
      remove(predfound, truepred)
    } else{
      if (!exists("predfound") & exists("truepred")) {
        a <- 0
        b <- specificity(truepred, numeric(0), ncol(simu$M_train))
        c <- 0
        d <- 0
        remove(truepred)
      } else{
        if (!exists("truepred") & exists("predfound")) {
          a <- 1
          b <- specificity(numeric(0), predfound, ncol(simu$M_train))
          c <- fdp(numeric(0), predfound)
          d <-
            estimatedR2(simu$M_test, predfound, simu$Y_test)$r.squared
          remove(predfound)
        } else{
          if (!exists("truepred") & !exists("predfound")) {
            a <- 1
            b <- 1
            c <- 0
            d <- 0
          }
        }
      }
    }
    predBMI_M[[k1]] <-
      c(
        predBMI_M[[k1]],
        sens = a,
        spec = b,
        fdp = c,
        R2_test = d
      )
    remove(a)
    remove(b)
    remove(c)
    remove(d)
    
  }
  
  ##########################################
  ##assessing performance linked to reverse causality between E and Y
  ##( sensibility to predicted exposures) for all methods
  truepred <- simu$exp_predicted$name
  for (k1 in 1:length(predBMI_E)) {
    truepred <- as.character(simu$exp_predicted$name)
    predfound <- as.character(predBMI_E[[k1]]$selected)
    if (exists("predfound") & exists("truepred")) {
      if (length(predfound) == 0) {
        print("no predictors found")
      }
      a <- sensitivity(truepred, predfound)
      b <- specificity(truepred, predfound, ncol(simu$E_train))
      c <- fdp(truepred, predfound)
      d <- estimatedR2(simu$E_test, predfound, simu$Y_test)$r.squared
      # print(a)
      # print(b)
      # print(c)
      # print(d)
      remove(predfound, truepred)
    } else{
      if (!exists("predfound") & exists("truepred")) {
        a <- 0
        b <- specificity(truepred, numeric(0), ncol(simu$E_train))
        c <- 0
        d <- 0
        remove(truepred)
      } else{
        if (!exists("truepred") & exists("predfound")) {
          a <- 1
          b <- specificity(numeric(0), predfound, ncol(simu$E_train))
          c <- fdp(numeric(0), predfound)
          d <-
            estimatedR2(simu$E_test, predfound, simu$Y_test)$r.squared
          remove(predfound)
        } else{
          if (!exists("truepred") & !exists("predfound")) {
            a <- 1
            b <- 1
            c <- 0
            d <- 0
          }
        }
      }
    }
    predBMI_E[[k1]] <-
      c(
        predBMI_E[[k1]],
        sens_rev_caus = a,
        spec_rev_caus = b,
        fdp_rev_caus = c,
        R2_test_rev_caus = d
      )
    remove(a)
    remove(b)
    remove
    remove(d)
  }
  
  ##########################################
  ##assessing performance between E and Y (false detection)
  truepred <- character()
  for (k1 in 1:length(predBMI_E)) {
    truepred <- character()
    predfound <- as.character(predBMI_E[[k1]]$selected)
    if (exists("predfound") & exists("truepred")) {
      if (length(predfound) == 0) {
        print("no predictors found")
      }
      a <- sensitivity(truepred, predfound)
      b <- specificity(truepred, predfound, ncol(simu$E_train))
      c <- fdp(truepred, predfound)
      d <- estimatedR2(simu$E_test, predfound, simu$Y_test)$r.squared
      # print(a)
      # print(b)
      # print(c)
      # print(d)
      remove(predfound, truepred)
    } else{
      if (!exists("predfound") & exists("truepred")) {
        a <- 0
        b <- specificity(truepred, numeric(0), ncol(simu$E_train))
        c <- 0
        d <- 0
        remove(truepred)
      } else{
        if (!exists("truepred") & exists("predfound")) {
          a <- 1
          b <- specificity(numeric(0), predfound, ncol(simu$E_train))
          c <- fdp(numeric(0), predfound)
          d <-
            estimatedR2(simu$E_test, predfound, simu$Y_test)$r.squared
          remove(predfound)
        } else{
          if (!exists("truepred") & !exists("predfound")) {
            a <- 1
            b <- 1
            c <- 0
            d <- 0
          }
        }
      }
    }
    predBMI_E[[k1]] <-
      c(
        predBMI_E[[k1]],
        sens = a,
        spec = b,
        fdp = c,
        R2_test = d
      )
    remove(a)
    remove(b)
    remove(c)
    remove(d)
    
  }
  print("performance characterized")
  
  ##building the list with datasets generated + results of methods +
  ##performance to return
  A <-
    list(
      simu = simu,
      predBMI_E = predBMI_E,
      predBMI_M = predBMI_M,
      nl_exp_select = n_exp_select
    )
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
Y_true <- readRDS("20190612_ZBMI_scaled_a_utiliser_pour_simu.rds")
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
c_n_yM <- c(10, 18, 25, 100)

c_BetayM <- c(0.0001, 0.001, 0.01, 0.1, 0.5, 2)
c_BetayE <- c(0.0001, 0.001, 0.01, 0.1, 0.5, 2)
c_n_yE <- c(1, 3, 10, 25)
n_mY = 0
R2_mY = 0
BetamY = 0

##initialization of table of results

comp_method <-
  data.frame(
    Methods = vector(),
    Association_tested = vector(),
    Nb_true_predictors_of_BMI_in_M = numeric(0),
    Nb_predicted_by_BMI_in_M = numeric(0),
    Nb_predicted_by_BMI_in_E = numeric(0),
    Total_variability_of_BMI_explained_by_M = numeric(0),
    Total_variability_of_BMI_explained_by_M_measured =
      numeric(0),
    Mean_variability_of_M_explained_by_Y = numeric(0),
    Mean_variability_of_E_explained_by_Y = numeric(0),
    Mean_SD_variability_of_M_explained_by_Y = numeric(0),
    Mean_SD_variability_of_E_explained_by_Y = numeric(0),
    
    
    Number_iterations = numeric(0),
    Mean_number_exp_selected_to_be_randomly_tested =
      numeric(0),
    Mean_number_predictors_found = numeric(0),
    Mean_sensitivity_rv = numeric(0),
    Mean_specificity_rv = numeric(0),
    Mean_fdp_rv = numeric(0),
    Mean_R2_test_rv = numeric(0),
    Mean_sensitivity_truepred = numeric(0),
    Mean_specificity_truepred = numeric(0),
    Mean_fdp_truepred = numeric(0),
    Mean_R2_test_truepred = numeric(0),
    SD_number_predictors_found = numeric(0),
    SD_sensitivity_rv = numeric(0),
    SD_specificity_rv = numeric(0),
    SD_fdp_rv = numeric(0),
    SD_R2_test_rv = numeric(0),
    SD_sensitivity_truepred = numeric(0),
    SD_specificity_truepred = numeric(0),
    SD_fdp_truepred = numeric(0),
    SD_R2_test_truepred = numeric(0),
    Which_iteration = numeric(0)
  )


n = 1


##looping on the different vectors of parameters
for (i2 in 1:length(c_BetayE)) {
  BetayE <- c_BetayE[i2]
  for (i3 in 1:length(c_n_yE)) {
    n_yE <- c_n_yE[i3]
    for (i4 in 1:length(c_BetayM)) {
      BetayM <- c_BetayM[i4]
      for (i1 in 1:length(c_n_yM)) {
        n_yM <- c_n_yM[i1]
        
        
        n_row = nrow(comp_method)
        simulated_data <- list()
        list_predBMI_E <- list()
        list_predBMI_M <- list()
        list_nl_exp_select <- list()
        start_time <- Sys.time()
        ##parallelization of f0
        cl <-
          makeCluster(getOption("cl.cores", round(detectCores())))
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
            "Y_true",
            "BetayE",
            "n_yE",
            "n_yM",
            "BetayM",
            "n_mY",
            "R2_mY",
            "BetamY"
          )
        )
        clusterEvalQ(cl, list(library("boot"), library("reshape"), 
                              library("glmnet"), library("DSA")))
        results_1_jeu <- clusterApply(cl, 1:n_iter, f0)
        stopCluster(cl)
        simulated_data <-
          lapply(results_1_jeu, function(x)
            x$simu)
        
        ##structure of results priorized by methods and not anymore
        ##priorized by datasets
        list_predBMI_E <-
          lapply(results_1_jeu, function(x)
            x$predBMI_E)
        list_predBMI_M <-
          lapply(results_1_jeu, function(x)
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
                length(X$cpg_predictors$betas))
            ))
          )
        param_simu[2, ] <-
          c("Nb_predicted_by_BMI_in_M", "Fixed", mean(unlist(
            lapply(simulated_data, function(X)
              length(X$cpg_predicted$betas))
          )))
        param_simu[3, ] <-
          c("Nb_predicted_by_BMI_in_E", "Fixed", mean(unlist(
            lapply(simulated_data, function(X)
              length(X$exp_predicted$betas))
          )))
        param_simu[4, ] <-
          c("Total_variability_of_BMI_explained_by_M",
            "Fixed",
            R2_mY)
        param_simu[5, ] <-
          c("Total_variability_of_BMI_explained_by_M",
            "Measured",
            mean(unlist(
              lapply(simulated_data, function(X)
                (X$R2_mY_measured))
            )))
        param_simu[6, ] <-
          c("Mean_variability_of_M_explained_by_Y",
            "Measured",
            mean(unlist(
              lapply(simulated_data, function(X)
                (X$R2_yM_mean))
            )))
        param_simu[7, ] <-
          c("Mean_variability_of_E_explained_by_Y",
            "Measured",
            mean(unlist(
              lapply(simulated_data, function(X)
                (X$R2_yE_mean))
            )))
        param_simu[8, ] <-
          c("Mean_SD_variability_of_M_explained_by_Y",
            "Measured",
            mean(unlist(
              lapply(simulated_data, function(X)
                (X$R2_yM_SD))
            )))
        param_simu[9, ] <-
          c("Mean_SD_variability_of_E_explained_by_Y",
            "Measured",
            mean(unlist(
              lapply(simulated_data, function(X)
                (X$R2_yE_SD))
            )))
        param_simu[10, ] <- c("Number_iterations", "Fixed", n_iter)
        param_simu[11, ] <-
          c(
            "Mean_number_exp_selected_to_be_randomly_tested",
            "Measured",
            mean(unlist(list_nl_exp_select))
          )
        
        ##summarizing each method performance by a line in comp_method datadrame 
        ##within this scenario
        for (k1 in (1:length(list_predBMI_M[[1]]))) {
          comp_method[n_row + k1, ] <-
            c(
              names(list_predBMI_M[[1]][k1]),
              "BMI - M",
              param_simu[1, 3],
              param_simu[2, 3],
              param_simu[3, 3],
              param_simu[4, 3],
              param_simu[5, 3],
              param_simu[6, 3],
              param_simu[7, 3],
              param_simu[8, 3],
              param_simu[9, 3],
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
              
              mean(unlist(
                lapply(list_predBMI_M, function(X)
                  (X[[k1]][[7]]))
              ), na.rm = TRUE),
              mean(unlist(
                lapply(list_predBMI_M, function(X)
                  (X[[k1]][[8]]))
              ), na.rm = TRUE),
              mean(unlist(
                lapply(list_predBMI_M, function(X)
                  (X[[k1]][[9]]))
              ), na.rm = TRUE),
              mean(unlist(
                lapply(list_predBMI_M, function(X)
                  (X[[k1]][[10]]))
              ), na.rm = TRUE),
              
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
              
              sd(unlist(
                lapply(list_predBMI_M, function(X)
                  (X[[k1]][[7]]))
              ), na.rm = TRUE),
              sd(unlist(
                lapply(list_predBMI_M, function(X)
                  (X[[k1]][[8]]))
              ), na.rm = TRUE),
              sd(unlist(
                lapply(list_predBMI_M, function(X)
                  (X[[k1]][[9]]))
              ), na.rm = TRUE),
              sd(unlist(
                lapply(list_predBMI_M, function(X)
                  (X[[k1]][[10]]))
              ), na.rm = TRUE),
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
              param_simu[3, 3],
              param_simu[4, 3],
              param_simu[5, 3],
              param_simu[6, 3],
              param_simu[7, 3],
              param_simu[8, 3],
              param_simu[9, 3],
              n_iter,
              param_simu[11, 3],
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
              
              NA,
              NA,
              NA,
              NA,
              
              
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
              NA,
              NA,
              NA,
              NA,
              n
            )
        }
        print(n)
        
        ##storing generated datasets and methods results for this scenario
        list_list_predBMI_E <-
          c(list_list_predBMI_E, list(list_predBMI_E))
        list_list_predBMI_M <-
          c(list_list_predBMI_M, list(list_predBMI_M))
        list_list_nl_exp_select <-
          c(list_list_nl_exp_select, list(list_nl_exp_select))
        end_time <- Sys.time()
        end_time - start_time
        ##saving
        saveRDS(comp_method, "comp_method_rev_caus_sans_mY.Rds")
        saveRDS(
          simulated_data,
          file = paste(
            n,
            '_simulated_data_scenario_iteration_rev_caus_sans_mY.Rds'
          )
        )
        saveRDS(list_list_predBMI_E,
                "list_list_predBMI_E_rev_caus_sans_mY.rds")
        saveRDS(list_list_predBMI_M,
                "list_list_predBMI_M_rev_caus_sans_mY.rds")
        saveRDS(list_list_nl_exp_select,
                "list_list_nl_exp_select_rev_caus_sans_mY.rds")
        
        remove(simulated_data)
        remove(list_predBMI_E)
        remove(list_predBMI_M)
        remove(list_nl_exp_select)
        
        n = n + 1
        
        print(n)
      }
    }
  }
}


##save
saveRDS(comp_method, "comp_method_rev_caus_sans_mY.Rds")
#saveRDS(list_simulated_data,"list_simulated_data_mediation.Rds")
saveRDS(list_list_predBMI_E,
        "list_list_predBMI_E_rev_caus_sans_mY.Rds")
saveRDS(list_list_predBMI_M,
        "list_list_predBMI_M_rev_caus_sans_mY.Rds")
saveRDS(list_list_nl_exp_select,
        "list_list_nl_exp_select_rev_caus_sans_mY.Rds")
