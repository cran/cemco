#' This function initialize the parameters for the EM solution
#' of the CemCO algorithm
#'
#' @param data numeric matrix of data to use to cluster
#'
#' @param G the desired number of clusters
#'
#' @param Y numeric matrix of data to use as covariates
#'
#' @param factor numeric value to initialize the covariates effect

InitEM <- function(data, G, Y, fator=1){

  P <- ncol(data)
  p <- ncol(Y)

  f <- list()
  for(i in 1:G){
    f[[i]] <- rep(list(rep(0, ncol(data))),p)
  }

  X_lag <- matrix(nrow = nrow(data), ncol=ncol(data))
  for(i in 1:P){
    df_model <- data.frame(cbind(data[,i],Y))
    names(df_model) <- paste0("X", paste0(1:(ncol(Y)+1)))
    lin_mod <- lm("X1~.", data=df_model)
    for(j in 1:p){
      for(k in 1:G){
        f[[k]][[j]][i] <- lin_mod$coefficients[j+1]*fator
      }
    }
    df_model <- data.frame(cbind(data[,i],Y*fator))
    names(df_model) <- paste0("X", paste0(1:(ncol(Y)+1)))
    X_lag[,i] <- data[,i] - predict(lin_mod, newdata=df_model)*fator - fator*lin_mod$coefficients[1]
  }

  fit <- Mclust(X_lag, G=G)

  alpha <- fit$parameters$pro

  mus <- lapply(apply(fit$parameters$mean, 2, list), unlist)

  variances <- list()
  vars <- fit$parameters$variance$sigma
  ini <- 1
  end <- P*P
  for(i in 1:G){
    variances[[i]] <- matrix(vars[ini:end], nrow=P)
    ini <- end+1
    end <- end+P*P
  }

  inicialization <- list(mu = mus, sig = variances, alpha = alpha, beta = f)
  return(inicialization)
}


#' E-Step : Calculating probabilities of each cluster
#' This function calculates the E step of the EM implementation
#' of the CemCO algorithm
#'
#' @param data numeric matrix of data to use to cluster
#'
#' @param Y numeric matrix of data to use as covariates
#'
#' @param phi list of all parameters fitted on the EM algorithm
#'
#' @param G the desired number of clusters
#' 
#' 
#' @export
EStep <- function(data, Y, phi, G) {
  prob <- matrix(nrow=nrow(data), ncol=G)
  for(grupo in 1:G){
    for(i in 1:nrow(data)){
      media <- phi$mu[[grupo]]
      for(p in 1:ncol(Y)){
        media <- media + phi$beta[[grupo]][[p]]*Y[i,p]
      }
      prob[i,grupo] <- mvtnorm::dmvnorm(matrix(data[i,], ncol=length(media)),
                               media,
                               round(phi$sig[[grupo]],10))+10^(-10)
    }
  }
  return(prob/rowSums(prob))
}

#' This function calculates the sum of all covariates effects
#' @param data numeric matrix of data to use to cluster
#'
#' @param Y numeric matrix of data to use as covariates
#'
#' @param phi list of all parameters fitted on the EM algorithm
#'
#' @param G the desired number of clusters

XBetaCalculus <- function(data, Y, phi, G){
  x_beta <- list()
  for(i in 1:G){
    x_beta[[i]] <- matrix(0,nrow(data), ncol(data))
  }
  for(i in 1:G){
    for(j in 1:ncol(Y)){
      x_beta[[i]] <- x_beta[[i]] + sweep(
        do.call(cbind, rep(list(Y[,j]),ncol(data))),
        MARGIN=2,
        phi$beta[[i]][[j]],
        `*`)
    }
  }
  return(x_beta)
}

#' This function calculates the M step of the EM implementation of the CemCO algorithm
#' @param data numeric matrix of data to use to cluster
#'
#' @param Y numeric matrix of data to use as covariates
#'
#' @param phi list of all parameters fitted on the EM algorithm
#'
#' @param probs n x G numeric matrix with the probability of
#'              each observation belong to each cluster
#'
#' @param G the desired number of clusters

MStep <- function(data, Y, phi, probs, G) {

  x_beta <- XBetaCalculus(data, Y, phi, G)

  for(i in 1:G){
    for(j in 1:ncol(Y)){
      phi$beta[[i]][[j]] <- colSums(
        (sweep(data, 2, phi$mu[[i]])-x_beta[[i]]+sweep(
          do.call(cbind, rep(list(Y[,j]),ncol(data))),
          MARGIN=2,
          phi$beta[[i]][[j]],
          `*`))*Y[,j]*as.vector(matrix(
            rep(probs[,i],ncol(data)), ncol=ncol(data), byrow = FALSE)))/
        sum(Y[,j]^(2)*probs[,i])
      x_beta <- XBetaCalculus(data, Y, phi, G)
    }
    phi$mu[[i]] <- colSums((data-x_beta[[i]])*probs[,i])/colSums(probs)[i]
  }

  covs <- lapply(1:G, function(i) cov.wt((data-x_beta[[i]]), probs[,i]))
  phi$sig <- lapply(covs, "[[", "cov")

  phi$alpha <- colMeans(probs)

  return(phi)
}

#' This function calculates the Log Likelihood of the CemCO algorithm
#' @param data numeric matrix of data to use to cluster
#'
#' @param Y numeric matrix of data to use as covariates
#'
#' @param phi list of all parameters fitted on the EM algorithm
#'
#' @param probs n x G numeric matrix with the probability of
#'              each observation belong to each cluster
#'
#' @param G the desired number of clusters
#' @export
LogLike <- function(data, Y, phi, G) {
  prob <- matrix(nrow=nrow(data), ncol=G)
  for(grupo in 1:G){
    for(i in 1:nrow(data)){
      media <- phi$mu[[grupo]]
      for(p in 1:ncol(Y)){
        media <- media + phi$beta[[grupo]][[p]]*Y[i,p]
      }
      prob[i,grupo] <- phi$alpha[grupo]*mvtnorm::dmvnorm(matrix(data[i,], ncol=length(media)),
                                                media,
                                                round(phi$sig[[grupo]],10))
    }
  }
  sum(log(rowSums(prob)))
}



#' This function calculates the CemCO algorithm using multiple threads
#' @param data numeric matrix of data to use to cluster
#'
#' @param y numeric matrix of data to use as covariates
#'
#' @param G the desired number of clusters
#'
#' @param max_iter numeric value with maximum number of iterations to consider in log likelihood convergence.
#'
#' @param n_start numeric value representing how many times the EM algorithm should be initialized with different start values.
#'
#' @param cores numeric value how many cores the fit should use
#' @export
CemCO <- function(data, y, G, max_iter=100, n_start=20, cores=4) {
  fatores <- seq(0,2,by=2/n_start)
  registerDoParallel(cores)
  phis <- foreach(start=1:n_start, .packages='mvtnorm') %dopar% {
    phi <- InitEM(data,G = G, Y=y, fator=fatores[start])
    likelihoods <- list()
    to_compare <- LogLike(data, y, phi, G)
    for(i in 1:max_iter) {
      oldphi <- phi
      probs <- EStep(data, y, phi, G)
      phi <- MStep(data, y, phi, probs, G)
      likelihoods[[i]] <- LogLike(data, y, phi, G)
      if((likelihoods[[i]] - to_compare) < 0.01)
        break
      to_compare <- likelihoods[[i]]
    }
    list(phi, likelihoods[[i]])
  }
  stopImplicitCluster()
  id_best <- which.max(unlist(lapply(phis, function(x) x[[2]])))
  return(phis[[id_best]])
}


#' CemCOVar - This function initialize the parameters for the EM solution of the CemCO algorithm
#'
#' @param data numeric matrix of data to use to cluster
#'
#' @param G the desired number of clusters
#'
#' @param Y numeric matrix of data to use as covariates
#'
#' @param factor numeric value to initialize the covariates effect

InitEMVar <- function(data, G, Y, fator){

  P <- ncol(data)
  p <- ncol(Y)

  f <- list()
  for(i in 1:G){
    f[[i]] <- rep(list(rep(0, ncol(data))),p)
  }

  X_lag <- matrix(nrow = nrow(data), ncol=ncol(data))
  for(i in 1:P){
    df_model <- data.frame(cbind(data[,i],Y))
    names(df_model) <- paste0("X", paste0(1:(ncol(Y)+1)))
    lin_mod <- lm("X1~.", data=df_model)
    for(j in 1:p){
      for(k in 1:G){
        f[[k]][[j]][i] <- lin_mod$coefficients[j+1]*fator
      }
    }
    df_model <- data.frame(cbind(data[,i],Y*fator))
    names(df_model) <- paste0("X", paste0(1:(ncol(Y)+1)))
    X_lag[,i] <- data[,i] - predict(lin_mod, newdata=df_model)*fator - fator*lin_mod$coefficients[1]
  }

  fit <- Mclust(X_lag, G=G, control = emControl(itmax=10000))

  alpha <- fit$parameters$pro

  mus <- lapply(apply(fit$parameters$mean, 2, list), unlist)

  par <- matrix(0.0001, ncol=P, nrow=P)
  par[,1] <- 1
  pars_list <- list()
  for(i in 1:G) pars_list[[i]] <- par


  E = matrix(c(1,1,1,1), nrow=2)

  variances <- list()
  vars <- fit$parameters$variance$sigma
  ini <- 1
  end <- P*P
  for(i in 1:G){
    variances[[i]] <- matrix(vars[ini:end], nrow=P)
    ini <- end+1
    end <- end+P*P
  }

  inicialization <- list(mu = mus, sig = variances, par = pars_list, alpha = alpha, beta = f)
  return(inicialization)
}


#' This function creates a P x P matrix representing the the variance effect for an observation
#'
#' @param P numeric value with the dimension of the data
#'
#' @param yi numeric value, covariate value for the observation
#'
#' @param matrix with the predicted value of the variance effect


CreateL <- function(P, yi, parameters){
  E <- matrix(0,nrow=P, ncol=P)
  for(i in 1:P){
    E[i,i] <- parameters[i,1]
    for(j in 2:ncol(parameters)){
      E[i,i] <- E[i,i]+parameters[i,j]*yi
    }
  }
  return(E)
}


#' E-Step : Calculating probabilities of each cluster
#' This function calculates the E step of the EM implementation of the CemCO algorithm
#' @param data numeric matrix of data to use to cluster
#'
#' @param Y numeric matrix of data to use as covariates
#'
#' @param phi list of all parameters fitted on the EM algorithm
#'
#' @param G the desired number of clusters
#'
#' @param y_cov numeric array of data to use as a covariate just for the variance effect
#' @export
EStepVar <- function(data, Y, phi, G, y_cov) {
  prob <- matrix(nrow=nrow(data), ncol=G)
  for(grupo in 1:G){
    for(i in 1:nrow(data)){
      media <- phi$mu[[grupo]]
      for(p in 1:ncol(Y)){
        media <- media + phi$beta[[grupo]][[p]]*Y[i,p]
      }
      L <- CreateL(ncol(data), y_cov[i], phi$par[[grupo]])
      vars <- L%*%phi$sig[[grupo]]%*%L
      prob[i,grupo] <- mvtnorm::dmvnorm(matrix(data[i,], ncol=length(media)),
                               media,
                               round(vars,10))+10^(-10)
    }
  }
  return(prob/rowSums(prob))
}

#' This function defines the objective function to be optimized in the M Step for the algorithm with variance effect.
#'
#' @param x parameter to be optimized
#'
#' @param parameters list with several internal values required to define the objective function

OptimFunctionVar <- function(x, parameters){

  E = parameters[[1]]
  y = parameters[[2]]
  X = parameters[[3]]
  beta = parameters[[5]]
  prob = parameters[[6]]
  Y = parameters[[7]]

  P <- ncol(X)
  f <- matrix(0,nrow=P, ncol=P)

  for(i in 1:length(y)){
    media = parameters[[4]]
    if(ncol(Y) > 1){
      for(p in 1:ncol(Y)){
        media <- media + beta[[p]]*Y[i,p]
      }
    } else {
      media <- media + unlist(beta)*Y[i,1]
    }


    D = matrix(unlist(X[i,] - media), ncol=1)
    C = D%*%t(D)
    L = CreateL(ncol(X), y[i], matrix(x, ncol=P, byrow = FALSE))

    E_inv <- ginv(E)
    L_inv <- ginv(L)
    f = f+prob[i]*(L-0.5*E_inv%*%L_inv%*%C-0.5*t(E_inv)%*%L_inv%*%C)
  }
  return(sum(f)^2)
}


#' This function estimate the variance effect by optimizing the log likelihood on the M Step
#'
#' @param data numeric matrix of data to use to cluster
#'
#' @param y_cov numeric array of data to use as a covariate just for the variance effect
#'
#' @param phi list of all parameters fitted on the EM algorithm
#'
#' @param probs n x G numeric matrix with the probability of each observation belong to each cluster
#'
#' @param G the desired number of clusters
#'
#' @param Y numeric matrix of data to use as covariates

VarEstimateVar <- function(data, y_cov, phi, probs, G, Y){
  vars <- list()
  for(l in 1:G){
    vars[[l]] <- matrix(optim(fn=OptimFunctionVar, par=phi$par[[l]], parameters = list(phi$sig[[l]], y_cov,
                              data, phi$mu[[l]],phi$beta[[l]], probs[,l], Y))$par,
                        ncol=ncol(data), byrow = FALSE)
  }
  return(vars)
}

#' This function calculates the sum of all covariates effects
#' @param data numeric matrix of data to use to cluster
#'
#' @param Y numeric matrix of data to use as covariates
#'
#' @param phi list of all parameters fitted on the EM algorithm
#'
#' @param G the desired number of clusters

XBetaCalculusVar <- function(data, Y, phi, G){
  x_beta <- list()
  for(i in 1:G){
    x_beta[[i]] <- matrix(0,nrow(data), ncol(data))
  }
  for(i in 1:G){
    for(j in 1:ncol(Y)){
      x_beta[[i]] <- x_beta[[i]] + sweep(
        do.call(cbind, rep(list(Y[,j]),ncol(data))),
        MARGIN=2,
        phi$beta[[i]][[j]],
        `*`)
    }
  }
  return(x_beta)
}

#' This function calculates the M step of the EM implementation of the CemCOVar algorithm
#'
#' @param data numeric matrix of data to use to cluster
#'
#' @param Y numeric matrix of data to use as covariates
#'
#' @param phi list of all parameters fitted on the EM algorithm
#'
#' @param probs n x G numeric matrix with the probability of each observation belong to each cluster
#'
#' @param G the desired number of clusters
#'
#' @param y_cov numeric array of data to use as a covariate just for the variance effect

MStepVar <- function(data, Y, phi, probs, G, y_cov) {

  x_beta <- XBetaCalculusVar(data, Y, phi, G)

  for(i in 1:G){
    for(j in 1:ncol(Y)){
      phi$beta[[i]][[j]] <- colSums(
        (sweep(data, 2, phi$mu[[i]])-x_beta[[i]]+sweep(
          do.call(cbind, rep(list(Y[,j]),ncol(data))),
          MARGIN=2,
          phi$beta[[i]][[j]],
          `*`))*Y[,j]*as.vector(matrix(
            rep(probs[,i],ncol(data)), ncol=ncol(data), byrow = FALSE)))/
        sum(Y[,j]^(2)*probs[,i])
      x_beta <- XBetaCalculusVar(data, Y, phi, G)
    }
    phi$mu[[i]] <- colSums((data-x_beta[[i]])*probs[,i])/colSums(probs)[i]
  }

  Ls <- list()
  Ds <- list()
  P <- ncol(data)

  mat <- list()
  for(i in 1:G){
    mat[[i]] <- matrix(0,nrow=ncol(data),ncol=ncol(data))
  }

  for(i in 1:nrow(data)){
    for(k in 1:G){
      Ls[[k]] <- CreateL(P, y_cov[i], phi$par[[k]])

      media <- phi$mu[[k]]
      for(j in 1:ncol(Y)){
        media <- media + phi$beta[[k]][[j]]*Y[i,j]
      }
      Ds[[k]] <- matrix(unlist(data[i,]-media),
                        ncol=1)

      m <- solve(Ls[[k]])%*%Ds[[k]]
      mat[[k]] <- mat[[k]] + probs[i,k]*((m)%*%t(m))
    }
  }
  for(j in 1:G){
    mat[[j]] <- mat[[j]]/sum(probs[,j])
  }

  phi$sig <- mat

  phi$par <- VarEstimateVar(data, y_cov, phi, probs, G, Y)

  phi$alpha <- colMeans(probs)

  return(phi)
}

#' This function calculates the Log Likelihood of the CemCOVar algorithm
#'
#' @param data numeric matrix of data to use to cluster
#'
#' @param Y numeric matrix of data to use as covariates
#'
#' @param phi list of all parameters fitted on the EM algorithm
#'
#' @param G the desired number of clusters
#'
#' @param y_cov numeric array of data to use as a covariate just for the variance effect
#' @export
LogLikeVar <- function(data, Y, phi, G, y_cov) {
  prob <- matrix(nrow=nrow(data), ncol=G)
  for(grupo in 1:G){
    for(i in 1:nrow(data)){
      media <- phi$mu[[grupo]]
      for(p in 1:ncol(Y)){
        media <- media + phi$beta[[grupo]][[p]]*Y[i,p]
      }
      L <- CreateL(ncol(data), y_cov[i], phi$par[[grupo]])
      var <- L%*%phi$sig[[grupo]]%*%L
      prob[i,grupo] <- phi$alpha[grupo]*mvtnorm::dmvnorm(matrix(data[i,], ncol=length(media)),
                                                media,
                                                round(var,10))
    }
  }
  sum(log(rowSums(prob)))
}

#' This function calculates the CemCOVar algorithm using multiple threads
#'
#' @param data numeric matrix of data to use to cluster
#'
#' @param y numeric matrix of data to use as covariates
#'
#' @param G the desired number of clusters
#'
#' @param y_cov numeric array of data to use as a covariate just for the variance effect
#'
#' @param max_iter numeric value with maximum number of iterations to consider in log likelihood convergence.
#'
#' @param n_start numeric value representing how many times the EM algorithm should be initialized with different start values.
#'
#' @param cores numeric value how many cores the fit should use
#' @export
CemCOVar <- function(data, y, G, y_cov, max_iter=100, n_start=20, cores=4) {
  fatores <- seq(0,2,by=2/n_start)
  registerDoParallel(cores)
  phis <- foreach(start=1:n_start, .packages='mvtnorm') %dopar% {
    phi <- InitEMVar(data,G = G, Y=y, fator=fatores[start])
    likelihoods <- list()
    to_compare <- LogLikeVar(data, y, phi, G, y_cov)
    for(i in 1:max_iter) {
      oldphi <- phi
      probs <- EStepVar(data, y, phi, G, y_cov)
      phi <- MStepVar(data, y, phi, probs, G, y_cov)
      likelihoods[[i]] <- LogLikeVar(data, y, phi, G, y_cov)
      if((likelihoods[[i]] - to_compare) < 0.01)
        break
      to_compare <- likelihoods[[i]]
    }
    if(i>1){
      list(oldphi, likelihoods[[i-1]])
    } else {
      list(oldphi, likelihoods[[i]])
    }
  }
  stopImplicitCluster()
  id_best <- which.max(unlist(lapply(phis, function(x) x[[2]])))[1]
  return(phis[[id_best]])
}


