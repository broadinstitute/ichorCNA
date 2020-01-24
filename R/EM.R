# file:   EM.R
# author: Gavin Ha, Ph.D.
#         Fred Hutchinson Cancer Research Center
# contact: <gha@fredhutch.org>
# website: https://GavinHaLab.org
#
# ichorCNA website: https://github.com/GavinHaLab/ichorCNA
# date:   January 6, 2020
# description: Hidden Markov model (HMM) to analyze Ultra-low pass whole genome sequencing (ULP-WGS) data.
# This script is the main script to run the HMM.

runEM <- function(copy, chr, chrInd, param, maxiter, verbose = TRUE, 
									estimateNormal = TRUE, estimatePloidy = TRUE, 
									estimateVar = TRUE, estimatePrecision = TRUE,
                  estimateTransition = TRUE, estimateInitDist = TRUE, estimateSubclone = TRUE,
									likChangeConvergence = 1e-3) {
  
  if (nrow(copy) != length(chr) || nrow(copy) != length(chrInd)) {
    stop("runEM: Length of inputs do not match for one of: copy, chr, chrInd")
  }
  
  if (is.null(param$ct) || is.null(param$lambda) || is.null(param$nu) ||
      is.null(param$kappa)) {
    stop("runEM: Parameter missing, ensure all parameters exist as columns in",
         "data frame: ct, lambda, nu, kappa")
  }
  
  S <- param$numberSamples
  K <- length(param$ct)
  Z <- sum(param$ct.sc) #number of subclonal states (# repeated copy number states)
  KS <- nrow(param$jointStates)#K ^ S
  N <- nrow(copy)
  rho <- matrix(0, KS, N)
  lambdas <- array(0, dim = c(K, S, maxiter))  # State Variances
  #covars <- rep(list(NULL), maxiter) #State covariance matrices
  vars <- array(0, dim = c(K, S, maxiter))
  phi <- matrix(NA, length(param$phi_0), maxiter)					 # Tumour ploidy
  n <- matrix(NA, S, maxiter)						 # Normal contamination
  sp <- matrix(NA, S, maxiter)     # cellular prevalence (tumor does not contain event)
  piG <- matrix(0, KS, maxiter)     # Initial state distribution
  A <- param$A
  converged <- FALSE               # Flag for convergence
  Zcounts <- matrix(0, KS, KS)
  loglik <- rep(0, maxiter)
  mus <- array(0, dim = c(KS, S, maxiter))     # State means
  py <- matrix(0, KS, N)            # Local evidence
  
  ptmTotal <- proc.time() # start total timer
  
  # SET UP
  # Set up the chromosome indices and make cell array of chromosome indicies
  chrs <- levels(chr)              # WARNING, gets resorted here...
  chrsI <- vector('list', length(chrs))
  # initialise the chromosome index and the init state distributions
  for(i in 1:length(chrs)) {
    chrsI[i] <- list(which(chr == chrs[i]))
  }
  
  # INITIALIZATION
  if (verbose) { message("runEM: Initialization") }
  i <- 1
  piG[, i] <- normalize(param$kappa)
  n[, i] <- param$n_0
  sp[, i] <- param$sp_0
  phi[, i] <- param$phi_0
  lambdas[, , i] <- param$lambda #matrix(param$lambda, nrow = K, ncol = S, byrow = TRUE)
  lambdasKS <- as.matrix(expand.grid(as.data.frame(lambdas[, , i]))) #as.matrix(expand.grid(rep(list(param$lambda), S)))
  #covars[[i]] <- param$covar
  if (param$likModel == "Gaussian"){
    vars[, , i] <- param$var
    varsKS <- as.matrix(expand.grid(as.data.frame(vars[, , i])))
  }
  mus[, , i] <- as.matrix(get2and3ComponentMixture(param$jointCNstates, param$jointSCstatus, n[, i], sp[, i], phi[, i]))
  
  # Likelihood #
  if (param$likModel == "t"){
    for (ks in 1:KS) {
      probs <- tdistPDF(copy, mus[ks, , i], lambdasKS[ks, ], param$nu)
      py[ks, ] <- apply(probs, 1, prod) # multiply across samples for each data point to get joint likelihood.
    }
  }else if (param$likModel == "Gaussian"){
    message("Using ", param$likModel, " emission model.")
    #py <- t(sapply(1:KS, function(ks) {
    #  mvGaussDistPDF(copy, mus[ks, , i], diag(varsKS[ks, ])) # not true covariance matrix
    #}))
    # py <- t(sapply(1:KS, function(ks) {
    #   probs <- normalpdf(copy, mus[ks, , i], varsKS[ks, ])
    #   apply(probs, 1, prod)
    # }))
    py <- getNormLik(copy, mus[, , i, drop = FALSE], varsKS, param$sw)
  }
  
  loglik[i] <- -Inf
  
  while(!converged && (i < maxiter)) {
    ptm <- proc.time()
    #if (verbose) { message(paste("runEM: EM iteration:", i, "Log likelihood:", loglik[i])) }
    i <- i + 1
    
    ################ E-step ####################
    if (verbose) { message("runEM iter", i-1 , ": Expectation") }
    Zcounts <- matrix(0, KS, KS)
    for (j in 1:length(chrsI)) {
      I <- intersect(chrsI[[j]], which(chrInd))
      if (length(I) > 0){
        output <- .Call("forward_backward", piG[, i - 1], param$A, py[, I],PACKAGE = "HMMcopy")
        rho[, I] <- output$rho
        loglik[i] <- loglik[i] + output$loglik
        Zcounts <- Zcounts + t(colSums(aperm(output$xi, c(3, 2, 1))))
      }
    }
    
    ## marginalize rho by sample ##
    #rhoS <- list()
    #for (s in 1:S){
    #  rhoS[[s]] <- matrix(0, nrow = K, ncol = N)
    #  for (k in 1:K){
    #    rhoS[[s]][k, ] <- colSums(rho[param$jointCNstates[, s] == k, chrInd])
    #  }
    #}
    
    ################ M-step ####################
    if (verbose) { message("runEM iter", i-1 , ": Maximization") }
    ptm.em <- proc.time()
    if (param$likModel == "Gaussian"){
      output <- estimateGaussianParamsMap(copy[chrInd, ], n[, i - 1], sp[, i - 1], phi[, i - 1], 
                                          varsKS, piG[, i - 1], A, param,
                                          rho[, chrInd], Zcounts, estimateNormal = estimateNormal, 
                                          estimatePloidy = estimatePloidy,
                                          estimateVar = estimateVar, 
                                          estimateTransition = estimateTransition,
                                          estimateInitDist = estimateInitDist, 
                                          estimateSubclone = estimateSubclone)
      #vars[, , i] <- output$var
      varsKS <- output$var
      if (verbose == TRUE) {
        for (s in 1:S){
          message("Sample", s, " n=", signif(output$n[s], digits = 4), ", sp=", signif(output$sp[s], digits = 4), 
                  ", phi=", signif(output$phi[s], digits = 4))
                  #", var=", paste0(signif(output$var[, s], digits=5), collapse = ","))
        }     
      }
    }else{ # Student's-t model
      ptm.em <- proc.time()
      output <- estimateTParamsMap(copy[chrInd, ], n[, i - 1], sp[, i - 1], phi[, i - 1], 
                                lambdas[, , i - 1], piG[, i - 1], A, param,
                                rho[, chrInd], Zcounts, 
                                estimateNormal = estimateNormal, estimatePloidy = estimatePloidy,
                                estimatePrecision = estimatePrecision, estimateTransition = estimateTransition,
                                estimateInitDist = estimateInitDist, estimateSubclone = estimateSubclone)
      
      lambdas[, , i] <- output$lambda
      if (verbose == TRUE) {
        
        for (s in 1:S){
          message("Sample", s, " n=", signif(output$n[s], digits = 4), ", sp=", signif(output$sp[s], digits = 4), 
                  ", phi=", signif(output$phi[s], digits = 4), 
                  ", lambda=", paste0(signif(output$lambda[, s], digits=5), collapse = ","),
                  ", F=", signif(output$F, digits=4))
        }     
      }
    }
    elapsedTime <- proc.time() - ptm.em
    message("runEM iter", i-1 , ": ", format(elapsedTime[3] / 60, digits = 2), "min.")
    n[, i] <- output$n
    sp[, i] <- output$sp
    phi[, i] <- output$phi
    piG[, i] <- output$piG
    A <- output$A
    estF <- output$F
    
    # Recalculate the likelihood
    varsKS <- as.matrix(expand.grid(as.data.frame(vars[, , i])))
    mus[, , i] <- as.matrix(get2and3ComponentMixture(param$jointCNstates, param$jointSCstatus, n[, i], sp[, i], phi[, i]))
    if (param$likModel == "t"){
      for (ks in 1:KS) {
        probs <- tdistPDF(copy, mus[ks, , i], lambdasKS[ks, ], param$nu)
        py[ks, ] <- apply(probs, 1, prod) # multiply across samples for each data point to get joint likelihood.
      }
    }else if (param$likModel == "Gaussian"){
      #py <- t(sapply(1:KS, function(ks) {
      #  mvGaussDistPDF(copy, mus[ks, , i], diag(varsKS[ks, ])) # not true covariance matrix
      #}))
      # py <- t(sapply(1:KS, function(ks) {
      #   probs <- normalpdf(copy, mus[ks, , i], varsKS[ks, ])
      #   apply(probs, 1, prod)
      # }))
      py <- getNormLik(copy, mus[, , i, drop = FALSE], varsKS, param$sw)
    }
    
    prior <- priorProbs(n[, i], sp[, i], phi[, i], lambdas[, , i], varsKS, piG[, i], A, param, 
                        estimateNormal = estimateNormal, estimatePloidy = estimatePloidy,
                        estimateVar = estimateVar, estimatePrecision = estimatePrecision, 
                        estimateTransition = estimateTransition, estimateInitDist = estimateInitDist, 
                        estimateSubclone = estimateSubclone)
    
    
    # check converence 
    loglik[i] <- loglik[i] + prior$prior
    elapsedTime <- proc.time() - ptm
    if (verbose) { 
      message(paste("runEM iter", i-1, " Log likelihood:", loglik[i])) 
      message("runEM iter", i-1, " Time: ", format(elapsedTime[3] / 60, digits = 2), " min.")
    }
    if ((abs(loglik[i] - loglik[i - 1]) / abs(loglik[i])) < likChangeConvergence){#} && loglik[i] > loglik[i - 1]) {
      message("runEM iter", i-1, " EM Converged")
      converged = 0
    }
    if (loglik[i] < loglik[i - 1]){
      message("LIKELIHOOD DECREASED!!!")
      #message("Using previous iteration ", i-2)
      i <- i - 1
      converged = 1
    }
  }
  
  if (converged) {
    # Perform one last round of E-step to get latest responsibilities
    #E-step
    if (verbose) { message("runEM iter", i-1 ,": Re-calculating responsibilties from converged parameters.") }   
    for (j in 1:length(chrsI)) {
      I <- chrsI[[j]]
      output <- .Call("forward_backward", piG[, i], A, py[, I], PACKAGE = "HMMcopy")
      rho[, I] <- output$rho
    }
  }
  
  if (verbose) {
    totalTime <- proc.time() - ptmTotal
    message("runEM: Using optimal parameters from iter", i-1)
    message("runEM: Total elapsed time: ", format(totalTime[3] / 60, digits = 2), "min.")
  }
  
  #### Return parameters ####
  n <- n[, 1:i, drop=FALSE]
  sp <- sp[, 1:i, drop=FALSE]
  phi <- phi[, 1:i, drop=FALSE]
  mus <- mus[, , 1:i, drop=FALSE]
  vars <- vars[, , 1:i, drop=FALSE]
  lambdas <- lambdas[, , 1:i, drop=FALSE]
  piG <- piG[, 1:i, drop=FALSE]
  loglik = loglik[1:i]
  
  output <- vector('list', 0);
  output$n <- n
  output$sp <- sp
  output$phi <- phi
  output$mus <- mus
  output$vars <- varsKS
  output$lambdas <- lambdas
  output$pi <- piG
  output$A <- A
  output$loglik <- loglik
  output$rho <- rho
  output$param <- param
  output$py <- py
  output$iter <- i
  
  return(output)
}

get2and3ComponentMixture <- function(ct, ct.sc.status, n, sp, phi){
  S <- length(n)
  cn <- 2
  mu <- matrix(NA, nrow = nrow(ct), ncol = ncol(ct))
  for (s in 1:S){
    #subclonal 3 component
    mu[ct.sc.status[, s] == TRUE, s] <- (((1 - n[s]) * (1 - sp[s]) * ct[ct.sc.status[, s] == TRUE, s]) + 
                                           ((1 - n[s]) * sp[s] * cn) + (n[s] * cn)) / ((1 - n[s]) * phi[s] + n[s] * cn)
    #clonal 2 component
    mu[ct.sc.status[ ,s] == FALSE, s] <- ((1 - n[s]) * ct[ct.sc.status[, s] == FALSE, s] + n[s] * cn) / ((1 - n[s]) * phi[s] + n[s] * cn)
  }
  return(log(mu))
}

get2ComponentMixture <- function(ct, n, phi){
  if (length(phi) > 1 && length(phi) != length(n)){
    stop("get2ComponentMixture: length(n) not equal to length(phi)")
  }
  S <- length(n)
  #if (length(phi) == 1){
  #  phi <- rep(phi, length(n))
  #}
  cn <- 2
  #mu <-  ((1 - n) * ct + n * cn) / ((1 - n) * phi + n * cn)
  mu <- NULL
  for (s in 1:S){
    mu <- cbind(mu, ((1 - n[s]) * ct[, s] + n[s] * cn) / ((1 - n[s]) * phi[s] + n[s] * cn))
  }
  return(log(mu))
}


# Dirichlet probability density function, returns the probability of vector
# x under the Dirichlet distribution with parameter vector alpha
# Author: David Ross http://www.cs.toronto.edu/~dross/code/dirichletpdf.m
dirichletpdf <- function(x, alpha) {
  if (any(x < 0)) {
    return(0);
  }
  if (abs(sum(x) - 1) > 1e-3) {
    stop("Dirichlet PDF: probabilities do not sum to 1")
    return(0);
  }
  p <- exp(lgamma(sum(alpha)) - sum(lgamma(alpha))) * prod(x ^ (alpha - 1))
  return(p)
}

dirichletpdflog <- function(x, k) {
  c <- lgamma(sum(k, na.rm = TRUE)) - sum(lgamma(k), na.rm = TRUE)  #normalizing constant
  l <- sum((k - 1) * log(x), na.rm = TRUE)  #likelihood
  y <- c + l
  return(y)
}

gammapdflog <- function(x, a, b) { #rate and scale parameterization
  c <- a * log(b) - lgamma(a)  # normalizing constant  
  l <- (a - 1) * log(x) + (-b * x)  #likelihood  
  y <- c + l
  return(y)
}

invgammapdflog <- function(x, a, b) {
  c <- a * log(b) - lgamma(a)
  l <- (-a - 1) * log(x) + (-b / x)
  y <- c + l
  return(y)
}

betapdflog <- function(x, a, b) {
  y = -lbeta(a, b) + (a - 1) * log(x) + (b - 1) * log(1 - x)
  return(y)
}

tdistPDF <- function(x, mu, lambda, nu) {
  S <- ncol(x)
  if (!is.null(S)){
    p <- NULL
    for (s in 1:S){
      tpdf <- (gamma(nu / 2 + 0.5)/gamma(nu / 2)) * ((lambda[s] / (pi * nu)) ^ (0.5)) *
        (1 + (lambda[s] * (x[, s] - mu[s]) ^ 2) / nu) ^ (-0.5 * nu - 0.5)
      p <- cbind(p, tpdf)
    }
  }else{
    nu <- rep(nu, length(mu))
    p <- (gamma(nu / 2 + 0.5)/gamma(nu / 2)) * ((lambda / (pi * nu)) ^ (0.5)) *
      (1 + (lambda * (x - mu) ^ 2) / nu) ^ (-0.5 * nu - 0.5)
  }
  p[is.na(p)] <- 1
  return(p)
}

normalpdf <- function(x, mu, var){
  # S <- ncol(x)
  # if (!is.null(S)){
  #   y <- NULL
  #   for (s in 1:S){
  #     c <- log(1/(sqrt(var[s]) * sqrt(2 * pi)))  #normalizing constant
  #     l <- -((x[, s] - mu[s])^2)/(2 * var[s])  #likelihood
  #     y <- cbind(y, c + l)
  #   }
  # }else{
    c <- log(1/(sqrt(var) * sqrt(2 * pi)))  #normalizing constant
    l <- -((x - mu)^2)/(2 * var)  #likelihood
    y <- c + l
  #}
  y[is.na(y)] <- 0
  y <- exp(y)
  return(y)
}

getNormLik <- function(copy, mus, varsKS, sw){
  S <- param$numberSamples
  K <- length(param$ct)
  Z <- sum(param$ct.sc) #number of subclonal states (# repeated copy number states)
  KS <- K ^ S
  py <- t(sapply(1:KS, function(ks) {
    y <- NULL
    for (s in 1:S){
      sl <- normalpdf(copy[, s], mus[ks, s, 1], varsKS[ks, s]) ^ sw[s]
      y <- cbind(y, sl)
    }
    probs <- apply(y, 1, prod)
  }))
  return(py)
}

mvGaussDistPDF <- function(x, mu, covar){
  # Multivariate Gaussian density function, return in exp scale #
  S <- ncol(x)
  T <- nrow(x)
  x <- t(as.matrix(x))
  mu <- matrix(mu, nrow = S, ncol = T, byrow = TRUE)
  c <- log(2 * pi) * (S/2) + log(det(covar)) * (1/2)
  l <-  1/2 * colSums(solve(covar) %*% (x - mu) * (x-mu))
  y <- -c - l
  y[is.na(y)] <- 0
  return(exp(y))
}

getCovarianceMatrix <- function(x){
  # sample covariance 
  covar <- cov(x, use = "pairwise.complete.obs")
  # correlation
  cor <- covar[1, 2] / sqrt(covar[1, 1] * covar[2, 2])
  return(list(covar = covar, cor = cor))
}

# Multivariate gamma function 
multiGammaFunlog <- function(p, a){
  if (p == 1){
    return(lgamma(a))
  }else{
    return(log(pi ^ ((p - 1) / 2) + multiGammaFunlog(p - 1, a) + multiGammaFunlog(1, a + (1 - p) / 2)))
  }
}

# Inverse Wishart density function in log scale 
invWishartpdflog <- function(covar, psi, nu){
  p <- ncol(covar)
  const <- log(det(psi)) * (nu / 2) - log(2) * (nu * p / 2) + multiGammaFunlog(p, nu / 2)
  lik <- log(det(covar)) * -((nu + p + 1) / 2) - 0.5 * sum(diag(psi %*% solve(covar)))
  y <- const + lik
  return(y)
}

estimateTParamsMap <- function(D, n_prev, sp_prev, phi_prev, lambda_prev, pi_prev, A_prev, 
                              params, rho, Zcounts, 
                              estimateNormal = TRUE, estimatePloidy = TRUE,
                              estimatePrecision = TRUE, estimateInitDist = TRUE, 
                              estimateTransition = TRUE, estimateSubclone = TRUE,
                              verbose = TRUE) {
  KS <- nrow(rho)
  K <- length(params$ct)
  T <- ncol(rho)
  S <- length(n_prev)
  #dimen <- ncol(D)
  intervalNormal <- c(1e-6, 1 - 1e-6)
  intervalSubclone <- c(1e-6, 1 - 1e-6)
  intervalPhi <- c(.Machine$double.eps, 10)
  intervalLambda <- c(1e-5, 1e4)
  
  # initialize params to be estimated
  n_hat <- n_prev
  sp_hat <- sp_prev
  phi_hat <- phi_prev
  lambda_hat <- lambda_prev
  pi_hat <- pi_prev
  A_hat <- A_prev
  F_lik <- 0
  # Update transition matrix A
  if (estimateTransition){
    for (k in 1:KS) {
      A_hat[k, ] <- Zcounts[k, ] + params$dirPrior[k, ]
      A_hat[k, ] <- normalize(A_hat[k, ])
    }
  }
  # map estimate for pi (initial state dist)
  if (estimateInitDist){
    pi_hat <- estimateMixWeightsParamMap(rho, params$kappa)
  }
  
  # map estimate for normal 
  if (estimateNormal){
    suppressWarnings(
      estNorm <- optim(n_prev, fn = completeLikelihoodFun, pType = rep("n", S), n = n_prev, sp = sp_prev, phi = phi_prev, 
                       lambda = lambda_prev, piG = pi_hat, A = A_hat,
                       params = params, D = D, rho = rho, Zcounts = Zcounts, 
                       estimateNormal = estimateNormal, estimatePloidy = estimatePloidy,
                       estimatePrecision = estimatePrecision, estimateInitDist = estimateInitDist,
                       estimateTransition = estimateTransition, estimateSubclone = estimateSubclone,
                       method = "L-BFGS-B",
                       lower = intervalNormal[1], upper = intervalNormal[2],
                       control = list(trace = 0, fnscale = -1))
    )
    n_hat <- estNorm$par
    F_lik <- estNorm$value
  }
  if (estimateSubclone){
    suppressWarnings(
      estSubclone <- optim(sp_prev, fn = completeLikelihoodFun, pType = rep("sp", S), n = n_hat, sp = sp_prev, 
                           phi = phi_prev, lambda = lambda_prev, piG = pi_hat, A = A_hat,
                           params = params, D = D, rho = rho, Zcounts = Zcounts, 
                           estimateNormal = estimateNormal, estimatePloidy = estimatePloidy,
                           estimatePrecision = estimatePrecision, estimateInitDist = estimateInitDist,
                           estimateTransition = estimateTransition, estimateSubclone = estimateSubclone,
                           method = "L-BFGS-B",
                           lower = intervalNormal[1], upper = intervalNormal[2],
                           control = list(trace = 0, fnscale = -1))
    )
    sp_hat <- estSubclone$par
    F_lik <- estSubclone$value
  }
  if (estimatePloidy){
    suppressWarnings(
      estPhi <- optim(phi_prev, fn = completeLikelihoodFun, pType = rep("phi", length(phi_prev)), n = n_hat, sp = sp_hat, 
                      phi = phi_prev, lambda = lambda_prev, piG = pi_hat, A = A_hat,
                      params = params, D = D, rho = rho, Zcounts = Zcounts, 
                      estimateNormal = estimateNormal, estimatePloidy = estimatePloidy,
                      estimatePrecision = estimatePrecision, estimateInitDist = estimateInitDist,
                      estimateTransition = estimateTransition, estimateSubclone = estimateSubclone,
                      method = "L-BFGS-B",
                      lower = intervalPhi[1], upper = intervalPhi[2],
                      control = list(trace = 0, fnscale = -1))
    )
    phi_hat <- estPhi$par
    F_lik <- estPhi$value
  }
  if (estimatePrecision){
    suppressWarnings(
      estLambda <- optim(c(lambda_prev), fn = completeLikelihoodFun, pType = rep("lambda", K*S), n = n_hat, sp = sp_hat,
                         phi = phi_hat, lambda = lambda_prev, piG = pi_hat, A = A_hat,
                         params = params, D = D, rho = rho, Zcounts = Zcounts, 
                         estimateNormal = estimateNormal, estimatePloidy = estimatePloidy,
                         estimatePrecision = estimatePrecision, estimateInitDist = estimateInitDist,
                         estimateTransition = estimateTransition, estimateSubclone = estimateSubclone,
                         method = "L-BFGS-B", lower = intervalLambda[1],
                         control = list(trace = 0, fnscale = -1))
    )
    lambda_hat <- matrix(estLambda$par, ncol = S, byrow = FALSE)
    F_lik <- estLambda$value
  }

  gc(verbose = FALSE, reset = TRUE)
  #if (verbose == TRUE) {
  #  message("n=", format(n_hat, digits = 2), ", phi=", format(phi_hat, digits = 4), ", lambda=", paste0(format(lambda_hat, digits=2), collapse = " ", sep = ","))
  #}
  return(list(n = n_hat, sp = sp_hat, phi = phi_hat, lambda = lambda_hat, piG = pi_hat, A = A_hat, F = F_lik))
}

estimateGaussianParamsMap <- function(D, n_prev, sp_prev, phi_prev, var_prev, pi_prev, A_prev, 
                              params, rho, Zcounts, 
                              estimateNormal = TRUE, estimatePloidy = TRUE,
                              estimateVar = TRUE, estimateInitDist = TRUE, 
                              estimateTransition = TRUE, estimateSubclone = TRUE,
                              verbose = TRUE) {
  KS <- nrow(rho)
  K <- length(params$ct)
  T <- ncol(rho)
  S <- length(n_prev)
  #dimen <- ncol(D)
  intervalNormal <- c(1e-6, 1 - 1e-6)
  intervalSubclone <- c(1e-6, 1 - 1e-6)
  intervalPhi <- c(.Machine$double.eps, 10)
  
  # initialize params to be estimated
  n_hat <- n_prev
  sp_hat <- sp_prev
  phi_hat <- phi_prev
  var_hat <- var_prev
  #varKS_prev <- as.matrix(expand.grid(as.data.frame(var_prev)))
  pi_hat <- pi_prev
  A_hat <- A_prev
  
  # Update transition matrix A
  if (estimateTransition){
    for (k in 1:KS) {
      A_hat[k, ] <- Zcounts[k, ] + params$dirPrior[k, ]
      A_hat[k, ] <- normalize(A_hat[k, ])
    }
  }
  # map estimate for pi (initial state dist)
  if (estimateInitDist){
    pi_hat <- estimateMixWeightsParamMap(rho, params$kappa)
  }
  
  # setup simplification variables for responsibilities and data
  a <- rowSums(rho)
  b <- rho %*% D
  c <- rho %*% (D^2)

  # map estimate for normal 
  if (estimateNormal){
    for (s in 1:S){
      funcN <- function(n) { nUpdateEqn(n, sp_prev[s], phi_prev[s], var_prev[, s], 
                                        params$jointStates[, s], params$jointCNstates[, s], params$jointSCstatus[, s], 
                                        params$alphaN, params$betaN, a, b[, s]) }
      n_hat[s] <- tryCatch({
        uniroot(funcN, intervalNormal, tol = 1e-15)$root
      }, error = function(x){ 
        message("ichorCNA updateParameters: Issue maximizing n, using previous iteration (", n_prev[s], ")")
        return(n_prev[s]) 
      })
    }
  }
  
  if (estimateSubclone){
    for (s in 1:S){
      funcS <- function(sp) { spUpdateEqn(sp, n_hat[s], phi_prev[s], var_prev[, s], 
                                         params$jointStates[, s], params$jointCNstates[, s], params$jointSCstatus[, s], 
                                         params$alphaSp, params$betaSp, a, b[, s]) }
      sp_hat[s] <- tryCatch({
        uniroot(funcS, intervalSubclone, tol = 1e-15)$root
      }, error = function(x){ 
        message("ichorCNA updateParameters: Issue maximizing sp, using previous iteration (", sp_prev[s], ")")
        return(sp_prev[s]) 
      })
    }
  } 

  if (estimatePloidy){
    for (s in 1:S){
      funcP <- function(phi) { phiUpdateEqn(phi, n_hat[s], sp_hat[s], var_prev[, s], 
                                            params$jointStates[, s], params$jointCNstates[, s], params$jointSCstatus[, s], 
                                            params$alphaSp, params$betaSp, a, b[, s]) }
      phi_hat[s] <- tryCatch({
        uniroot(funcP, intervalPhi, tol = 1e-15)$root
      }, error = function(x){ 
        message("ichorCNA updateParameters: Issue maximizing phi, using previous iteration (", phi_prev[s], ")")
        return(phi_prev[s]) 
      })
    }
  } 
  
  if (estimateVar){
    for (s in 1:S){
      var_hat[, s] <- varUpdateEqn(n_hat[s], sp_hat[s], phi_hat[s], params$jointStates[, s], 
                                 params$jointCNstates[, s], params$jointSCstatus[, s], 
                                 params$alphaVar[, s], params$betaVar[s], a, b[, s], c[, s])
    }
    #covar_hat <- covarIWUpdateEqn(n_hat, sp_hat, phi_hat, params, T, a, b, c)
  } 
  
  #gc(verbose = FALSE, reset = TRUE)
  #if (verbose == TRUE) {
  #  message("n=", format(n_hat, digits = 2), ", phi=", format(phi_hat, digits = 4), ", lambda=", paste0(format(lambda_hat, digits=2), collapse = " ", sep = ","))
  #}
  return(list(n = n_hat, sp = sp_hat, phi = phi_hat, var = var_hat, piG = pi_hat, A = A_hat))
}

# Update equation for normal parameter (Gaussian)
nUpdateEqn <- function(n, sp, phi, var, jointStates, jointCNstates, jointSCstatus, alphaN, betaN, a, b){
  mus <- as.matrix(get2and3ComponentMixture(data.frame(jointCNstates), data.frame(jointSCstatus), n, sp, phi))
  KS <- nrow(mus)
  S <- length(n)
  cn <- 2
  G <- jointStates
  K <- length(unique(G))
  
  #dQ_dn <- rep(0, KS)
  #for (k in 1:KS){ 
  #  dmu_dn <- (cn - sp * cn - (1 - sp) * ct[k]) / (n * cn + (1 - n) * sp * cn + (1 - n) * (1 - sp) * ct[k]) - 
  #    (cn - phi)/(n * cn + (1 - n) * phi)
  #  dQ_dn[k] <-  var[k]^-1 * (b[k] - a[k] * mus[k,1]) * dmu_dn 
  #}
  ind <- which(!jointSCstatus)
  ct <- jointCNstates[ind]
  dmu_dn <- (cn - ct) / (n * cn + (1 - n) * ct) - (cn - phi) / (n * cn + (1 - n) * phi)
  dQ_dn.noSC <- var[ind]^-1 * (b[ind] - a[ind] * mus[ind]) * dmu_dn
  
  ind <- which(jointSCstatus)
  ct.s <- jointCNstates[ind]
  if (length(ct.s) > 0){
    dmu_dn.s <- (cn - sp * cn - (1 - sp) * ct.s) / (n * cn + (1 - n) * sp * cn + (1 - n) * (1 - sp) * ct.s) - 
      (cn - phi)/(n * cn + (1 - n) * phi)
    dQ_dn.SC <- var[ind]^-1 * (b[ind] - a[ind] * mus[ind]) * dmu_dn.s
  }else{
    dQ_dn.SC <- 0
  }
  dQ_dn <- sum(c(dQ_dn.noSC, dQ_dn.SC))
  dbeta_dn <- (alphaN - 1) / n - (betaN - 1) / (1 - n)
  f <- dQ_dn + dbeta_dn
  return(f)
}

# Update equation for cellular prevalence parameter (Gaussian)
spUpdateEqn <- function(sp, n, phi, var, jointStates, jointCNstates, jointSCstatus, alphaSp, betaSp, a, b){
  mus <- as.matrix(get2and3ComponentMixture(data.frame(jointCNstates), data.frame(jointSCstatus), n, sp, phi))
  KS <- nrow(mus)
  S <- length(n)
  cn <- 2
  ct <- jointCNstates
  dmu_ds <- ((1 - n) * cn - (1 - n) * ct) / (n * cn + (1 - n) * sp * cn + (1 - n) * (1 - sp) * ct)
  #dQ_dmu <- b %*%  - a * mus %*% solve(covar)
  dQ_ds <- var^-1 * (b - a * mus) * dmu_ds
  dbeta_ds <- (alphaSp - 1) / sp - (betaSp - 1) / (1 - sp)
  f <- sum(dQ_ds) + dbeta_ds
  return(f)
}

# Update equation for tumor ploidy parameter (Gaussian)
phiUpdateEqn <- function(phi, n, sp, var, jointStates, jointCNstates, jointSCstatus, alphaPhi, betaPhi, a, b){
  mus <- as.matrix(get2and3ComponentMixture(data.frame(jointCNstates), data.frame(jointSCstatus), n, sp, phi))
  KS <- nrow(mus)
  S <- length(n)
  cn <- 2
  ct <- jointCNstates
  dmu_dphi <- -(1 - n) / (n * cn + (1 - n) * phi)
  #dQ_dmu <- b %*%  - a * mus %*% solve(covar)
  dQ_dphi <- var^-1 * (b - a * mus) * dmu_dphi
  dgamma_dphi <- (alphaPhi - 1) / phi - (betaPhi - 1) / (1 - phi)
  f <- sum(dQ_dphi) + dgamma_dphi
  return(f)
}

# Update equation for variance parameter (MVN)
varUpdateEqn <- function(n, sp, phi, jointStates, jointCNstates, jointSCstatus, alphaVar, betaVar, a, b, c){
  mus <- as.matrix(get2and3ComponentMixture(data.frame(jointCNstates), data.frame(jointSCstatus), n, sp, phi))
  KS <- nrow(mus)
  K <- nrow(param$var)
  S <- length(n)
  cn <- 2
  ct <- jointCNstates
  # consolidate to K states from K^S expand states
  # term1 <- rep(0, K)
  # term2 <- rep(0, K)
  # for (i in 1:K){
  #   k <- jointStates == i
  #   term1[i] <- -sum((c[k] - 2*b[k]*mus[k, 1] + a[k]*mus[k, 1]^2)) 
  #   term2[i] <- -sum(a[k])
  # }
  # term1 <- term1 - 2 * betaVar
  # term2 <- term2 - 2 * (alphaVar + 1)
  term1 <- c - 2*b*mus + a*mus^2 + 2*betaVar
  term2 <- a + 2*(alphaVar + 1)
  var_hat <- term1 / term2
  return(var_hat)
}

# Update equation for covariance parameter (Gaussian)
covarUpdateEqn <- function(covar, n, sp, phi, params, a, b, c){
  mus <- as.matrix(get2and3ComponentMixture(params$jointCNstates, params$jointSCstatus, n, sp, phi))
  KS <- nrow(mus)
  S <- length(n)
  cn <- 2
  ct <- params$jointCNstates
  dQ_dcovar <- 0
  for (ks in 1:KS){
    tmp <- a[ks] * solve(covar) - solve(covar) * 
      (c[ks, ] - b[ks, ] * mus[ks, ] - mus[ks, ] * b[ks, ] + a[ks] * mus[ks, ] * mus[ks, ]) * solve(covar)
    dQ_dcovar <- dQ_dcovar + tmp
  }
  dig_dcovar <- ((params$nu + S + 1) / 2) * sum(diag(covar)) + 
    ((params$psi) / 2) * sum(diag(solve(covar) * solve(covar)))
  f <- dQ_dcovar - dig_dcovar
  return(f)
}

# Update equation for covariance parameter (MVN + IW)
covarIWUpdateEqn <- function(n, sp, phi, params, T, a, b, c){
  mus <- as.matrix(get2and3ComponentMixture(params$jointCNstates, params$jointSCstatus, n, sp, phi))
  KS <- nrow(mus)
  S <- length(n)
  cn <- 2
  ct <- params$jointCNstates
  covar_k <- (c[k] - 2*b[k]*mus[k, ] + a[k]*mus[k, ]^2) / T
  covar_hat <- (b[k]*params$psi + T*covar_k) / (params$nu + T + S + 1)
  return(covar_hat)
}

# Student's t likelihood function #
stLikelihood <- function(n, sp, phi, lambda, params, D, rho){
  KS <- nrow(rho)
  lik <- 0
  # Recalculate the likelihood
  lambdaKS <- as.matrix(expand.grid(as.data.frame(lambda)))
  mus <- as.matrix(get2and3ComponentMixture(params$jointCNstates, params$jointSCstatus, n, sp, phi))
  lik <- sapply(1:KS, function(ks){
    probs <- log(tdistPDF(D, mus[ks, ], lambdaKS[ks, ], params$nu))
    # multiply across samples for each data point to get joint likelihood.
    l <- rho[ks, ] %*% rowSums(as.data.frame(probs)) 
  })
  return(sum(lik))
}

## length of x will depend on which pararmeters are being estimated
## x values should be same as n, phi, lambda, pi, but may not necessarily
## include all of those variables
completeLikelihoodFun <- function(x, pType, n, sp, phi, lambda, piG, A, params, D, rho, Zcounts,
                                  estimateNormal = TRUE, estimatePloidy = TRUE,
                                  estimatePrecision = TRUE, estimateTransition = FALSE,
                                  estimateInitDist = TRUE, estimateSubclone = TRUE){
  KS <- nrow(rho)
  N <- ncol(rho)
  S <- params$numberSamples
  K <- length(params$ct)
  
  lambda <- matrix(lambda, ncol = S, byrow = FALSE)
  
  if (estimatePrecision && sum(pType == "lambda") > 0){
    lambda <- matrix(x[pType == "lambda"], ncol = S, byrow = FALSE)
  } 
  if (estimateNormal && sum(pType == "n") > 0){
    n <- x[pType == "n"]
  }  
  if (estimateSubclone && sum(pType == "sp") > 0){
    sp <- x[pType == "sp"]
  }  
  if (estimatePloidy && sum(pType == "phi") > 0){
    phi <- x[pType == "phi"]
  }
  
  ## prior probabilities ##
  prior <- priorProbs(n, sp, phi, lambda, var = NA, piG, A, params, 
                      estimateNormal = estimateNormal, estimatePloidy = estimatePloidy,
                      estimateVar = estimateVar, estimatePrecision = estimatePrecision, 
                      estimateTransition = estimateTransition, estimateInitDist = estimateInitDist, 
                      estimateSubclone = estimateSubclone)
  ## likelihood terms ##
  likObs <- stLikelihood(n, sp, phi, lambda, params, D, rho)
  likInit <- rho[, 1] %*% log(piG)
  likTxn <- sapply(1:KS, function(ks){
    Zcounts[ks, ] %*% log(A[ks, ])
  })
 
  ## sum together ##
  lik <- likObs + likInit + sum(likTxn)
  prior <- prior$prior
  f <- lik + prior
  return(f)
}

priorProbs <- function(n, sp, phi, lambda, var, piG, A, params, 
                       estimateNormal = TRUE, estimatePloidy = TRUE,
                       estimateVar = TRUE, estimatePrecision = TRUE, 
                       estimateTransition = TRUE, estimateInitDist = TRUE, 
                       estimateSubclone = TRUE){
  S <- params$numberSamples
  K <- length(params$ct)
  KS <- K ^ S
  ## prior terms ##
  priorA <- rep(0, KS)
  if (estimateTransition){
    priorA <- sapply(1:KS, function(ks){
      dirichletpdflog(A[ks, ], params$dirPrior[ks, ])
    })
  }
  priorA <- sum(priorA)
  priorN <- 0
  if (estimateNormal){
    priorN <- sum(betapdflog(n, params$alphaN, params$betaN))
  }
  priorSP <- 0
  if (estimateSubclone){
    priorSP <- sum(betapdflog(sp, params$alphaSp, params$betaSp))
  }
  priorPhi <- 0
  if (estimatePloidy){
    priorPhi <- sum(gammapdflog(phi, params$alphaPhi, params$betaPhi))
  }
  priorLambda <- 0
  priorVar <- 0
  if (params$likModel == "Gaussian"){
    if (estimateVar){
      for (s in 1:S){
        priorVar <- priorVar + 
          sum(invgammapdflog(as.data.frame(var[, s]), params$alphaVar[, s], params$betaVar[s]))
      }
    }
  }else{ # param$likModel == "t"
    if (estimatePrecision){
      for (s in 1:S){
        priorLambda <- priorLambda + 
          sum(gammapdflog(as.data.frame(lambda)[, s], params$alphaLambda, params$betaLambda[, s]))
      }
      priorLambda <- as.numeric(priorLambda)
    }
  }
  priorPi <- 0
  if (estimateInitDist){
    priorPi <- dirichletpdflog(piG, params$kappa)
  }
  prior <- priorA + priorLambda + priorVar + priorN + priorSP + priorPhi + priorPi  
  return(list(prior = prior, priorA = priorA, priorLambda = priorLambda, 
              priorN = priorN, priorSP = priorSP, priorPhi = priorPhi, priorPi = priorPi))
}

estimatePrecisionParamMap <- function(lambda, n, phi, params, D, rho){
  mu <- get2ComponentMixture(params$ct, n, phi)
  mu_0 <- get2ComponentMixture(params$ct, params$n_0, params$phi_0)
  yr <- t(matrix(D, length(D), length(mu)))        # Vectorize parameter
  u <- (1 + params$nu) / (((yr - mu) ^ 2) * lambda + params$nu); # scale parameter
  
  # Calculate the precision
  lambda_hat <- (rowSums(rho, na.rm = TRUE) + params$alphaLambda + 1) /
    (rowSums(rho * u * ((yr - mu) ^ 2), na.rm = TRUE) +
       (params$eta * (mu - mu_0) ^ 2 + params$betaLambda))
  return(lambda_hat)
}

estimateMixWeightsParamMap <- function(rho, kappa) {
  K <- nrow(rho)
  #pi <- (rho[, 1] + kappa - 1) / (sum(rho[, 1]) + sum(kappa) - K)
  pi <- (rowSums(rho, na.rm = TRUE) + kappa - 1) / 
    #	(ncol(rho) + sum(kappa) - K)
    (sum(rowSums(rho, na.rm = TRUE)) + sum(kappa) - K)
  
  return(pi)
} 

