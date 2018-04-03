# file:   EM.R
# author: Gavin Ha, Ph.D.
#         Justin Rhoades
#               Dana-Farber Cancer Institute
#               Broad Institute
# contact: <gavinha@broadinstitute.org>
# ULP-WGS website: http://www.broadinstitute.org/~gavinha/ULP-WGS/
# HMMcopy website: http://compbio.bccrc.ca/software/hmmcopy/ and https://www.bioconductor.org/packages/release/bioc/html/HMMcopy.html
# date:   Oct 26, 2016
# description: Hidden Markov model (HMM) to analyze Ultra-low pass whole genome sequencing (ULP-WGS) data.
# This script is the main script to run the HMM.

runEM <- function(copy, chr, chrTrain, param, maxiter, verbose = TRUE, 
									estimateNormal = TRUE, estimatePloidy = TRUE, estimatePrecision = TRUE,
                  estimateTransition = TRUE, estimateInitDist = TRUE, estimateSubclone = TRUE,
                  likChangeConvergence = 1e-3) {
  
  if (nrow(copy) != length(chr) || nrow(copy) != length(chrTrain)) {
    stop("runEM: Length of inputs do not match for one of: copy, chr, chrTrain")
  }
  
  if (is.null(param$ct) || is.null(param$lambda) || is.null(param$nu) ||
      is.null(param$kappa)) {
    stop("runEM: Parameter missing, ensure all parameters exist as columns in",
         "data frame: ct, lambda, nu, kappa")
  }
  
  S <- param$numberSamples
  K <- length(param$ct)
  Z <- sum(param$ct.sc) #number of subclonal states (# repeated copy number states)
  KS <- K ^ S
  N <- nrow(copy)
  rho <- matrix(0, KS, N)
  py <- matrix(0, KS, N)            # Local evidence
  mus <- array(0, dim = c(KS, S, maxiter))     # State means
  lambdas <- array(0, dim = c(K, S, maxiter))  # State Variances
  phi <- matrix(NA, length(param$phi_0), maxiter)					 # Tumour ploidy
  n <- matrix(NA, S, maxiter)						 # Normal contamination
  sp <- matrix(NA, S, maxiter)     # cellular prevalence (tumor does not contain event)
  piG <- matrix(0, KS, maxiter)     # Initial state distribution
  converged <- FALSE               # Flag for convergence
  Zcounts <- matrix(0, KS, KS)
  loglik <- rep(0, maxiter)
  
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
  mus[, , i] <- as.matrix(get2and3ComponentMixture(param$jointCNstates, param$jointSCstatus, n[, i], sp[, i], phi[, i]))
  
  # Likelihood #
  for (ks in 1:KS) {
    probs <- tdistPDF(copy, mus[ks, , i], lambdasKS[ks, ], param$nu)
    py[ks, ] <- apply(probs, 1, prod) # multiply across samples for each data point to get joint likelihood.
  }
  
  # initialize transition prior #
  A <- normalize(param$A)
  A_prior <- A
  dirPrior <- param$dirPrior
  
  loglik[i] <- -Inf
  
  while(!converged && (i < maxiter)) {
    ptm <- proc.time()
    #if (verbose) { message(paste("runEM: EM iteration:", i, "Log likelihood:", loglik[i])) }
    i <- i + 1
    
    ################ E-step ####################
    if (verbose) { message("runEM iter", i-1 , ": Expectation") }
    Zcounts <- matrix(0, KS, KS)
    for (j in 1:length(chrsI)) {
      I <- intersect(chrsI[[j]], which(chrTrain))
      if (length(I) > 0){
        output <- .Call("forward_backward", piG[, i - 1], A, py[, I],PACKAGE = "HMMcopy")
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
    #    rhoS[[s]][k, ] <- colSums(rho[param$jointCNstates[, s] == k, chrTrain])
    #  }
    #}
    
    ################ M-step ####################
    if (verbose) { message("runEM iter", i-1 , ": Maximization") }
    output <- estimateParamsMap(copy[chrTrain, ], n[, i - 1], sp[, i - 1], phi[, i - 1], 
                                lambdas[, , i - 1], piG[, i - 1], A, param,
                                rho[, chrTrain], Zcounts, 
                                estimateNormal = estimateNormal, estimatePloidy = estimatePloidy,
                                estimatePrecision = estimatePrecision, estimateTransition = estimateTransition,
                                estimateInitDist = estimateInitDist, estimateSubclone = estimateSubclone)
    if (verbose == TRUE) {
      for (s in 1:S){
        message("Sample", s, " n=", signif(output$n[s], digits = 4), ", sp=", signif(output$sp[s], digits = 4), 
                ", phi=", signif(output$phi[s], digits = 4), 
                ", lambda=", paste0(signif(output$lambda[, s], digits=5), collapse = ","),
                ", F=", signif(output$F, digits=4))
      }     
    }
    n[, i] <- output$n
    sp[, i] <- output$sp
    phi[, i] <- output$phi
    lambdas[, , i] <- output$lambda
    piG[, i] <- output$piG
    A <- output$A
    estF <- output$F
    
    # Recalculate the likelihood
    lambdasKS <- as.matrix(expand.grid(as.data.frame(lambdas[, , i])))
    mus[, , i] <- as.matrix(get2and3ComponentMixture(param$jointCNstates, param$jointSCstatus, n[, i], sp[, i], phi[, i]))
    for (ks in 1:KS) {
      probs <- tdistPDF(copy, mus[ks, , i], lambdasKS[ks, ], param$nu)
      py[ks, ] <- apply(probs, 1, prod) # multiply across samples for each data point to get joint likelihood.
    }
    
    prior <- priorProbs(n[, i], sp[, i], phi[, i], lambdas[, , i], piG[, i], A, param, 
                        estimateNormal = estimateNormal, estimatePloidy = estimatePloidy,
                        estimatePrecision = estimatePrecision, estimateTransition = estimateTransition,
                        estimateInitDist = estimateInitDist, estimateSubclone = estimateSubclone)
    
    
    # check converence 
    loglik[i] <- loglik[i] + prior$prior
    elapsedTime <- proc.time() - ptm
    if (verbose) { 
      message(paste("runEM iter", i-1, " Log likelihood:", loglik[i])) 
      message("runEM iter", i-1, " Time: ", format(elapsedTime[3] / 60, digits = 2), " min.")
    }
    if ((abs(loglik[i] - loglik[i - 1]) / abs(loglik[i])) < likChangeConvergence){#} && loglik[i] > loglik[i - 1]) {
      converged = 1
    }
    if (loglik[i] < loglik[i - 1]){
      #message("LIKELIHOOD DECREASED!!!")
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
  lambdas <- lambdas[, , 1:i, drop=FALSE]
  piG <- piG[, 1:i, drop=FALSE]
  loglik = loglik[1:i]
  
  output <- vector('list', 0);
  output$n <- n
  output$sp <- sp
  output$phi <- phi
  output$mus <- mus
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


estimateParamsMap <- function(D, n_prev, sp_prev, phi_prev, lambda_prev, pi_prev, A_prev, 
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
  }
  
  #}
  #}
  
  #map estimate for phi (ploidy) 
  # if (estimatePloidy) {
  #    suppressWarnings(
  #   estPhi <- optim(phi_prev, fn = phiLikelihoodFun, n = n_prev, lambda = lambda_prev, 
  #  								params = params, D = D, rho = rhoS, method = "L-BFGS-B",
  #                 control = list(fnscale = -1, ndeps = 1e-5), lower = intervalPhi[1])#, upper = intervalPhi[2])
  #  )
  #  phi_hat <- estPhi$par
  #}
  
  # map estimate for lambda (Student's t precision)
  # if (estimatePrecision){
  #for (s in 1:S){
  #  paramsS <- params
  #  paramsS$betaLambda <- paramsS$betaLambda[s]
  #  for (k in 1:K){
  #		  suppressWarnings(
  #	  estLambda <- optim(lambda_prev[k, s], fn = lambdaLikelihoodFun, n = n_prev[s], phi = phi_prev, 
  #										params = paramsS, D = D[, s], rho = rhoS[[s]][k, , drop = FALSE], k = k,
  #										control = list(fnscale = -1), method = "L-BFGS-B",
  #										lower = intervalLambda[1])#, upper = intervalLambda[2])
  
  
  #lambda_hat <- matrix(estLambda$par, ncol = S, byrow = TRUE)
  #  }
  #}
  #}
  #lambda_hat <- estimatePrecisionParamMap(lambda_prev, n_prev, phi_prev, params, D, rho)
  
  
  #rm(a, b, c, d, e)
  gc(verbose = FALSE, reset = TRUE)
  #if (verbose == TRUE) {
  #  message("n=", format(n_hat, digits = 2), ", phi=", format(phi_hat, digits = 4), ", lambda=", paste0(format(lambda_hat, digits=2), collapse = " ", sep = ","))
  #}
  return(list(n = n_hat, sp = sp_hat, phi = phi_hat, lambda = lambda_hat, piG = pi_hat, A = A_hat, F = estLambda$value))
}

# Student's t likelihood function #
stLikelihood <- function(n, sp, phi, lambda, params, D, rho){
  KS <- nrow(rho)
  lik <- 0
  # Recalculate the likelihood
  lambdaKS <- as.matrix(expand.grid(as.data.frame(lambda)))
  mus <- as.matrix(get2and3ComponentMixture(params$jointCNstates, params$jointSCstatus, n, sp, phi))
  for (ks in 1:KS) {
    probs <- log(tdistPDF(D, mus[ks, ], lambdaKS[ks, ], params$nu))
    # multiply across samples for each data point to get joint likelihood.
    l <- rho[ks, ] %*% rowSums(as.data.frame(probs)) 
    lik <- lik + as.numeric(l)
  }
  return(lik)
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
  prior <- priorProbs(n, sp, phi, lambda, piG, A, params, 
                      estimateNormal = estimateNormal, estimatePloidy = estimatePloidy,
                      estimatePrecision = estimatePrecision, estimateTransition = estimateTransition,
                      estimateInitDist = estimateInitDist, estimateSubclone = estimateSubclone)
  ## likelihood terms ##
  likObs <- stLikelihood(n, sp, phi, lambda, params, D, rho)
  likInit <- rho[, 1] %*% log(piG)
  likTxn <- 0
  for (ks in 1:KS){
    likTxn <- likTxn + Zcounts[ks, ] %*% log(A[ks, ])
  }
  
  ## sum together ##
  lik <- likObs + likInit + likTxn
  prior <- prior$prior
  f <- lik + prior
  return(f)
}

priorProbs <- function(n, sp, phi, lambda, piG, A, params, 
                       estimateNormal = TRUE, estimatePloidy = TRUE,
                       estimatePrecision = TRUE, estimateTransition = TRUE,
                       estimateInitDist = TRUE, estimateSubclone = TRUE){
  S <- params$numberSamples
  K <- length(params$ct)
  KS <- K ^ S
  ## prior terms ##
  priorA <- 0
  if (estimateTransition){
    for (ks in 1:KS){
      priorA <- priorA + dirichletpdflog(A[ks, ], params$dirPrior[ks, ])
    }
  }
  priorLambda <- 0
  if (estimatePrecision){
    for (s in 1:S){
      for (k in 1:K){
        priorLambda <- priorLambda + 
          gammapdflog(as.data.frame(lambda)[k, s], params$alphaLambda[k], params$betaLambda[k, s])
      }
    } 
  }
  priorLambda <- as.numeric(priorLambda)
  priorN <- 0
  if (estimateNormal){
    priorN <- sum(betapdflog(n, params$alphaN, params$betaN))
  }
  priorSP <- 0
  if (estimateNormal){
    priorSP <- sum(betapdflog(sp, params$alphaSp, params$betaSp))
  }
  priorPhi <- 0
  if (estimatePloidy){
    priorPhi <- sum(gammapdflog(phi, params$alphaPhi, params$betaPhi))
  }
  priorPi <- 0
  if (estimateInitDist){
    priorPi <- dirichletpdflog(piG, params$kappa)
  }
  prior <- priorA + priorLambda + priorN + priorSP + priorPhi + priorPi  
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

