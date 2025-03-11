library(glmnet)

# Load the data
#
#   Parameters:
# filename - string, name of the file containing the data.
#
#   Return value:
# A data frame with categorical variables represented as factors. Response
# variable is a factor 'disoutcome' taking on 3 different values:
# 1 - survival without sequelae;
# 2 - survival with sequelae;
# 3 - death.
fLoadData <- function(filename) {
  d <- read.csv(filename, header = TRUE, row.names = 1)
  # d$comasc <- factor(d$comasc) # TODO: should it be numeric, ordinal or factor?
  if(!is.null(d$admretin)) d$admretin <- factor(d$admretin)
  d$conv <- factor(d$conv)
  d$glu <- factor(d$glu)
  if(!is.null(d$hiv_indicator)) d$hiv_indicator <- factor(d$hiv_indicator)
  if(!is.null(d$hiv_ind)) {
    d$hiv_indicator <- factor(d$hiv_ind)
    d$hiv_ind <- NULL
  }
  d$disoutcome <- ordered(d$disoutcome)
  attr(d, "units") = c(age = "month", wt = "kg", ht = "cm", admgluc = "mmol/L",
                       comasc = "", disoutcome = "", admpta = "?", abldwbc = "?",
                       platelet = "?", comadur = "hour", admretin = "",
                       conv = "", glu = "", hiv_indicator = "",
                       Tmax = "degree C", lactate = "?", edema = "?", HRP2 = "?")
  d
}

# Sigmoid function f(x) = 1 / (1 + exp(-x))
#
#   Parameters:
# x - value, matrix-like.
#
#   Return value:
# Value of the sigmoid function at x, has the same shape as x.
sigmoid <- function(x) {
  1 / (1 + exp(-x))
}

# Inverse of the sigmoid function f(x) = 1 / (1 + exp(-x)): f(y) = ln(y / (1 - y))
#
#   Parameters:
# x - value, matrix-like.
#
#   Return value:
# Value of the inverse sigmoid function at x, has the same shape as x.
sigmoidInverse <- function(x) {
  log(x / (1 - x))
}

dropAttributes <- function(object, attributes) {
  res <- object
  for(a in attributes) {
    attr(res, a) <- NULL
  }
  res
}

############################################################################
# Baseline Category Logit Model
#
# Logits are given (for a 3-class problem) as
# \log(\frac{\pi_1}{\pi_3}) = \alpha_1 + \beta_1^Tx = L_1(x)
# \log(\frac{\pi_2}{\pi_3}) = \alpha_2 + \beta_2^Tx = L_2(x)
# This is a standard multinomial regression model implemented in statistical
# software packages.
#    For the model that satisfies the prevailing theory, this model has a form
# \log(\frac{\pi_1}{\pi_3}) = \alpha_1 + \beta^Tx = L_1(x)
# \log(\frac{\pi_2}{\pi_3}) = \alpha_2 + \gamma \beta^Tx = L_2(x)
# where \gamma \in (0, 1)
#   We use vglm function in VGAM package to estimate unrestricted model.
#   For the restricted model that satisfies the prevailing theory we use
# our own implementation.
#
#   Parameters:
# formula - a formula to be passed to the vglm function, describing the model
#           to estimate
# data - a data frame with the data to use for model estimation
# restricted - boolean, indicator if we want to fit a restricted model (TRUE),
#              the one that satisfies the prevailing theory assumption or
#              unrestricted one (FALSE).
# maxit - maximal number of iterations in the optim function (L-BFGS-B
#         algorithm); ignored if restricted = FALSE.
#
#   Return value:
# An object returned by VGAM::vglm function.
bclm <- function(formula, data, restricted = TRUE, maxit = 1000) {
  if(restricted) {
    #cl is not really needed, we only adds this call to the return object
    cl <- match.call()
    res <- bclmRestricted(formula, data, maxit = maxit)
    res$call = cl
    return(res)
  } else {
    VGAM::vglm(formula = formula, data = data, family = VGAM::multinomial)
  }
}

### functions used in optimization in bclmRestricted

# Compute values of the logits L1(x) and L2(x), see bclm function description, restricted model.
# 
#   Parameters:
# b - a vector of parameters, first element corresponds to \alpha_1,
#     second element - to \alpha_2, third element - to \gamma, and all the rest
#     to \beta.
# x - design matrix of  size n x p, where p is the number of linear predictors.
#
#   Return value:
# A matrix of size nobs x 2 of computed logits. Also has two attributes:
# betax: betax = x %*% beta
# params: list of alpha1, alpha2, gamma and beta
bclmLogits <- function(b, x) {
  alpha1 <- b[1]
  alpha2 <- b[2]
  gamma <- b[3]
  beta <- b[-(1:3)]
  betax <- as.vector(x %*% beta)
  L1 <- alpha1 + betax
  L2 <- alpha2 + gamma * betax
  res <- cbind(L1, L2)
  colnames(res) <- c("log(mu[, 1] / mu[, 3])",  "log(mu[, 2] / mu[, 3])")
  attr(res, "betax") <- betax
  attr(res, "params") <- list(alpha1 = alpha1, alpha2 = alpha2, gamma = gamma, beta = beta)
  res
  
}

# Compute probabilities for the baseline category model, see bclm function
# description, restricted model.
# 
#   Parameters:
# b - a vector of parameters, first element corresponds to \alpha_1,
#     second element - to \alpha_2, third element - to \gamma, and all the rest
#     to \beta.
# x - design matrix of  size n x p, where p is the number of linear predictors.
#
#   Return value:
# A matrix of size nobs x 3 of computed probabilities. Also has two attributes:
# betax: betax = x %*% beta
# params: list of alpha1, alpha2, gamma and beta
bclmProbabilities <- function(b, x) {
  L <- bclmLogits(b, x)
  v1 <- exp(L[, 1])
  v2 <- exp(L[, 2])
  p1 <- v1 / (v1 + v2 + 1)
  p2 <- v2 / (v1 + v2 + 1)
  p3 <- 1 / (v1 + v2 + 1)
  res <- cbind(p1, p2, p3)
  colnames(res) <- 1:3
  attr(res, "betax") <- attr(L, "betax")
  attr(res, "params") <- attr(L, "params")
  res
}

# Log-likelihood function.
#
#   Parameters:
# b - a vector of parameters, first element corresponds to \alpha_1,
#     second element - to \alpha_2, third element - to \gamma, and all the rest
#     to \beta.
# x - design matrix of  size n x p, where p is the number of linear predictors.
# y - vector of responses, assume dto be integers in \{1, 2, 3\}.
#
#   Return value:
# Value of the log-likelihood.
bclmF <- function(b, x, y) {
  p <- bclmProbabilities(b, x)
  sum((y == 1) * log(p[, 1]) + (y == 2) * log(p[, 2]) + (y == 3) * log(p[, 3]))
}

# Gradient of the log-likelihood function.
#
#   Parameters:
# b - a vector of parameters, first element corresponds to \alpha_1,
#     second element - to \alpha_2, third element - to \gamma, and all the rest
#     to \beta.
# x - design matrix of  size n x p, where p is the number of linear predictors.
# y - vector of responses, assume dto be integers in \{1, 2, 3\}.
#
#   Return value:
# Vector of gradients for the log-likelihood function.
bclmG <- function(b, x, y) {
  p <- bclmProbabilities(b, x)
  betax <- attr(p, "betax")
  gamma <- attr(p, "params")$gamma
  res <- c(sum((y == 1) - p[, 1]),
           sum((y == 2) - p[, 2]),
           sum(betax * ((y == 2) - p[, 2])),
           colSums(x * ((y == 1) - p[, 1] + gamma * ((y == 2) - p[, 2]))))
  names(res) <- names(b)
  res
}

# Second derivative (Hessian) of the log-likelihood function.
#
#   Parameters:
# b - a vector of parameters, first element corresponds to \alpha_1,
#     second element - to \alpha_2, third element - to \gamma, and all the rest
#     to \beta.
# x - design matrix of  size n x p, where p is the number of linear predictors.
# y - vector of responses, assume dto be integers in \{1, 2, 3\}.
#
#   Return value:
# Matrix of second derivatives for the log-likelihood function.
bclmH <- function(b, x, y) {
  alpha1 <- b[1]
  alpha2 <- b[2]
  gamma <- b[3]
  beta <- b[-(1:3)]
  betax <- as.vector(x %*% beta)
  L1 <- alpha1 + betax
  L2 <- alpha2 + gamma * betax
  v1 <- exp(L1)
  v2 <- exp(L2)
  p1 <- v1 / (v1 + v2 + 1)
  p2 <- v2 / (v1 + v2 + 1)
  p3 <- 1 / (v1 + v2 + 1)
  a11 <- -sum(p1 * (1 - p1))
  a12 <- sum(p1 * p2)
  a13 <- sum(p1 * p2 * betax)
  a1n <- -colSums(x * (p1 * (1 - p1 - gamma * p2)))
  a22 <- -sum(p2 * (1 - p2))
  a23 <- -sum(p2 * (1 - p2) * betax)
  a2n <- -colSums(x * (p2 * (gamma - p1 - gamma * p2)))
  a33 <- -sum(betax * betax * p2 * (1 - p2))
  a3n <- colSums(x * ((y == 2) - p2 - betax * p2 * (gamma - p1 - gamma * p2)))
  ann <- -crossprod(x * sqrt(p1 * (1 - p1 - gamma * p2) + gamma * p2 * (gamma - p1 - gamma * p2)))
  a <- rbind(c(a11, a12, a13, a1n),
             c(a12, a22, a23, a2n),
             c(a13, a23, a33, a3n))
  res <- rbind(a,
               cbind(t(a[, -(1:3)]), ann))
  colnames(res) <- names(b)
  rownames(res) <- names(b)
  res
}

# Helper function for estimation of the restricted model
bclmRestricted <- function(formula, data, weights = NULL, maxit = 1000, model = TRUE) {
  # prepare a call for actual non-standard evaluation
  mf <- match.call(expand.dots = FALSE)
  #message("Call with ellipses not expanded: ")
  #note that there are no ellipses in the function arguments for now, 
  #but you might want to change that later
  #print(mf)
  #turn weights into symbol if character is passed
  if (is.character(mf$weights)) mf$weights <- as.symbol(mf$weights)
  m <- match(c("formula", "data", "weights"), names(mf), 0L)
  #message("Position of formula, data, weights and maxit in the call:")
  #print(m)
  mf <- mf[c(1L, m)]
  #message("New call that only contains what is needed:")
  #print(mf)
  mf$drop.unused.levels <- TRUE 
  #message("Call with argument added:")
  #print(mf)
  mf[[1L]] <- quote(stats::model.frame) 
  #message("Change call to a call to model.frame:")
  #print(mf)
  mf <- eval(mf, parent.frame()) #evaluate call
  wt <- as.vector(model.weights(mf))
  stopifnot(is.null(wt)) # we do not support non-uniform weights
  y <- model.response(mf) # this is the response
  stopifnot(length(levels(y)) == 3) # only 3-level response is supported
  y <- as.numeric(y) # we convert response to numeric values
  mt <- attr(mf, "terms")
  x <- model.matrix(mt, mf) # this is the design matrix, it may include an intercept
  stopifnot(attr(mt, "intercept") == 1)
  if(attr(mt, "intercept") == 1) { # we drop the intercept column - it will be implicitly used
    x <- x[, colnames(x) != "(Intercept)", drop = FALSE]
  }
  b <- c(0, 0, 0.5, rep(0, ncol(x))) # initialize the parameter vector
  names(b) <- c("(Intercept):1", "(Intercept):2", "gamma", colnames(x))
  ll <- rep(-Inf, length(b))
  ul <- rep(Inf, length(b))
  ll[3] <- 0
  ul[3] <- 1
  oresG <- optim(par = b, fn = bclmF, gr = bclmG,
                  x = x, y = y,
                  method = "L-BFGS-B",
                  lower = ll, upper = ul,
                  control = list(fnscale = -1, maxit = maxit))
  if(oresG$convergence != 0){
    warning(sprintf("L-BFGS-B didn't converge: %s", oresG$message))
  }
  # compute fitted logit values L1 = log(mu[, 1] / mu[, 3]) and L2 = log(mu[, 2] / mu[, 3])
  oresG$fitted.logits <- dropAttributes(bclmLogits(oresG$par, x), c("betax", "params"))
  # compute fitted probabilities
  oresG$fitted.values <- dropAttributes(bclmProbabilities(oresG$par, x), c("betax", "params"))
  # add the info for predict function
  oresG$xlevels <- .getXlevels(mt, mf)
  oresG$terms <- mt
  if (model) 
    oresG$model <- mf
  # set the class of the result - necessary for logLik, coef and predict functions
  class(oresG) <- c("bclmRestricted", class(oresG))
  #names(ores$par) <- names(b)
  oresG
}

logLik.bclmRestricted <- function(object) {
  res <- object$value
  attr(res, "df") <- length(object$par)
  class(res) <- c("logLik", class(res))
  res
}

coef.bclmRestricted <- function(object) {
  object$par
}

predict.bclmRestricted <- function(object, newdata = NULL, type = c("link", "response", "terms")) {
  type <- match.arg(type)
  if(is.null(newdata)) {
    if(type == "link") {
      return(object$fitted.logits)
    }
    else if(type == "response") {
      return(object$fitted.values)
    }
    else stop(sprintf("Unsupported type: %s", type))
  }
  else {
    Terms <- delete.response(object$terms)
    m <- model.frame(Terms, newdata, xlev = object$xlevels)
    X <- model.matrix(Terms, m)
    if(attr(Terms, "intercept") == 1) { # we drop the intercept column - it will be implicitly used
      X <- X[, colnames(X) != "(Intercept)", drop = FALSE]
    }
    if(type == "link") {
      return(dropAttributes(bclmLogits(object$par, X), c("betax", "params")))
    }
    else if(type == "response") {
      return(dropAttributes(bclmProbabilities(object$par, X), c("betax", "params")))
    }
    else stop(sprintf("Unsupported type: %s", type))
  }
}

print.bclmRestricted <- function(object, digits = max(3, getOption("digits") - 3)) {
  cat("Call:\n")
  print(object$call)
  cat("\nCoefficients:\n")
  print(coef(object))
  cat("\nLog-likelihood: ")
  cat(round(logLik(object), digits = digits))
  cat(sprintf(" (df=%s)\n", attr(logLik(object), "df")))
  invisible(object)
}

############################################################################
# Adjacent Categories Logit Model
#
# Logits are given (for a 3-class problem) as
# \log(\frac{\pi_1}{\pi_2}) = \alpha_1 + \beta_1^Tx = L_1(x)
# \log(\frac{\pi_2}{\pi_3}) = \alpha_2 + \beta_2^Tx = L_2(x)
#   This model is a reparameterization of the Baseline Category Logit Model,
# where
# \log(\frac{\pi_1}{\pi_3}) = L_1(x) + L2(x)
# \log(\frac{\pi_2}{\pi_3}) = L_2(x)
#   The standard Adjacent Categories Logit Model assumes \beta_1 = \beta_2
# and this case corresponds to our prevailing theory and the unrestricted
# case without this constraint is equivalent to the standard multinomial
# logistic regression model (Baseline Category Logit Model), but has different
# parameterization. The restricted version (family = acat(parallel = TRUE) in
# vglm) is a special case of the restricted Baseline Category Logit Model
# with gamma = 0.5.
#   We use vglm function in VGAM package to estimate both restricted
# (prevailing theory) and unrestricted models.
#
#   Parameters:
# formula - a formula to be passed to the vglm function, describing the model
#           to estimate
# data - a data frame with the data to use for model estimation
# restricted - boolean, indicator if we want to fit a restricted model (TRUE),
#              the one that satisfies the prevailing theory assumption or
#              unrestricted one (FALSE).
#
#   Return value:
# An object returned by VGAM::vglm function.
aclm <- function(formula, data, restricted = TRUE) {
  if(restricted) {
    VGAM::vglm(formula = formula, data = data, family = VGAM::acat(parallel = TRUE))
  } else {
    VGAM::vglm(formula = formula, data = data, family = VGAM::acat(parallel = FALSE))
  }
}

############################################################################
# Proportional-Odds Cumulative Logit Model
#
# Logits are given (for a 3-class problem) as
# \log(\frac{\pi_1}{\pi_2 + \pi_3}) = \alpha_1 + \beta_1^Tx = L_1(x)
# \log(\frac{\pi_1 + \pi_2}{\pi_3}) = \alpha_2 + \beta_2^Tx = L_2(x)
#   This model is __NOT__ a reparameterization of the Baseline Category Logit
# Model, moreover, the unrestricted version of it, given by the above logits,
# may produce probabilities outside of the [0, 1] range.
#   The standard Proportional-Odds Cumulative Logit Model assumes
# \beta_1 = \beta_2 and this case corresponds to our prevailing theory. In
# addition, this restriction ensures that estimated probabilities (with proper
# estimation procedure) are in [0, 1] range.
#   We use vglm funciton in VGAM package to estimate both restricted
# (prevailing theory) and unrestricted models.
#
#   Parameters:
# formula - a formula to be passed to the vglm function, describing the model
#           to estimate
# data - a data frame with the data to use for model estimation
# restricted - boolean, indicator if we want to fit a restricted model (TRUE),
#              the one that satisfies the prevailing theory assumption or
#              unrestricted one (FALSE).
#
#   Return value:
# An object returned by VGAM::vglm function.
poclm <- function(formula, data, restricted = TRUE) {
  if(restricted) {
    VGAM::vglm(formula = formula, data = data, family = VGAM::cumulative(parallel = TRUE))
  } else {
    VGAM::vglm(formula = formula, data = data, family = VGAM::cumulative(parallel = FALSE))
  }
}

############################################################################
# Sequence of Binomial Logit Models. Special case for our 3-class problem.
#
# Since we assume a 3-class model with classes encoded as 1, 2 and 3:
# 1 - survival without sequelae,
# 2 - survival with sequelae,
# 3 - death,
# We first fit a logistic (binomial) regression model of death vs. survival,
# giving us estimates probabilities p_3 = \pi_3; and then fit a separate 
# (conditional) model of survival with sequelae vs survival without sequelae,
# giving us a probabilities p_2. These probabilities can be used to compute
# class probabilities:
# \pi_1 = (1 - p_3) * (1 - p_2)
# \pi_2 = (1 - p_3) * p_2
# \pi_3 = p_3.
# Model that does not make any assumptions on the coefficients of the two
# regressions corresponds to H_1 and the one that restricts them - to H_0.
# Logits for the unconstrained model are of the form
# L_1 = \alpha_1 + \beta_1^Tx # for p_2 estimation.
# L_2 = \alpha_2 + \beta_2^Tx # for p_3 estimation.
# In case of the constrained model:
# L_1 = \alpha_1 + \beta^Tx # for p_2 estimation.
# L_2 = \alpha_2 + \gamma \beta^Tx # for p_3 estimation.
# where \gamma > 0.
# We use our own implementation for both models. Unconstrained implementation
# fits two glm models and constrained model uses L-BFGS-B on the likelihood
# function to estimate the parameters.
#
#   Parameters:
# formula - a formula to be passed to the vglm function, describing the model
#           to estimate
# data - a data frame with the data to use for model estimation
# restricted - boolean, indicator if we want to fit a restricted model (TRUE),
#              the one that satisfies the prevailing theory assumption or
#              unrestricted one (FALSE).
# maxit - positive integer, maximal number of iterations for the L-BFGS-B
#         algorithm (optim function); ignored if restricted = FALSE.
# lambda - non-negative value of the regularization parameter, only used for
#          the restricted model; intercepts are not rgularized. Small values
#          of this parameter in [0.0001, 0.01] range help to deal with
#          degenerate cases when probabilities are estimated as 0/1.
#
#   Return value:
# An object returned either by soblmRestricted or soblmUnrestricted functions.
soblm <- function(formula, data, restricted = TRUE, maxit = 2000, lambda = 0) {
  #cl is not really needed, we only adds this call to the return object
  cl <- match.call()
  if(restricted) {
    res <- soblmRestricted(formula = formula, data = data, maxit = maxit, lambda = lambda)
  } else {
    res <- soblmUnrestricted(formula = formula, data = data)
  }
  res$call <- cl
  res
}

soblmUnrestricted <- function (formula, data, weights = NULL) {
  # prepare for the actual non-standard evaluation
  mf <- match.call(expand.dots = FALSE)
  #message("Call with ellipses not expanded: ")
  #note that there are no ellipses in the function arguments for now, 
  #but you might want to change that later
  #print(mf)
  #turn weights into symbol if character is passed
  if (is.character(mf$weights)) mf$weights <- as.symbol(mf$weights)
  m <- match(c("formula", "data", "weights", "maxit"), names(mf), 0L)
  #message("Position of formula, data, weights and maxit in the call:")
  #print(m)
  mf <- mf[c(1L, m)]
  #message("New call that only contains what is needed:")
  #print(mf)
  mf$drop.unused.levels <- TRUE 
  #message("Call with argument added:")
  #print(mf)
  mf[[1L]] <- quote(stats::model.frame) 
  #message("Change call to a call to model.frame:")
  #print(mf)
  mf <- eval(mf, parent.frame()) #evaluate call
  wt <- as.vector(model.weights(mf))
  stopifnot(is.null(wt)) # we do not support non-uniform weights
  y <- model.response(mf) # this is the response
  stopifnot(length(levels(y)) == 3) # only 3-level response is supported
  y <- as.numeric(y) # we convert response to numeric values 1, 2 and 3
  mt <- attr(mf, "terms")
  x <- model.matrix(mt, mf) # this is the design matrix, it may include an intercept
  # fit model for death vs. survival
  c3vsOther <- glm.fit(x, y == 3, family = binomial())
  # fit model for survival with sequelae vs. survival without sequelae
  c2vsc1 <- glm.fit(x[y != 3, , drop = FALSE], (y == 2)[y != 3], family = binomial())
  # compute estimated probabilities
  #p3 <- c3vsOther$fitted.values
  L <- cbind(as.vector(x %*% c2vsc1$coefficients),
             as.vector(x %*% c3vsOther$coefficients))
  colnames(L) <- c("log(mu[, 2] / mu[, 1])", "log(mu[, 3] / (1 - mu[, 3]))")
  p3 <- sigmoid(L[, 2])
  p2 <- sigmoid(L[, 1])
  fitted.values <- cbind((1 - p2) * (1 - p3),
                         p2 * (1 - p3),
                         p3)
  colnames(fitted.values) <- 1:3
  # compute log-likelihood
  ll <- sum((y == 1) * log(fitted.values[, 1]) + 
            (y == 2) * log(fitted.values[, 2]) +
            (y == 3) * log(fitted.values[, 3]))
  attr(ll, "df") <- 2 * ncol(x)
  class(ll) <- c("logLik", class(ll))
  # prepare the value to return
  res <- list(c3vsOther = c3vsOther, c2vsc1 = c2vsc1, L = L, fitted.values = fitted.values, logLik = ll)
  # set the class of the result - necessary for logLik, coef and predict functions
  class(res) <- c("soblmUnrestricted", class(res))
  res
}

logLik.soblmUnrestricted <- function(object) {
  object$logLik
}

coef.soblmUnrestricted <- function(object) {
  v2 <- object$c3vsOther$coefficients
  v1 <- object$c2vsc1$coefficients
  names(v1) <- paste(names(v1), 1, sep = ":")
  names(v2) <- paste(names(v2), 2, sep = ":")
  structure(as.vector(rbind(v1, v2)),
            names = as.vector(rbind(names(v1), names(v2))))
}

predict.soblmUnrestricted <- function(object, newdata = NULL, type = c("link", "response", "terms")) {
  type <- match.arg(type)
  if(is.null(newdata)) {
    if(type == "link") {
      return(object$L)
    }
    else if(type == "response") {
      return(object$fitted.values)
    }
    else stop(sprintf("Unsupported type: %s", type))
  }
  else {
    if(type == "link") {
      L1 <- predict(object$c2vsc1, newdata, type = "link")
      L2 <- predict(object$c3vsOther, newdata, type = "link")
      res <- cbind(L1, L2)
      colnames(res) <- c("log(mu[, 2] / mu[, 1])", "log(mu[, 3] / (1 - mu[, 3]))")
      return(res)
    }
    else if(type == "response") {
      p2 <- predict(object$c2vsc1, newdata, type = "response")
      p3 <- predict(object$c3vsOther, newdata, type = "response")
      res <- cbind((1 - p2) * (1 - p3),
                   p2 * (1 - p3),
                   p3)
      colnames(res) <- 1:3
      return(res)
    }
    else stop(sprintf("Unsupported type: %s", type))
  }
}

print.soblmUnrestricted <- function(object, digits = max(3, getOption("digits") - 3)) {
  cat("Call:\n")
  print(object$call)
  cat("\nCoefficients:\n")
  print(coef(object))
  cat("\nLog-likelihood: ")
  cat(round(logLik(object), digits = digits))
  cat(sprintf(" (df=%s)\n", attr(logLik(object), "df")))
  invisible(object)
}

soblmRestrictedLogits <- function(b, x) {
  alpha1 <- b[1]
  alpha2 <- b[2]
  gamma <- b[3]
  beta <- b[-(1:3)]
  betax <- as.vector(x %*% beta)
  L <- cbind(alpha1 + betax, alpha2 + gamma * betax)
  colnames(L) <- c("log(mu[, 2] / mu[, 1])", "log(mu[, 3] / (1-mu[, 3]))")
  list(alpha1 = alpha1,
       alpha2 = alpha2,
       gamma = gamma,
       beta = beta,
       betax = betax,
       L = L)
}

soblmRestrictedProbabilities <- function(b, x) {
  logits <- soblmRestrictedLogits(b, x)
  logits$p2 <- sigmoid(logits$L[, 1])
  logits$p3 <- sigmoid(logits$L[, 2])
  logits
}

soblmRestrictedFinalProbabilities <- function(b, x) {
  probs <- soblmRestrictedProbabilities(b, x)
  pis <- cbind((1 - probs$p2) * (1 - probs$p3),
               probs$p2 * (1 - probs$p3),
               probs$p3)
  colnames(pis) <- 1:3
  probs$pis <- pis
  probs
}

# Compute the log-likelihood function for restricted Sequence-of-Binomials model.
#
#   Parameters:
# b - vector of parameters, first two components correspond to alpha_1 and _alpha2,
#     third - to gamma and the rest to beta.
# x - n x p design matrix.
# y - responds vector of integers in {1, 2, 3} set.
# lambda - regularization parameter, small value of it (e.g. 0.0001~0.01) help to deal with
#          degenerate cases when some probabilities are estimated as 0/1.
#
#   Return value:
# Value of the (penalized) log-likelihood function. To get un-penalized
# log-likelihood, set lambda = 0.
soblmRestrictedF <- function(b, x, y, lambda = 0) {
  logits <- soblmRestrictedLogits(b, x)
  v1 <- log1p(exp(logits$L[, 1]))
  v2 <- log1p(exp(logits$L[, 2]))
  idx <- v1 == Inf
  v1[idx] <- logits$L[idx, 1]
  idx <- v2 == Inf
  v2[idx] <- logits$L[idx, 2]
  res <- sum(-v2 - (y !=3) * v1 + (y == 2) * logits$L[, 1] + I(y == 3) * logits$L[, 2])
           - lambda * (sum(logits$beta^2) + logits$gamma^2)
  res
}

# Compute the gradient of the log-likelihood function for restricted Sequence-of-Binomials model.
#
#   Parameters:
# b - vector of parameters, first two components correspond to alpha_1 and _alpha2,
#     third - to gamma and the rest to beta.
# x - n x p design matrix.
# y - responds vector of integers in {1, 2, 3} set.
# lambda - regularization parameter, small value of it (e.g. 0.0001~0.01) help to deal with
#          degenerate cases when some probabilities are estimated as 0/1.
#
#   Return value:
# Value of the (penalized) log-likelihood function.
soblmRestrictedG <- function(b, x, y, lambda = 0) {
  probs <- soblmRestrictedProbabilities(b, x)
  res <- c(sum((y == 2) - (y != 3) * probs$p2), # df/dalpha1
           sum((y == 3) - probs$p3), # df/dalpha2
           sum(((y == 3) - probs$p3) * probs$betax) - 2 * lambda * probs$gamma, # df/dgamma
           colSums(x * (probs$gamma * ((y == 3) - probs$p3) + (y == 2) - (y != 3) * probs$p2)) 
             - 2 * lambda * probs$beta # df/dbeta
          )
  names(res) <- names(b)
  res
}

soblmRestricted <- function (formula, data, weights = NULL, maxit = 2000, lambda = 0, model = TRUE) {
  #cl is not really needed, we only adds this call to the return object
  cl <- match.call()
  # prepare for non-standard evaluation of the call
  mf <- match.call(expand.dots = FALSE)
  #message("Call with ellipses not expanded: ")
  #note that there are no ellipses in the function arguments for now, 
  #but you might want to change that later
  #print(mf)
  #turn weights into symbol if character is passed
  if (is.character(mf$weights)) mf$weights <- as.symbol(mf$weights)
  m <- match(c("formula", "data", "weights"), names(mf), 0L)
  #message("Position of formula, data, weights and maxit in the call:")
  #print(m)
  mf <- mf[c(1L, m)]
  #message("New call that only contains what is needed:")
  #print(mf)
  mf$drop.unused.levels <- TRUE 
  #message("Call with argument added:")
  #print(mf)
  mf[[1L]] <- quote(stats::model.frame) 
  #message("Change call to a call to model.frame:")
  #print(mf)
  mf <- eval(mf, parent.frame()) #evaluate call
  wt <- as.vector(model.weights(mf))
  stopifnot(is.null(wt)) # we do not support non-uniform weights
  y <- model.response(mf) # this is the response
  stopifnot(length(levels(y)) == 3) # only 3-level response is supported
  y <- as.numeric(y) # we convert response to numeric values
  mt <- attr(mf, "terms")
  x <- model.matrix(mt, mf) # this is the design matrix, it may include an intercept
  stopifnot(attr(mt, "intercept") == 1)
  if(attr(mt, "intercept") == 1) { # we drop the intercept column - it will be implicitly used
    x <- x[, colnames(x) != "(Intercept)", drop = FALSE]
  }
  b <- c(0, 0, 1, rep(0, ncol(x))) # initialize the parameter vector
  names(b) <- c("(Intercept):1", "(Intercept):2", "gamma", colnames(x))
  ll <- rep(-Inf, length(b))
  ul <- rep(Inf, length(b))
  ll[3] <- 0
  oresG <- optim(par = b, fn = soblmRestrictedF, gr = soblmRestrictedG,
                 x = x, y = y, lambda = lambda,
                 method = "L-BFGS-B",
                 lower = ll, upper = ul,
                 control = list(fnscale = -1, maxit = maxit))
  if(oresG$convergence != 0){
    warning(sprintf("L-BFGS-B didn't converge: %s", oresG$message))
  }
  # compute fitted probabilities
  probs <- soblmRestrictedFinalProbabilities(oresG$par, x)
  oresG$fitted.values <- probs$pis
  oresG$fitted.logits <- probs$L
  oresG$lambda = lambda
  # add the info for predict function
  oresG$xlevels <- .getXlevels(mt, mf)
  oresG$terms <- mt
  if (model) 
    oresG$model <- mf
  # set the class of the result - necessary for logLik, coef and predict functions
  class(oresG) <- c("soblmRestricted", class(oresG))
  #names(ores$par) <- names(b)
  oresG
}

logLik.soblmRestricted <- function(object) {
  res <- object$value + object$lambda * sum(object$par[-(1:2)]^2)
  attr(res, "df") <- length(object$par)
  class(res) <- c("logLik", class(res))
  res
}

coef.soblmRestricted <- function(object) {
  object$par
}

predict.soblmRestricted <- function(object, newdata = NULL, type = c("link", "response", "terms")) {
  type <- match.arg(type)
  if(is.null(newdata)) {
    if(type == "link") {
      return(object$fitted.logits)
    }
    else if(type == "response") {
      return(object$fitted.values)
    }
    else stop(sprintf("Unsupported type: %s", type))
  }
  else {
    Terms <- delete.response(object$terms)
    m <- model.frame(Terms, newdata, xlev = object$xlevels)
    X <- model.matrix(Terms, m)
    if(attr(Terms, "intercept") == 1) { # we drop the intercept column - it will be implicitly used
      X <- X[, colnames(X) != "(Intercept)", drop = FALSE]
    }
    if(type == "link") {
      return(soblmRestrictedLogits(object$par, X)$L)
    }
    else if(type == "response") {
      return(soblmRestrictedFinalProbabilities(object$par, X)$pis)
    }
    else stop(sprintf("Unsupported type: %s", type))
  }
}

print.soblmRestricted <- function(object, digits = max(3, getOption("digits") - 3)) {
  cat("Call:\n")
  print(object$call)
  cat("\nCoefficients:\n")
  print(coef(object))
  cat("\nLog-likelihood: ")
  cat(round(logLik(object), digits = digits))
  cat(sprintf(" (df=%s)\n", attr(logLik(object), "df")))
  cat("lambda:", object$lambda, "\n")
  invisible(object)
}

#############################################################################
# Perform a likelihood ratio test of H_0: prevailing theory holds vs.
# H_1: general form of multinomial regression is appropriate. Use both 
# parametric bootstrap and Wilks' theorem.
#
#   Parameters:
# formula - a formula for a model to fit.
# data - a data frame with the data.
# h0ModelFunc - a function to fit the restricted model (the one that
#               satisfies the prevailing theorem, e.g. bclm, aclm or
#               poclm with restricted = TRUE)
# h1ModelFunc - a function to fit the unrestricted model (the one that
#               satisfies the alternative hypothesis, e.g. bclm, aclm or
#               poclm with restricted = FALSE)
# B - number of bootstrap samples to use for p-value estimation.
#
#   Return value:
# A list containing 3 components:
LRTest <- function(formula, data, h0ModelFunc, h1ModelFunc, B = 10) {
  # Helper function to compute Likelihood ratio statistic.
  #
  #   Parameters:
  # h0 - model fitted for the null hypothesis.
  # h1 - model fitted for the alternative hypothesis.
  #
  #   Return value:
  # An object of "LogLikStatistic" class, which is a single numer with an attribute
  # "df" - number of degrees of freedom to be used in the Wilks theorem.
  LRTestStatistic <- function(h0, h1) {
    # Test statistic value is equal to
    # -2 ln \frac{\sup_{\theta \in \Theta_0}L(\theta)}{\sup_{\theta \in \Theta_1}L(\theta)} = 
    # -2(ln(\sup_{\theta \in \Theta_0}L(\theta)) - ln(\sup_{\theta \in \Theta_1}L(\theta)))
    # and it has (asymptotically) \chi^2 distribution with p degrees of freedom,
    # where p = (num_parameters in H1) - (num_parameters in H0).
    res <- -2 * (logLik(h0) - logLik(h1))
    attr(res, "df") <- length(coef(h1)) - length(coef(h0))
    class(res) <- c("LogLikStatistic", class(res))
    res
  }
  # estimate the constrained model
  h0 <- h0ModelFunc(formula, data)
  h1 <- h1ModelFunc(formula, data)
  # compute predicted (estimated) probabilities - we shall use them for 
  # for bootstrap data generation
  p <- predict(h0, type = "response")
  statistic <- LRTestStatistic(h0, h1)
  # generate bootstrap statistics
  bootstrapStatistics <- sapply(1:B,
                                function(i) {
                                  y <- apply(p, 1, function(x) sample(1:3, 1, prob = x))
                                  if(any(is.na(y))) {
                                    print(y)
                                    stop("NAs in generated response y.")
                                  }
                                  data$disoutcome <- ordered(y)
                                  h0 <- h0ModelFunc(formula, data)
                                  h1 <- h1ModelFunc(formula, data)
                                  LRTestStatistic(h0, h1)
                                })
  # return the result
  res <- list(statistic = statistic,
              bootstrapStatistics = bootstrapStatistics,
              pvalue = list(Wilks = pchisq(q = statistic, df = attr(statistic, "df"), lower.tail = FALSE),
                            bootstrap = mean(bootstrapStatistics >= statistic)),
              h0 = h0,
              h1 = h1)
  class(res) <- c("LRTest", class(res))
  res
}

print.LRTest <- function(o) {
  pvalues <- c(o$pvalue$Wilks, o$pvalue$bootstrap, attr(o$statistic, "df"), length(o$bootstrapStatistics))
  # TODO: improve output formatting
  names(pvalues) <- c("Wilks", "bootstrap", "df", "B")
  print(pvalues)
  invisible(pvalues)
}

##########################################################################
# An example of running the experiment.
##########################################################################
#
#   Parameters:
# dataFile - string, name of the file where data are stored, this parameter is passed to fLoadData function.
# B - positive integer, number of bootstrap samples to use when estimating the p-value.
#
#   Return value:
# An object returned by LRTest function.
example <- function(dataFile = "data/new_data_1.csv", B = 2000){ # use 2000 bootstrap samples to estimate p-value.
  # load data from disk
  d <- fLoadData(dataFile)
  # run the bootstrap likelihood ratio test; it also produces p-value according
  # to the Wilks theorem; this call will take about 200 seconds.
  res <- LRTest(disoutcome ~ age + comasc + log(comadur), # formula describing the model 
                data = d, # this is the 
                # function for fitting a model under H_0, we use bclm function, but any other could be used
                h0ModelFunc = function(formula, data) soblm(formula, data, restricted = TRUE, lambda = 0.0001),
                # function for fitting a model under H_1, we use bclm function, but any other could be used
                h1ModelFunc = function(formula, data) soblm(formula, data, restricted = FALSE),
                B = B) # use 2000 bootstrap samples to estimate p-value.
  res
}