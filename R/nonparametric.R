##' Estimates non-parametric MIDAS regression
##'
##' Estimates non-parametric MIDAS regression accodring Breitung et al.
##'
##' @title Estimate non-parametric MIDAS regression
##' @param formula formula specifying MIDAS regression
##' @param data a named list containing data with mixed frequencies
##' @param lambda smoothing parameter, defaults to \code{NULL}, which means that it is chosen by minimising AIC.
##' @return a \code{midas_r_np} object
##' @author Vaidotas Zemlys
##' @references Breitung J, Roling C, Elengikal S (2013). \emph{Forecasting inflation rates using daily data: A nonparametric MIDAS approach} Working paper, URL http://www.ect.uni-bonn.de/mitarbeiter/joerg-breitung/npmidas.
##' @export
##' @import Matrix
##' @importFrom stats optimize
##' @examples
##' data("USunempr")
##' data("USrealgdp")
##' y <- diff(log(USrealgdp))
##' x <- window(diff(USunempr),start=1949)
##' trend <- 1:length(y)
##' midas_r_np(y~trend+fmls(x,12,12))
midas_r_np <- function(formula, data, lambda = NULL) {
  Zenv <- new.env(parent = environment(formula))

  if (missing(data)) {
    ee <- NULL
  }
  else {
    ee <- data_to_env(data)
    parent.env(ee) <- parent.frame()
  }

  assign("ee", ee, Zenv)
  formula <- as.formula(formula)
  cl <- match.call()
  mf <- match.call(expand.dots = FALSE)
  mf$formula <- formula
  m <- match(c("formula", "data"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf[[1L]] <- as.name("model.frame")
  mf[[3L]] <- as.name("ee")
  mf[[4L]] <- as.name("na.omit")
  names(mf)[c(2, 3, 4)] <- c("formula", "data", "na.action")

  mf <- eval(mf, Zenv)
  mt <- attr(mf, "terms")

  terms.lhs <- as.list(attr(mt, "variables"))[-2:-1]
  term.labels <- attr(mt, "term.labels")

  rfd <- vector("list", length(terms.lhs))
  yname <- all.vars(mt[[2]])

  for (i in 1:length(rfd)) {
    fr <- terms.lhs[[i]]
    # This a behaviour of R we rely on. It might be non-standard one.
    fun <- as.character(fr)[1]
    rfd[[i]] <- if (fun %in% c("fmls", "mls", "dmls")) {
      lags <- eval(fr[[3]], Zenv)
      nm <- as.character(fr[[2]])
      nol <- switch(fun,
        fmls = lags + 1,
        dmls = lags + 1,
        mls = length(lags)
      )
      lagstruct <- switch(fun,
        fmls = 0:lags,
        dmls = 0:lags,
        mls = lags
      )

      if (nol < 3 & nm != yname) stop("For nonparametric MIDAS you need at least 3 high frequency lags")

      wlab <- ifelse(nm == yname, "", nm)
      list(
        weight = function(p) p,
        term_name = nm,
        gradient = NULL,
        start = NULL,
        weight_name = "non-parametric weight",
        frequency = eval(fr[[4]], Zenv),
        lag_structure = lagstruct
      )
    } else {
      list(
        weight = function(p) p,
        term_name = term.labels[i],
        gradient = NULL,
        start = NULL,
        weight_name = "",
        frequency = 1,
        lag_structure = 0
      )
    }
  }

  if (attr(mt, "intercept") == 1) {
    rfd <- c(list(list(
      weight = function(p) p,
      term_name = "(Intercept)",
      gradient = NULL,
      start = NULL,
      weight_name = "",
      frequency = 1,
      lag_structure = 0
    )), rfd)
  }


  names(rfd) <- sapply(rfd, "[[", "term_name")

  ## Note this is a bit of misnomer. Variable weight_names is actualy a vector of term names which have MIDAS weights.
  ## It *is not* the same as actual name of weight function. This is a left-over from the old code. Grabbed from prepmidas_r

  weight_names <- sapply(rfd, "[[", "weight_name")
  weight_inds <- which(weight_names != "")
  weight_names <- names(rfd)[weight_names != ""]

  lengths <- sapply(lapply(rfd, "[[", "lag_structure"), length)

  build_indices <- function(ci, nm) {
    inds <- cbind(c(1, ci[-length(ci)] + 1), ci)
    inds <- apply(inds, 1, function(x) list(x[1]:x[2]))
    inds <- lapply(inds, function(x) x[[1]])
    names(inds) <- nm
    inds
  }

  pinds <- build_indices(cumsum(lengths), names(rfd))

  if (length(weight_names) > 1) stop("Only one non-autoregressive mixed frequency term is currently supported")

  resplace <- pinds[[weight_names]][1]
  rno <- length(rfd[[weight_names]]$lag_structure)

  y <- model.response(mf, "numeric")
  X <- model.matrix(mt, mf)

  k <- ncol(X)
  if (k < nrow(X)) {
    unrestricted <- lm(y ~ . - 1, data = data.frame(cbind(y, X), check.names = FALSE))
  } else {
    unrestricted <- NULL
  }

  D <- bandSparse(rno - 2, k, resplace - 1 + c(0, 1, 2), diagonals = list(rep(1, rno - 2), rep(-2, rno - 2), rep(1, rno - 2)))

  DD <- crossprod(D)
  ol <- opt_lambda(y, X, DD, lambda)

  fit <- X %*% ol$beta
  res <- y - fit

  cf <- as.numeric(ol$beta)
  names(cf) <- names(unlist(pinds))

  term_info <- rfd
  names(term_info) <- sapply(term_info, "[[", "term_name")
  term_info <- mapply(function(term, pind, xind) {
    term$start <- NULL
    term$coef_index <- pind
    term$midas_coef_index <- pind
    term
  }, term_info, pinds[names(term_info)], SIMPLIFY = FALSE)


  out <- list(
    coefficients = cf,
    midas_coefficients = cf,
    model = cbind(y, X),
    unrestricted = unrestricted,
    call = cl,
    terms = mt,
    fitted.values = as.numeric(fit),
    residuals = as.numeric(res),
    term_info = term_info,
    lambda = ol$lambda,
    klambda = ol$klambda,
    AIC = ol$AIC,
    opt = ol$opt,
    Zenv = Zenv
  )
  class(out) <- "midas_r_np"
  out
}

##' @importFrom stats AIC
##' @export
##' @method AIC midas_r_np
AIC.midas_r_np <- function(object, ..., k) {
  object$AIC(object$lambda)
}

##' @importFrom stats BIC
##' @export
##' @method BIC midas_r_np
BIC.midas_r_np <- function(object, ..., k) {
  # for compatibility
  object$AIC(object$lambda)
}

##' @export
##' @method forecast midas_r_np
forecast.midas_r_np <- forecast.midas_r

##' @export
##' @method predict midas_r_np
predict.midas_r_np <- predict.midas_r

##' @export
##' @method print midas_r_np
print.midas_r_np <- function(x, ...) {
  cat("Nonparametric MIDAS regression model", paste0("(", nrow(x$model)), "low frequency observations)")
  cat("\nFormula: ", deparse(formula(terms(x))))
  cat("\nThe smoothing parameter: ", x$lambda)
  cat("\nThe effective number of parameters:", x$klambda)
  cat("\nAIC of the model: ", AIC(x))
  cat("\nRoot mean squared error: ", sqrt(mean(residuals(x)^2)), "\n")
}

##' @export
##' @method summary midas_r_np
summary.midas_r_np <- function(object, ...) {
  print(object, ...)
}

opt_lambda <- function(y, X, DD, lambda) {
  n <- length(y)
  XX <- crossprod(X)
  Xy <- crossprod(X, y)
  tX <- t(X)
  AIC <- function(lambda) {
    Qlambda <- XX + lambda * n * DD
    klambda <- sum(X * t(solve(Qlambda, tX)))
    beta <- solve(Qlambda, Xy)
    res <- y - X %*% beta
    log(sum(res^2)) + 2 * (klambda + 1) / (n - klambda - 2)
  }
  if (is.null(lambda)) {
    opt <- optimize(AIC, c(0, 100))
    lambda <- opt$minimum
  }
  else {
    opt <- NULL
  }
  Qlambda <- XX + lambda * n * DD
  klambda <- sum(X * t(solve(Qlambda, tX)))
  beta <- solve(Qlambda, Xy)
  AICl <- AIC(lambda)
  list(beta = beta, klambda = klambda, lambda = lambda, AIC = AIC, opt = opt)
}
