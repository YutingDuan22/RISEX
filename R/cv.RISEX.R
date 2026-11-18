
#' K-fold cross-validation for tuning \code{zeta} and \code{zeta0} in RISE-X
#'
#' Selects the spatial smoothing penalty \code{zeta} and the sparsity penalty
#' \code{zeta0} of \code{\link{RISEX}} by K-fold cross-validation
#' over subjects. For each fold, the model is fitted on the training subjects
#' for every \code{(zeta, zeta0)} pair on the grid, and the fitted model is
#' scored on the held-out subjects with a robust (truncated) weighted squared
#' error whose voxel weights are fixed within each fold.
#'
#' @param img Numeric matrix (\eqn{V \times N}) of images, one column per subject.
#' @param voxels Numeric matrix (\eqn{V \times 3}) of voxel coordinates.
#' @param cov_mat Numeric matrix (\eqn{N \times p}) of covariates.
#' @param psi_mat Numeric matrix (\eqn{V \times K}) of eigenfunctions, e.g. from
#'   \code{\link{voxel.wise.fit}}. Used for every fold unless
#'   \code{refit_basis = TRUE}.
#' @param lambda Eigenvalues corresponding to \code{psi_mat}.
#' @param zeta_grid Numeric vector of candidate \code{zeta} values. Default is
#'   \code{c(0.5, 1, 2, 4, 8) / sigma0^2}, which brackets the
#'   \code{RISEX} default \code{2 / sigma0^2}.
#' @param zeta0_grid Numeric vector of candidate \code{zeta0} values. Default is
#'   \code{c(0.01, 0.05, 0.1, 0.5, 1) / sigma0^2}, which brackets the
#'   \code{RISEX} default \code{0.1 / sigma0^2}.
#' @param tau Smoothness threshold passed to \code{RISEX}. Either a
#'   numeric scalar / length-\eqn{p} vector used for every fit, or a function
#'   \code{function(beta0, adj_matrix)} evaluated on each fold's initial
#'   coefficients (e.g. \code{calculate_tau}).
#' @param nfolds Number of folds. Default \code{5}.
#' @param foldid Optional integer vector of length \eqn{N} giving fold
#'   membership of each subject. Overrides \code{nfolds}.
#' @param adj Adjacency list from \code{\link{coordinates_to_adjacency_matrix.nf}}.
#'   Computed once from \code{voxels} if not supplied.
#' @param sigma0 Global noise scale used to build the default grids. If
#'   \code{NULL}, estimated as the SD of residuals from \code{voxel.wise.fit}-style
#'   OLS on the full data (only used for scaling the grids).
#' @param zeta_gamma Outlier threshold passed to \code{RISEX} and also
#'   used as the truncation level of the CV loss. Default \code{2.5^2}.
#' @param tau0 Sparsity scale passed to \code{RISEX}. Default
#'   \code{(1e-3)^2}.
#' @param loss Held-out loss. \code{"robust"} (default) and \code{"mse"} are
#'   described below. \code{"wmse"} is WMSE exactly as written in the paper:
#'   for each held-out subject, \eqn{\sum_s \{Y_i(s) - \hat Y_i(s)\}^2 /
#'   \hat\sigma^2(s)}, averaged over held-out subjects, where \eqn{\hat\beta},
#'   \eqn{\hat\sigma(\cdot)} come from the training fit and
#'   \eqn{\hat\eta_i = \Psi\hat\alpha_i} uses the \eqn{\alpha} update.
#'   Because the weights \eqn{\hat\sigma^{2(-k)}} change with the tuning values
#'   (heavier penalties inflate \eqn{\hat\sigma}), it can favour all-zero fits;
#'   it is kept for comparison. \code{"robust"} uses the truncated
#'   standardized squared error \eqn{\min(r^2/\sigma^2, \zeta_\gamma)}, matching
#'   the RISE-X data-fit term so extreme signals in test subjects do not dominate
#'   the CV error; \code{"mse"} uses the plain standardized squared error.
#' @param refit_basis Logical; if \code{TRUE}, re-estimate \code{psi_mat} and
#'   \code{lambda} on each training fold with \code{\link{voxel.wise.fit}}
#'   (avoids using held-out subjects in the basis, but is slow). Default
#'   \code{FALSE}.
#' @param refit Logical; if \code{TRUE} (default), refit \code{RISEX}
#'   on all subjects at the selected tuning values and return it.
#' @param rule Which selected pair to use for the refit: \code{"min"} (lowest
#'   CV error) or \code{"1se"} (the most regularized pair within one standard
#'   error of the minimum; ties broken by larger \code{zeta0}, then larger
#'   \code{zeta}).
#' @param n_cores Number of cores for \code{parallel::mclapply} over
#'   fold-by-grid jobs. Forced to \code{1} on Windows.
#' @param seed Optional seed for fold assignment.
#' @param verbose Logical; print progress. Default \code{TRUE}.
#' @param ... Further arguments passed to \code{RISEX}
#'   (e.g. \code{max_it}, \code{tol.beta}, \code{update_sigma}).
#'
#' @return An object of class \code{"cv.risex"}, a list containing:
#' \describe{
#'   \item{cv_table}{Data frame with one row per \code{(zeta, zeta0)} pair:
#'     mean CV error \code{cvm}, its standard error \code{cvsd} across folds,
#'     and the average number of nonzero coefficients \code{nonzero}.}
#'   \item{cv_folds}{Matrix of per-fold CV errors (grid pairs \eqn{\times} folds).}
#'   \item{zeta.min, zeta0.min}{Pair minimizing the CV error.}
#'   \item{zeta.1se, zeta0.1se}{Pair chosen by the one-standard-error rule.}
#'   \item{foldid}{Fold assignment used.}
#'   \item{fit}{Full-data \code{RISEX} fit at the chosen pair
#'     (if \code{refit = TRUE}).}
#' }
#'
#' @details
#' Folds are formed over subjects (columns of \code{img}). For a held-out
#' subject \eqn{i} with residual \eqn{r_i = y_i - \hat\beta x_i}, the
#' subject-specific latent scores are predicted by the same penalized
#' (BLUP-type) update used inside \code{RISEX},
#' \eqn{\hat\alpha_i = (\Psi^\top \Sigma^{-1} \Psi + \Lambda^{-1})^{-1}
#' \Psi^\top \Sigma^{-1} (r_i - \hat\gamma_i)}, iterated with hard-threshold
#' outlier trimming of \eqn{\hat\gamma_i}. The loss is computed on
#' \eqn{r_i - \Psi \hat\alpha_i}.
#'
#' For \code{loss = "robust"} or \code{"mse"}, to make errors comparable
#' across tuning values, residuals are standardized
#' by a fixed per-fold reference scale (robust MAD of training-fold OLS
#' residuals), not by each fit's own \code{sigma} estimate.
#'
#' Prediction-oriented CV tends to favour less sparse fits than are optimal for
#' activation-region recovery; \code{rule = "1se"} is a common remedy.
#'
#' @examples
#' \dontrun{
#' init <- voxel.wise.fit(img, cov_mat, voxels)
#' adj  <- coordinates_to_adjacency_matrix.nf(voxels)
#' cv   <- cv.RISEX(img, voxels, cov_mat,
#'                           psi_mat = init$psi_mat, lambda = init$lambda,
#'                           sigma0 = init$sigma0, tau = calculate_tau,
#'                           adj = adj, nfolds = 5, n_cores = 4)
#' cv$cv_table
#' plot(cv)
#' beta_hat <- cv$fit$beta
#' }
#'
#' @importFrom stats sd median mad
#' @importFrom parallel mclapply
#' @export
cv.RISEX <- function(img,
                     voxels,
                     cov_mat,
                     psi_mat,
                     lambda,
                     zeta_grid = NULL,
                     zeta0_grid = NULL,
                     tau,
                     nfolds = 5,
                     foldid = NULL,
                     adj = NULL,
                     sigma0 = NULL,
                     zeta_gamma = 2.5^2,
                     tau0 = (1e-3)^2,
                     loss = c("robust", "wmse", "mse"),
                     refit_basis = FALSE,
                     refit = TRUE,
                     rule = c("min", "1se"),
                     n_cores = 1,
                     seed = NULL,
                     verbose = TRUE,
                     ...) {
  loss <- match.arg(loss)
  rule <- match.arg(rule)
  cov_mat <- as.matrix(cov_mat)
  img <- as.matrix(img)
  N <- ncol(img)
  V <- nrow(img)
  if (nrow(cov_mat) != N) stop("nrow(cov_mat) must equal ncol(img).")
  ## ---- input checks (common causes of NA CV losses) ----
  if (anyNA(img)) stop("`img` contains NA/NaN values (e.g. out-of-mask voxels). Remove those voxels (rows) first.")
  if (anyNA(cov_mat)) stop("`cov_mat` contains NA values. Remove or impute those subjects first.")
  if (anyNA(psi_mat) || anyNA(lambda)) stop("`psi_mat` or `lambda` contains NA values.")
  if (nrow(voxels) != V) stop("nrow(voxels) must equal nrow(img).")
  if (nrow(psi_mat) != V) stop("nrow(psi_mat) must equal nrow(img).")
  if (missing(tau)) stop("Please supply `tau` (a number/vector or a function(beta0, adj_matrix)).")
  lambda <- as.numeric(lambda)
 
  if (is.null(adj)) adj <- coordinates_to_adjacency_matrix.nf(voxels)
  if (sum(adj$adj_x) + sum(adj$adj_y) + sum(adj$adj_z) == 0) {
    stop("No neighbouring voxels found. `voxels` must be integer grid indices ",
         "(neighbours differ by 1), not mm/MNI coordinates.")
  }
 
  ## ---- default grids (relative to sigma0) ----
  if (is.null(zeta_grid) || is.null(zeta0_grid)) {
    if (is.null(sigma0)) {
      sigma0 <- sd(c(.ols_fit(img, cov_mat)$resid))
      if (verbose) message(sprintf("sigma0 not supplied; using OLS residual SD = %.4g for default grids.", sigma0))
    }
    if (is.null(zeta_grid))  zeta_grid  <- c(0.5, 1, 2, 4, 8) / sigma0^2
    if (is.null(zeta0_grid)) zeta0_grid <- c(0.01, 0.05, 0.1, 0.5, 1) / sigma0^2
  }
  grid <- expand.grid(zeta = sort(zeta_grid), zeta0 = sort(zeta0_grid))
  n_grid <- nrow(grid)
 
  ## ---- folds over subjects ----
  if (is.null(foldid)) {
    if (!is.null(seed)) set.seed(seed)
    if (nfolds < 2 || nfolds > N) stop("nfolds must be between 2 and N.")
    foldid <- sample(rep(seq_len(nfolds), length.out = N))
  } else {
    if (length(foldid) != N) stop("foldid must have length ncol(img).")
    foldid <- as.integer(factor(foldid))
  }
  K <- max(foldid)
 
  ## ---- per-fold quantities shared across the grid ----
  fold_info <- lapply(seq_len(K), function(k) {
    tr <- which(foldid != k)
    te <- which(foldid == k)
    ols <- .ols_fit(img[, tr, drop = FALSE], cov_mat[tr, , drop = FALSE])
    # reference scale: robust per-voxel MAD of OLS residuals, shrunk to the mean
    s_ref <- apply(ols$resid, 1, mad)
    s_ref[!is.finite(s_ref) | s_ref <= 0] <- mean(s_ref[is.finite(s_ref) & s_ref > 0])
    s_ref <- 0.5 * s_ref + 0.5 * mean(s_ref)
    if (refit_basis) {
      init <- voxel.wise.fit(img[, tr, drop = FALSE], cov_mat[tr, , drop = FALSE], voxels)
      psi_k <- init$psi_mat
      lam_k <- as.numeric(init$lambda)
      sig_init <- init$sigma1
    } else {
      psi_k <- psi_mat
      lam_k <- lambda
      sig_init <- NULL
    }
    tau_k <- if (is.function(tau)) tau(ols$beta, adj$adj_matrix) else tau
    list(tr = tr, te = te, beta0 = ols$beta, s_ref = s_ref,
         psi = psi_k, lambda = lam_k, sigma = sig_init, tau = tau_k)
  })
 
  ## ---- fit + score every (fold, grid) job ----
  jobs <- expand.grid(g = seq_len(n_grid), k = seq_len(K))
  run_job <- function(j) {
    g <- jobs$g[j]; k <- jobs$k[j]
    fi <- fold_info[[k]]
    fit <- tryCatch(
      RISEX(
        img = img[, fi$tr, drop = FALSE],
        voxels = voxels,
        cov_mat = cov_mat[fi$tr, , drop = FALSE],
        psi_mat = fi$psi,
        lambda = fi$lambda,
        adj = adj,
        beta0 = fi$beta0,
        sigma = fi$sigma,
        zeta = grid$zeta[g],
        zeta0 = grid$zeta0[g],
        tau = fi$tau,
        tau0 = tau0,
        zeta_gamma = zeta_gamma,
        trace = FALSE,
        plot = FALSE,
        ...
      ),
      error = function(e) e
    )
    if (inherits(fit, "error")) {
      return(list(err = NA_real_, nonzero = NA_real_, msg = conditionMessage(fit)))
    }
    err <- .cv_heldout_loss(
      y = img[, fi$te, drop = FALSE],
      X = cov_mat[fi$te, , drop = FALSE],
      beta = fit$beta, psi = fi$psi, lambda = fi$lambda,
      sigma_fit = fit$sigma, s_ref = fi$s_ref,
      zeta_gamma = zeta_gamma, loss = loss
    )
    msg <- if (!is.finite(err)) {
      sprintf("held-out loss not finite (NA in fitted sigma: %s; NA in beta: %s)",
              anyNA(fit$sigma), anyNA(fit$beta))
    } else NULL
    list(err = if (is.finite(err)) err else NA_real_,
         nonzero = sum(fit$beta != 0), msg = msg)
  }
 
  if (.Platform$OS.type == "windows") n_cores <- 1
  if (verbose) message(sprintf("Fitting %d grid points x %d folds = %d fits on %d core(s)...",
                               n_grid, K, nrow(jobs), n_cores))
  res <- if (n_cores > 1) {
    mclapply(seq_len(nrow(jobs)), run_job, mc.cores = n_cores, mc.preschedule = FALSE)
  } else {
    lapply(seq_len(nrow(jobs)), function(j) {
      out <- run_job(j)
      if (verbose) message(sprintf("  fold %d/%d, zeta = %.3g, zeta0 = %.3g: CV loss = %.5g",
                                   jobs$k[j], K, grid$zeta[jobs$g[j]], grid$zeta0[jobs$g[j]], out$err))
      out
    })
  }
 
  msgs <- unique(unlist(lapply(res, `[[`, "msg")))
  if (length(msgs)) warning("Some fits/losses failed: ", paste(msgs, collapse = " | "))
 
  cv_folds <- matrix(sapply(res, `[[`, "err"), n_grid, K)
  nz_folds <- matrix(sapply(res, `[[`, "nonzero"), n_grid, K)
  if (all(is.na(cv_folds))) stop("All CV losses are NA: ", paste(msgs, collapse = " | "))
 
  cvm  <- rowMeans(cv_folds, na.rm = TRUE)
  cvsd <- apply(cv_folds, 1, sd, na.rm = TRUE) / sqrt(rowSums(!is.na(cv_folds)))
  cv_table <- data.frame(grid, cvm = cvm, cvsd = cvsd,
                         nonzero = rowMeans(nz_folds, na.rm = TRUE))
 
  i_min <- which.min(cvm)
  ok <- which(cvm <= cvm[i_min] + cvsd[i_min])
  ok <- ok[order(-grid$zeta0[ok], -grid$zeta[ok])]
  i_1se <- ok[1]
 
  out <- list(
    cv_table = cv_table,
    cv_folds = cv_folds,
    zeta.min = grid$zeta[i_min], zeta0.min = grid$zeta0[i_min],
    zeta.1se = grid$zeta[i_1se], zeta0.1se = grid$zeta0[i_1se],
    foldid = foldid, loss = loss, rule = rule
  )
 
  ## ---- refit on all subjects ----
  if (refit) {
    i_sel <- if (rule == "min") i_min else i_1se
    if (verbose) message(sprintf("Refitting on all subjects at zeta = %.4g, zeta0 = %.4g (%s rule).",
                                 grid$zeta[i_sel], grid$zeta0[i_sel], rule))
    ols_all <- .ols_fit(img, cov_mat)
    tau_all <- if (is.function(tau)) tau(ols_all$beta, adj$adj_matrix) else tau
    out$fit <- RISEX(
      img = img, voxels = voxels, cov_mat = cov_mat,
      psi_mat = psi_mat, lambda = lambda, adj = adj,
      beta0 = ols_all$beta,
      zeta = grid$zeta[i_sel], zeta0 = grid$zeta0[i_sel],
      tau = tau_all, tau0 = tau0, zeta_gamma = zeta_gamma,
      trace = FALSE, plot = FALSE, ...
    )
  }
 
  class(out) <- "cv.risex"
  out
}
 
 
## Closed-form voxel-wise OLS: beta (V x p) and residuals (V x N)
.ols_fit <- function(img, cov_mat) {
  XtX_inv_Xt <- solve(crossprod(cov_mat), t(cov_mat))   # p x N
  beta <- img %*% t(XtX_inv_Xt)                          # V x p
  list(beta = beta, resid = img - beta %*% t(cov_mat))
}
 
 
## Held-out loss for one fitted model
.cv_heldout_loss <- function(y, X, beta, psi, lambda, sigma_fit, s_ref,
                             zeta_gamma, loss, n_iter = 10) {
  if (length(sigma_fit) == 1) sigma_fit <- rep(sigma_fit, nrow(y))
  r <- y - beta %*% t(X)                                 # V x n_test
  # BLUP-type predictor for latent scores, same form as inside RISEX
  A_alpha <- crossprod(psi / sigma_fit) + diag(1 / lambda, length(lambda))
  A_alpha <- solve(A_alpha, t(psi / sigma_fit^2))       # K x V
  gamma <- matrix(0, nrow(r), ncol(r))
  for (it in seq_len(n_iter)) {
    alpha <- A_alpha %*% (r - gamma)
    e <- r - psi %*% alpha
    gamma_new <- e
    gamma_new[e^2 <= zeta_gamma * sigma_fit^2] <- 0
    if (max(abs(gamma_new - gamma)) < 1e-8) { gamma <- gamma_new; break }
    gamma <- gamma_new
  }
  if (loss == "wmse") {
    # paper (Sec 2.2.3): alpha update without trimming, weights 1 / sigma_fit^2,
    # sum over voxels, average over held-out subjects
    alpha <- A_alpha %*% r
    return(sum(((r - psi %*% alpha) / sigma_fit)^2) / ncol(r))
  }
  alpha <- A_alpha %*% (r - gamma)
  z2 <- ((r - psi %*% alpha) / s_ref)^2
  if (loss == "robust") z2 <- pmin(z2, zeta_gamma)
  mean(z2)
}
 
 
#' @export
print.cv.risex <- function(x, ...) {
  cat("RISE-X cross-validation (", length(unique(x$foldid)), "-fold, ",
      x$loss, " loss)\n", sep = "")
  cat(sprintf("  min rule: zeta = %.4g, zeta0 = %.4g\n", x$zeta.min, x$zeta0.min))
  cat(sprintf("  1se rule: zeta = %.4g, zeta0 = %.4g\n", x$zeta.1se, x$zeta0.1se))
  if (!is.null(x$fit)) cat("  Refit (", x$rule, " rule) stored in $fit\n", sep = "")
  invisible(x)
}
 
 
#' Plot cross-validation curves for RISE-X
#'
#' Plots the mean CV error (with \eqn{\pm 1} SE bars) against \code{zeta0}
#' on a log scale, one curve per value of \code{zeta}.
#'
#' @param x A \code{"cv.risex"} object from \code{\link{cv.RISEX}}.
#' @param ... Unused.
#' @importFrom graphics matplot arrows legend points
#' @export
plot.cv.risex <- function(x, ...) {
  tab <- x$cv_table
  zs <- sort(unique(tab$zeta)); z0 <- sort(unique(tab$zeta0))
  # expand.grid order: zeta varies fastest -> rows = zeta0, cols = zeta after t()
  M  <- t(matrix(tab$cvm,  length(zs), length(z0)))
  S  <- t(matrix(tab$cvsd, length(zs), length(z0)))
  cols <- seq_along(zs)
  matplot(z0, M, type = "b", pch = 19, lty = 1, col = cols, log = "x",
          ylim = range(c(M - S, M + S), na.rm = TRUE),
          xlab = expression(zeta[0]), ylab = "CV error")
  for (j in seq_along(zs)) {
    arrows(z0, M[, j] - S[, j], z0, M[, j] + S[, j], angle = 90, code = 3,
           length = 0.03, col = cols[j])
  }
  points(x$zeta0.min, min(tab$cvm), pch = 4, cex = 2, lwd = 2)
  legend("topleft", legend = sprintf("zeta = %.3g", zs), col = cols, lty = 1,
         pch = 19, bty = "n", cex = 0.8)
  invisible(x)
}
