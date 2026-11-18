#' K-fold cross-validation for tuning \code{zeta} and \code{zeta0} in RISE-X
#'
#' Selects the smoothness penalty \code{zeta} and the sparsity penalty
#' \code{zeta0} of \code{\link{RISEX}} by K-fold cross-validation over
#' subjects. Within each fold, the procedure is:
#' \enumerate{
#'   \item fit \code{\link{voxel.wise.fit}} on the training subjects to get
#'     \code{psi_mat}, \code{lambda}, the initial \code{beta0}, \code{sigma}
#'     (\code{sigma1}) and \code{sigma0};
#'   \item set \code{tau = calculate_tau(beta0, adj$adj_matrix, c_tau)};
#'   \item fit \code{RISEX} on the training subjects for every
#'     \code{(zeta, zeta0)} on the grid;
#'   \item for the test subjects, compute \eqn{\hat\eta = Y - \hat\beta X},
#'     \eqn{\hat\alpha = (\Psi^\top\Psi + \Lambda^{-1})^{-1}\Psi^\top\hat\eta},
#'     the prediction \eqn{\hat Y = \hat\beta X + \Psi\hat\alpha}, and the
#'     CV errors.
#' }
#' The voxel-wise fit in step 1 does not depend on \code{(zeta, zeta0)}, so it
#' is computed once per fold and shared across the grid.
#'
#' @param img Numeric matrix (\eqn{V \times N}) of images, one column per subject.
#' @param voxels Numeric matrix (\eqn{V \times 3}) of integer grid coordinates.
#' @param cov_mat Numeric matrix (\eqn{N \times p}) of covariates.
#' @param zeta_grid Numeric vector of candidate \code{zeta} values.
#' @param zeta0_grid Numeric vector of candidate \code{zeta0} values.
#' @param nfolds Number of folds. Default \code{10}.
#' @param foldid Optional fold assignment (length \eqn{N}); overrides
#'   \code{nfolds} and \code{seed}.
#' @param seed Seed used right before drawing the folds, as
#'   \code{sample(rep(1:nfolds, length.out = N))}. Default \code{123}.
#' @param adj Adjacency list from \code{\link{coordinates_to_adjacency_matrix.nf}}.
#'   Built from \code{voxels} with \code{connectivity} if \code{NULL}.
#' @param connectivity Connectivity used when \code{adj} is \code{NULL}.
#' @param c_tau Multiplier passed to \code{\link{calculate_tau}}. Default \code{0.9}.
#' @param zeta_gamma Outlier threshold passed to \code{RISEX}. Default \code{2.5^2}.
#' @param psi_mat,lambda Optional fixed eigenbasis. If \code{NULL} (default),
#'   it is re-estimated on each training fold by \code{voxel.wise.fit}
#'   (GAM smoothing + SVD).
#' @param select_by CV error used to choose the tuning values:
#'   \describe{
#'     \item{\code{"wmse"}}{(default) \code{mean((Y - Yhat)^2 / sigma^2)},
#'       where \code{sigma} is the voxel-wise noise scale from the training
#'       \code{RISEX} fit.}
#'     \item{\code{"mse"}}{\code{mean((Y - Yhat)^2)}.}
#'     \item{\code{"robust"}}{\code{mean(pmin((Y - Yhat)^2 / sigma1^2, zeta_gamma))},
#'       with \code{sigma1} from the training-fold \code{voxel.wise.fit}, which
#'       does not change across the grid.}
#'   }
#'   All three are always computed and returned.
#' @param rule \code{"min"} (lowest CV error) or \code{"1se"} (the most
#'   regularized pair within one standard error of the minimum; ties broken by
#'   larger \code{zeta0}, then larger \code{zeta}).
#' @param refit Logical; refit on all subjects at the chosen pair. Default \code{TRUE}.
#' @param n_cores Cores for \code{parallel::mclapply} over (fold, grid) fits.
#'   Forced to \code{1} on Windows.
#' @param verbose Logical; print progress.
#' @param ... Further arguments passed to \code{RISEX}, e.g. \code{max_it},
#'   \code{tol.beta}, \code{tol.obj}.
#'
#' @return An object of class \code{"cv.risex"}:
#' \describe{
#'   \item{cv_table}{Data frame with one row per \code{(zeta, zeta0)}: mean
#'     and SE across folds of \code{wmse}, \code{mse} and \code{robust}, and the
#'     average number of nonzero coefficients.}
#'   \item{cv_mse, cv_mseSigma, cv_robust}{Matrices (grid \eqn{\times} folds)
#'     of per-fold errors; \code{cv_mseSigma} is the WMSE.}
#'   \item{zeta.min, zeta0.min, zeta.1se, zeta0.1se}{Selected pairs.}
#'   \item{foldid}{Fold assignment used.}
#'   \item{fit}{Full-data \code{RISEX} fit at the chosen pair (if \code{refit}).}
#' }
#'
#' @examples
#' \dontrun{
#' adj <- coordinates_to_adjacency_matrix.nf(voxels, connectivity = 6)
#' cv  <- cv.RISEX(img, voxels, cov_mat,
#'                 zeta_grid = 0.3,
#'                 zeta0_grid = c(0.01, 0.03, 0.05, 0.1, 0.15, 0.2, 0.3),
#'                 adj = adj, nfolds = 10, n_cores = 8)
#' cv$cv_table
#' plot(cv)
#' }
#'
#' @importFrom stats sd
#' @importFrom parallel mclapply
#' @export
cv.RISEX <- function(img,
                     voxels,
                     cov_mat,
                     zeta_grid,
                     zeta0_grid,
                     nfolds = 10,
                     foldid = NULL,
                     seed = 123,
                     adj = NULL,
                     connectivity = 6,
                     c_tau = 0.9,
                     zeta_gamma = 2.5^2,
                     psi_mat = NULL,
                     lambda = NULL,
                     select_by = c("wmse", "mse", "robust"),
                     rule = c("min", "1se"),
                     refit = TRUE,
                     n_cores = 1,
                     verbose = TRUE,
                     ...) {
  select_by <- match.arg(select_by)
  rule <- match.arg(rule)
  img <- as.matrix(img)
  cov_mat <- as.matrix(cov_mat)
  voxels <- as.matrix(voxels)
  N <- ncol(img)
  V <- nrow(img)

  ## ---- input checks (common causes of NA CV errors) ----
  if (nrow(cov_mat) != N) stop("nrow(cov_mat) must equal ncol(img).")
  if (nrow(voxels) != V) stop("nrow(voxels) must equal nrow(img).")
  if (anyNA(img)) stop("`img` contains NA/NaN values (e.g. out-of-mask voxels). Remove those voxels (rows) first.")
  if (anyNA(cov_mat)) stop("`cov_mat` contains NA values. Remove or impute those subjects first.")
  if (!is.null(psi_mat) && (is.null(lambda) || nrow(psi_mat) != V)) {
    stop("If `psi_mat` is supplied, `lambda` must be too and nrow(psi_mat) must equal nrow(img).")
  }
  if ("tol" %in% names(list(...))) {
    stop("RISEX() has no `tol` argument; use `tol.beta` and/or `tol.obj`.")
  }

  if (is.null(adj)) adj <- coordinates_to_adjacency_matrix.nf(voxels, connectivity = connectivity)
  if (sum(adj$adj_x) + sum(adj$adj_y) + sum(adj$adj_z) == 0) {
    stop("No neighbouring voxels found. `voxels` must be integer grid indices ",
         "(neighbours differ by 1), not mm/MNI coordinates.")
  }

  grid <- expand.grid(zeta = sort(unique(zeta_grid)), zeta0 = sort(unique(zeta0_grid)))
  n_grid <- nrow(grid)

  ## ---- folds over subjects (same draw as the simulation scripts) ----
  if (is.null(foldid)) {
    if (nfolds < 2 || nfolds > N) stop("nfolds must be between 2 and N.")
    if (!is.null(seed)) set.seed(seed)
    foldid <- sample(rep(seq_len(nfolds), length.out = N))
  } else {
    if (length(foldid) != N) stop("foldid must have length ncol(img).")
    foldid <- as.integer(factor(foldid))
  }
  K <- max(foldid)

  ## ---- step 1-2: voxel-wise fit + tau, once per fold ----
  if (verbose) message(sprintf("Voxel-wise initialization on %d training folds...", K))
  fold_info <- lapply(seq_len(K), function(k) {
    tr <- which(foldid != k)
    te <- which(foldid == k)
    vw <- suppressMessages(voxel.wise.fit(img[, tr, drop = FALSE], cov_mat[tr, , drop = FALSE],
                                          voxels, psi_mat = psi_mat, lambda = lambda))
    list(tr = tr, te = te,
         psi = vw$psi_mat, lambda = as.numeric(vw$lambda),
         beta0 = vw$beta, sigma1 = vw$sigma1, sigma0 = vw$sigma0,
         tau = calculate_tau(vw$beta, adj$adj_matrix, c_tau))
  })

  ## ---- step 3-4: RISEX fit + held-out errors for every (fold, grid) ----
  jobs <- expand.grid(g = seq_len(n_grid), k = seq_len(K))
  na_out <- function(msg) list(wmse = NA_real_, mse = NA_real_, robust = NA_real_,
                               nonzero = NA_real_, msg = msg)
  run_job <- function(j) {
    g <- jobs$g[j]; k <- jobs$k[j]
    fi <- fold_info[[k]]
    fit <- tryCatch(
      RISEX(img = img[, fi$tr, drop = FALSE], voxels = voxels,
            cov_mat = cov_mat[fi$tr, , drop = FALSE],
            psi_mat = fi$psi, lambda = fi$lambda, adj = adj,
            beta0 = fi$beta0, sigma = fi$sigma1, sigma0 = fi$sigma0,
            zeta = grid$zeta[g], zeta0 = grid$zeta0[g], tau = fi$tau,
            zeta_gamma = zeta_gamma, trace = FALSE, plot = FALSE, ...),
      error = function(e) e)
    if (inherits(fit, "error")) return(na_out(paste("RISEX error:", conditionMessage(fit))))

    err <- .cv_errors(y = img[, fi$te, drop = FALSE], X = cov_mat[fi$te, , drop = FALSE],
                      beta = fit$beta, psi = fi$psi, lambda = fi$lambda,
                      sigma_fit = fit$sigma, sigma1 = fi$sigma1, zeta_gamma = zeta_gamma)
    msg <- if (!all(is.finite(unlist(err)))) {
      sprintf("non-finite CV error (NA in fitted sigma: %s; NA in beta: %s)",
              anyNA(fit$sigma), anyNA(fit$beta))
    } else NULL
    c(err, list(nonzero = sum(fit$beta != 0), msg = msg))
  }

  if (.Platform$OS.type == "windows") n_cores <- 1
  if (verbose) message(sprintf("Fitting RISEX: %d grid points x %d folds = %d fits on %d core(s)...",
                               n_grid, K, nrow(jobs), n_cores))
  res <- if (n_cores > 1) {
    mclapply(seq_len(nrow(jobs)), run_job, mc.cores = n_cores, mc.preschedule = FALSE)
  } else {
    lapply(seq_len(nrow(jobs)), function(j) {
      out <- run_job(j)
      if (verbose) message(sprintf("  fold %d/%d, zeta = %g, zeta0 = %g: WMSE = %.5g, MSE = %.5g",
                                   jobs$k[j], K, grid$zeta[jobs$g[j]], grid$zeta0[jobs$g[j]],
                                   out$wmse, out$mse))
      out
    })
  }
  # mclapply returns a try-error object if a worker dies
  res <- lapply(res, function(r) if (inherits(r, "try-error")) na_out(as.character(r)) else r)

  msgs <- unique(unlist(lapply(res, `[[`, "msg")))
  if (length(msgs)) warning("Some CV fits failed: ", paste(msgs, collapse = " | "))

  to_mat <- function(name) matrix(sapply(res, `[[`, name), n_grid, K)
  cv_mseSigma <- to_mat("wmse")
  cv_mse      <- to_mat("mse")
  cv_robust   <- to_mat("robust")
  nz_folds    <- to_mat("nonzero")
  sel_mat <- switch(select_by, wmse = cv_mseSigma, mse = cv_mse, robust = cv_robust)
  if (all(is.na(sel_mat))) stop("All CV errors are NA: ", paste(msgs, collapse = " | "))

  m  <- function(M) rowMeans(M, na.rm = TRUE)
  se <- function(M) apply(M, 1, sd, na.rm = TRUE) / sqrt(rowSums(!is.na(M)))
  cv_table <- data.frame(grid,
                         wmse = m(cv_mseSigma), wmse_se = se(cv_mseSigma),
                         mse = m(cv_mse), mse_se = se(cv_mse),
                         robust = m(cv_robust), robust_se = se(cv_robust),
                         nonzero = m(nz_folds))

  cvm <- m(sel_mat); cvsd <- se(sel_mat)
  i_min <- which.min(cvm)
  ok <- which(cvm <= cvm[i_min] + cvsd[i_min])
  ok <- ok[order(-grid$zeta0[ok], -grid$zeta[ok])]
  i_1se <- ok[1]

  out <- list(cv_table = cv_table,
              cv_mse = cv_mse, cv_mseSigma = cv_mseSigma, cv_robust = cv_robust,
              zeta.min = grid$zeta[i_min], zeta0.min = grid$zeta0[i_min],
              zeta.1se = grid$zeta[i_1se], zeta0.1se = grid$zeta0[i_1se],
              foldid = foldid, select_by = select_by, rule = rule)

  ## ---- refit on all subjects ----
  if (refit) {
    i_sel <- if (rule == "min") i_min else i_1se
    if (verbose) message(sprintf("Refitting on all subjects at zeta = %g, zeta0 = %g (%s, %s rule).",
                                 grid$zeta[i_sel], grid$zeta0[i_sel], select_by, rule))
    vw <- suppressMessages(voxel.wise.fit(img, cov_mat, voxels, psi_mat = psi_mat, lambda = lambda))
    out$voxel_wise_fit <- vw
    out$fit <- RISEX(img = img, voxels = voxels, cov_mat = cov_mat,
                     psi_mat = vw$psi_mat, lambda = as.numeric(vw$lambda), adj = adj,
                     beta0 = vw$beta, sigma = vw$sigma1, sigma0 = vw$sigma0,
                     zeta = grid$zeta[i_sel], zeta0 = grid$zeta0[i_sel],
                     tau = calculate_tau(vw$beta, adj$adj_matrix, c_tau),
                     zeta_gamma = zeta_gamma, trace = FALSE, plot = FALSE, ...)
  }

  class(out) <- "cv.risex"
  out
}


## Held-out CV errors for one fitted model (matches the simulation script)
.cv_errors <- function(y, X, beta, psi, lambda, sigma_fit, sigma1, zeta_gamma) {
  lambda <- as.numeric(lambda)
  if (length(sigma_fit) == 1) sigma_fit <- rep(sigma_fit, nrow(y))
  eta_est <- y - beta %*% t(X)
  alpha <- solve(crossprod(psi) + diag(1 / lambda, length(lambda)), crossprod(psi, eta_est))
  r <- y - (beta %*% t(X) + psi %*% alpha)
  list(wmse   = mean(r^2 / sigma_fit^2),
       mse    = mean(r^2),
       robust = mean(pmin(r^2 / sigma1^2, zeta_gamma)))
}


#' @export
print.cv.risex <- function(x, ...) {
  cat("RISE-X cross-validation (", length(unique(x$foldid)), "-fold, selected by ",
      x$select_by, ")\n", sep = "")
  cat(sprintf("  min rule: zeta = %g, zeta0 = %g\n", x$zeta.min, x$zeta0.min))
  cat(sprintf("  1se rule: zeta = %g, zeta0 = %g\n", x$zeta.1se, x$zeta0.1se))
  if (!is.null(x$fit)) cat("  Refit (", x$rule, " rule) stored in $fit\n", sep = "")
  invisible(x)
}


#' Plot cross-validation curves for RISE-X
#'
#' Plots the mean CV error (\eqn{\pm 1} SE) against \code{zeta0}, one curve per
#' \code{zeta}.
#'
#' @param x A \code{"cv.risex"} object from \code{\link{cv.RISEX}}.
#' @param which Which error to plot: \code{"wmse"}, \code{"mse"} or
#'   \code{"robust"}. Defaults to the one used for selection.
#' @param ... Unused.
#' @importFrom graphics matplot arrows legend points
#' @export
plot.cv.risex <- function(x, which = x$select_by, ...) {
  tab <- x$cv_table
  zs <- sort(unique(tab$zeta)); z0 <- sort(unique(tab$zeta0))
  # expand.grid order: zeta varies fastest -> rows = zeta0, cols = zeta after t()
  M <- t(matrix(tab[[which]], length(zs), length(z0)))
  S <- t(matrix(tab[[paste0(which, "_se")]], length(zs), length(z0)))
  cols <- seq_along(zs)
  matplot(z0, M, type = "b", pch = 19, lty = 1, col = cols,
          ylim = range(c(M - S, M + S), na.rm = TRUE),
          xlab = expression(zeta[0]), ylab = paste("CV", which))
  for (j in seq_along(zs)) {
    arrows(z0, M[, j] - S[, j], z0, M[, j] + S[, j], angle = 90, code = 3,
           length = 0.03, col = cols[j])
  }
  i <- which.min(tab[[which]])
  points(tab$zeta0[i], tab[[which]][i], pch = 4, cex = 2, lwd = 2)
  legend("topleft", legend = sprintf("zeta = %g", zs), col = cols, lty = 1,
         pch = 19, bty = "n", cex = 0.8)
  invisible(x)
}


