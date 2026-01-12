#' Voxel-wise initialization for RISE-X
#'
#' Performs voxel-wise linear regression followed by spatially smoothed
#' eigenfunction estimation to initialize model components in the RISE-X
#' framework. Residuals are smoothed using tensor-product splines and
#' decomposed via SVD to obtain low-dimensional latent structure.
#'
#' @param img Numeric matrix of dimension \eqn{V \times N}, where \eqn{V}
#'   is the number of voxels and \eqn{N} is the number of subjects. Each row
#'   corresponds to one voxel.
#' @param cov_mat Numeric matrix of dimension \eqn{N \times p} containing
#'   subject-level covariates. An intercept should be included if desired.
#' @param voxels Numeric matrix of dimension \eqn{V \times 3} giving spatial
#'   coordinates \code{(x, y, z)} for each voxel.
#' @param psi_mat Optional numeric matrix of eigenfunctions. If \code{NULL},
#'   eigenfunctions are estimated from smoothed residuals.
#' @param lambda Optional numeric vector or matrix of eigenvalues corresponding
#'   to \code{psi_mat}. Must be supplied if \code{psi_mat} is provided.
#' @param cut_alpha Proportion of variance explained used to select the number
#'   of eigenfunctions. Default is \code{0.95}.
#'
#' @return A list containing:
#' \describe{
#'   \item{beta}{Estimated voxel-wise regression coefficients (\eqn{V \times p}).}
#'   \item{r_squared}{Vector of voxel-wise \eqn{R^2} values.}
#'   \item{psi_mat}{Estimated eigenfunctions (\eqn{V \times K}).}
#'   \item{lambda}{Eigenvalues corresponding to \code{psi_mat}.}
#'   \item{alpha}{Latent coefficient matrix (\eqn{K \times N}).}
#'   \item{eta}{Low-rank signal component (\eqn{V \times N}).}
#'   \item{eps}{Residual noise component (\eqn{V \times N}).}
#'   \item{sigma0}{Global noise standard deviation.}
#'   \item{sigma1}{Voxel-specific noise standard deviations.}
#' }
#'
#' @details
#' The function first fits voxel-wise linear models to estimate fixed effects
#' and residuals. If eigenfunctions are not supplied, residuals are spatially
#' smoothed using tensor-product splines via \code{\link[mgcv]{gam}} and
#' decomposed using singular value decomposition (SVD).
#'
#' @importFrom stats lm coef residuals sd
#' @importFrom mgcv gam te
#' @export
#'
voxel.wise.fit <- function(img, cov_mat, voxels, psi_mat = NULL, lambda = NULL, cut_alpha = 0.95) {
  # Step 1: Voxel-wise linear model
  fit_1 <- lm(t(img) ~ 0 + cov_mat)
  beta <- t(coef(fit_1))
  eta_est <- t(residuals(fit_1))

  # Step 2: R-squared per voxel
  r_squared <- apply(img, 1, function(v) {
    summary(lm(v ~ cov_mat))$r.squared
  })

  # Step 3: Estimate eigenfunctions if not given
  if (is.null(psi_mat) || is.null(lambda)) {
    message("Estimating eigenfunctions with GAM smoothing...")

    coords_df <- as.data.frame(voxels)
    colnames(coords_df) <- c("x", "y", "z")

    smooth_eta <- matrix(NA, nrow = nrow(eta_est), ncol = ncol(eta_est))

    smooth_single_eta <- function(i) {
      df <- data.frame(eta = eta_est[, i], coords_df)
      mod <- gam(eta ~ te(x, y, z, k = c(5,5,5)),  data = df)
      fitted(mod)
    }

    # Auto-detect cores
    detected_cores <- 1#max(1, parallelly::availableCores() - 1)  # physical cores
    #use_cores <- min(detected_cores, max_cores)

    if (detected_cores > 1) {
      cl <- makeCluster(detected_cores)
      clusterExport(cl, varlist = c("eta_est", "coords_df", "smooth_single_eta"), envir = environment())
      clusterEvalQ(cl, library(mgcv))
      smooth_list <- parLapply(cl, seq_len(ncol(eta_est)), smooth_single_eta)
      stopCluster(cl)
    } else {
      smooth_list <- lapply(seq_len(ncol(eta_est)), smooth_single_eta)
    }

    smooth_eta <- do.call(cbind, smooth_list)

    # Step 4: SVD
    #smooth_eta <- eta_est
    svd_res <- svd(smooth_eta)
    cut_off <- min(which(cumsum(svd_res$d) / sum(svd_res$d) > cut_alpha))
    psi_mat <- svd_res$u[, 1:cut_off, drop = FALSE]
    lambda <- matrix(svd_res$d[1:cut_off]^2, cut_off, 1)
  }
  alpha=solve(t(psi_mat)%*%psi_mat+diag(1/lambda),t(psi_mat)%*%eta_est)
  eta = psi_mat%*%alpha
  eps = eta_est-eta

  sigma0 = sd(c(eps)) #median(abs(eps - median(eps))) #
  sigma1 = apply(eps,1,sd)
  return(list(beta = beta, r_squared = r_squared,
              psi_mat = psi_mat, lambda = lambda,
              alpha = alpha, eta = eta, eps = eps,
              sigma0 = sigma0, sigma1 = sigma1))
}






norm.2 = function(x) {return(sqrt(sum(x^2)))}



#' Construct voxel adjacency matrices from spatial coordinates
#'
#' Builds sparse adjacency matrices encoding neighborhood structure among
#' voxels based on 3D grid coordinates. Neighborhoods are defined using
#' 6-, 18-, or 26-connectivity, commonly used in volumetric neuroimaging
#' analyses. Direction-specific adjacency matrices along the x-, y-, and
#' z-axes are also returned.
#'
#' @param xgrid Numeric matrix or data frame of dimension \eqn{V \times 3},
#'   where each row gives the spatial coordinates \code{(x, y, z)} of a voxel.
#' @param connectivity Integer specifying neighborhood connectivity.
#'   Must be one of \code{6}, \code{18}, or \code{26}. Default is \code{6}.
#'
#' @return A list containing:
#' \describe{
#'   \item{adj_matrix}{Sparse \eqn{V \times V} adjacency matrix encoding all
#'   voxel neighbors under the specified connectivity.}
#'   \item{adj_x}{Sparse adjacency matrix for neighbors differing by \eqn{\pm 1}
#'   in the x-direction only.}
#'   \item{adj_y}{Sparse adjacency matrix for neighbors differing by \eqn{\pm 1}
#'   in the y-direction only.}
#'   \item{adj_z}{Sparse adjacency matrix for neighbors differing by \eqn{\pm 1}
#'   in the z-direction only.}
#' }
#'
#' @details
#' Two voxels are considered neighbors if their coordinate difference matches
#' one of the allowed offsets under the chosen connectivity. The function
#' constructs directed adjacency matrices using sparse representations for
#' computational efficiency in high-dimensional imaging settings.
#'
#' @importFrom Matrix sparseMatrix
#' @importFrom data.table data.table setkey :=
#' @export
coordinates_to_adjacency_matrix.nf <- function(xgrid, connectivity = c(6, 18, 26)) {
  if (!is.matrix(xgrid)) xgrid <- as.matrix(xgrid)
  n <- nrow(xgrid)

  connectivity <- match.arg(as.character(connectivity), choices = c("6", "18", "26"))
  connectivity <- as.integer(connectivity)

  dt <- data.table(x = xgrid[,1], y = xgrid[,2], z = xgrid[,3], idx = 1:n)
  setkey(dt, x, y, z)

  neighbor_offsets_6  <- rbind(
    c(1,0,0), c(-1,0,0),
    c(0,1,0), c(0,-1,0),
    c(0,0,1), c(0,0,-1)
  )

  neighbor_offsets_18 <- rbind(
    neighbor_offsets_6,
    c(1,1,0), c(1,-1,0), c(-1,1,0), c(-1,-1,0),
    c(1,0,1), c(1,0,-1), c(-1,0,1), c(-1,0,-1),
    c(0,1,1), c(0,1,-1), c(0,-1,1), c(0,-1,-1)
  )

  neighbor_offsets_26 <- rbind(
    neighbor_offsets_18,
    c(1,1,1), c(1,1,-1), c(1,-1,1), c(1,-1,-1),
    c(-1,1,1), c(-1,1,-1), c(-1,-1,1), c(-1,-1,-1)
  )

  offsets <- switch(
    as.character(connectivity),
    "6" = neighbor_offsets_6,
    "18" = neighbor_offsets_18,
    "26" = neighbor_offsets_26
  )

  total_i <- integer(0)
  total_j <- integer(0)

  x_i <- integer(0)
  x_j <- integer(0)

  y_i <- integer(0)
  y_j <- integer(0)

  z_i <- integer(0)
  z_j <- integer(0)

  for (k in seq_len(nrow(offsets))) {
    shift <- offsets[k, ]
    neighbors <- dt[, .(x = x + shift[1], y = y + shift[2], z = z + shift[3])]
    neighbors[, nidx := .I]  # index in neighbors

    # Join neighbors with dt on (x,y,z) to get actual indices
    matched <- neighbors[dt, on=.(x,y,z), nomatch=0]
    # 'matched' has:
    # x, y, z from neighbors
    # nidx from neighbors
    # idx from dt (base points)

    # For each neighbor point (nidx), matched.idx is base point idx
    # So edges: from matched.idx (base) to neighbors.nidx

    # Edges from base point (idx) to neighbor point (nidx)
    i_indices <- matched$idx        # base points
    j_indices <- matched$nidx       # neighbor points

    total_i <- c(total_i, i_indices)
    total_j <- c(total_j, j_indices)

    dx <- shift[1]; dy <- shift[2]; dz <- shift[3]

    if (abs(dx) == 1 && dy == 0 && dz == 0) {
      x_i <- c(x_i, i_indices)
      x_j <- c(x_j, j_indices)
    }
    if (abs(dy) == 1 && dx == 0 && dz == 0) {
      y_i <- c(y_i, i_indices)
      y_j <- c(y_j, j_indices)
    }
    if (abs(dz) == 1 && dx == 0 && dy == 0) {
      z_i <- c(z_i, i_indices)
      z_j <- c(z_j, j_indices)
    }
  }

  adj_total <- sparseMatrix(i = total_i, j = total_j, dims = c(n, n), symmetric = FALSE)
  adj_x <- sparseMatrix(i = x_i, j = x_j, dims = c(n, n), symmetric = FALSE)
  adj_y <- sparseMatrix(i = y_i, j = y_j, dims = c(n, n), symmetric = FALSE)
  adj_z <- sparseMatrix(i = z_i, j = z_j, dims = c(n, n), symmetric = FALSE)

  return(list(
    adj_matrix = adj_total,
    adj_x = adj_x,
    adj_y = adj_y,
    adj_z = adj_z
  ))
}




obj.value = function(img,
                     cov_mat,
                     psi_mat,
                     lambda,
                     zeta,
                     tau,
                     zeta0,
                     tau0,
                     zeta_gamma,
                     beta,
                     alpha,
                     sigma,
                     num_obs,
                     num_voxel,
                     neig_id) {
  residual = img - beta %*% t(cov_mat) - psi_mat %*% alpha

  return(
    sum(pmin(residual ^ 2 / sigma ^ 2, zeta_gamma)) / num_obs
    +  zeta * sum(pmin((beta[neig_id[, "i"], ] - beta[neig_id[, "j"], ]) ^ 2, rep(tau, each = dim(neig_id)[1])))
    + sum(rowSums(alpha ^ 2) / lambda) / num_obs
    +  zeta0 * sum(pmin(beta ^ 2 / tau0, 1)) #num_obs *
  )
}



#' RISE-X main fitting routine with adjacency-based spatial regularization
#'
#' Fits the RISE-X model for an image-on-scalar regression problem with
#' (i) voxel-wise robust residual trimming for extreme signals,
#' (ii) low-rank latent structure for subject-specific deviations via a
#' fixed eigenbasis \code{psi_mat}, and
#' (iii) spatial smoothing and sparsity regularization on regression
#' coefficients using voxel adjacency.
#'
#' The algorithm alternates updates of the outlier component \code{gamma},
#' voxel-wise noise scales \code{sigma} (optional), low-rank latent component
#' \code{eta}, and regression coefficients \code{beta} solved via a quadratic
#' subproblem (C++ backend).
#'
#' @param img Numeric matrix of dimension \eqn{V \times N}, where \eqn{V} is the
#'   number of voxels and \eqn{N} is the number of subjects/observations.
#' @param voxels Numeric matrix of dimension \eqn{V \times 3} giving spatial
#'   coordinates \code{(x, y, z)} for each voxel.
#' @param cov_mat Numeric matrix of dimension \eqn{N \times p} containing
#'   subject-level covariates. An intercept should be included if desired.
#' @param psi_mat Numeric matrix of dimension \eqn{V \times K} containing
#'   eigenfunctions/basis functions for the low-rank component.
#' @param lambda Numeric vector (length \eqn{K}) or \eqn{K \times 1} matrix of
#'   eigenvalues corresponding to \code{psi_mat}.
#'
#' @param num_obs Integer; number of subjects/observations \eqn{N}. Defaults to
#'   \code{ncol(img)}.
#' @param num_covs Integer; number of covariates \eqn{p}. Defaults to
#'   \code{ncol(cov_mat)}.
#' @param num_voxel Integer; number of voxels \eqn{V}. Defaults to
#'   \code{nrow(img)}.
#' @param num_lambda Integer; number of eigenfunctions \eqn{K}. Defaults to
#'   \code{length(lambda)}.
#'
#' @param max_it Maximum number of outer iterations. Default is \code{1000}.
#' @param adj A list of sparse adjacency matrices as returned by
#'   \code{coordinates_to_adjacency_matrix()}, containing \code{adj_matrix},
#'   \code{adj_x}, \code{adj_y}, and \code{adj_z}. Defaults to constructing
#'   adjacency from \code{voxels}.
#'
#' @param beta0 Optional numeric matrix of initial regression coefficients
#'   of dimension \eqn{V \times p}. If \code{NULL}, initialized to zeros.
#' @param sigma Optional initial voxel-wise noise scale(s). May be a scalar or a
#'   length-\eqn{V} numeric vector. If \code{NULL}, initialized to \code{1}.
#' @param sigma0 Optional global noise scale used in default regularization
#'   parameters \code{zeta} and \code{zeta0}. Must be non-NULL if you rely on the
#'   defaults \code{zeta = 2 / sigma0^2} and \code{zeta0 = 0.1 / sigma0^2}.
#' @param update_sigma Logical; whether to update voxel-wise \code{sigma} during
#'   iterations using a robust MAD-type estimator. Default is \code{TRUE}.
#'
#' @param zeta Nonnegative smoothing penalty parameter controlling spatial
#'   differences of \code{beta} across neighbors. Default is \code{2 / sigma0^2}.
#' @param tau Baseline smoothness scale used in the smoothing penalty. May be a
#'   scalar or covariate-specific vector. Default uses
#'   \code{calculate_tau(beta0, adj$adj_matrix)}.
#' @param zeta0 Nonnegative sparsity penalty parameter for \code{beta}. Default
#'   is \code{0.1 / sigma0^2}.
#' @param tau0 Positive scale for the sparsity penalty. Default is \code{(1e-3)^2}.
#' @param zeta_gamma Positive threshold parameter controlling outlier detection
#'   via \code{gamma^2 > zeta_gamma * sigma^2}. Default is \code{2.5^2}.
#'
#' @param trace Logical; whether to print iteration diagnostics and timings.
#'   Default is \code{TRUE}.
#' @param plot Logical; whether to generate diagnostic plots during iterations.
#'   Default is \code{FALSE}.
#' @param true_beta Optional \eqn{V \times p} matrix of true coefficients for
#'   simulation diagnostics (RMSE/FNR/FDR). If \code{NULL}, diagnostic metrics are
#'   not returned.
#' @param true_eta Optional \eqn{V \times N} matrix of true low-rank component
#'   for visualization when \code{plot = TRUE}.
#' @param true_sigma Optional length-\eqn{V} vector of true noise scales for
#'   visualization when \code{plot = TRUE}.
#'
#' @param tol.beta Convergence tolerance for maximum absolute change in
#'   \code{beta}. Default is \code{1e-3}.
#' @param tol.obj Convergence tolerance for relative change in objective value.
#'   Default is \code{1e-4}.
#'
#' @return A list containing (at minimum):
#' \describe{
#'   \item{beta}{Estimated regression coefficients (\eqn{V \times p}).}
#'   \item{outlier_detected}{Logical matrix (\eqn{V \times N}) indicating detected
#'     extreme signals in the final iteration.}
#'   \item{alpha}{Estimated latent coefficients (\eqn{K \times N}).}
#'   \item{sigma}{Estimated voxel-wise noise scales (scalar or length-\eqn{V}).}
#'   \item{eta}{Estimated low-rank component (\eqn{V \times N}).}
#'   \item{total_time}{Elapsed time (seconds).}
#'   \item{AIC}{Information criterion based on the fitted noise scales and model
#'     complexity (nonzero \code{beta}).}
#'   \item{BIC}{Bayesian information criterion analogue.}
#'   \item{likelihood}{Log-likelihood proxy used in AIC/BIC calculations.}
#'   \item{adjAIC}{Adjacency-adjusted AIC using effective degrees of freedom
#'     from spatial penalties.}
#'   \item{adjBIC}{Adjacency-adjusted BIC analogue.}
#' }
#' If \code{true_beta} is provided, additional summary metrics are returned:
#' \code{RMSEsig}, \code{FNRsig}, and \code{FDR}.
#'
#' @details
#' The \code{beta} update is carried out by solving a quadratic program using
#' \code{quadratic_cpp()}, which must be available via the package's compiled
#' code. Adjacency matrices are represented using sparse matrices for scalability.
#'
#' @importFrom stats rnorm median sd
#' @importFrom Matrix Matrix Diagonal sparseMatrix kronecker summary
#' @export
smooth.fit.adj = function(img,
                          voxels,
                          cov_mat,
                          psi_mat,
                          lambda,
                          num_obs = ncol(img),
                          num_covs = ncol(cov_mat),
                          num_voxel = nrow(img),
                          num_lambda = length(lambda),
                          max_it = 1000,
                          adj = coordinates_to_adjacency_matrix(voxels),
                          beta0 = NULL,
                          sigma = NULL,
                          sigma0 = NULL,
                          update_sigma = T,
                          zeta = 2 / sigma0^2,
                          tau = calculate_tau(beta0, adj$adj_matrix),  
                          zeta0 = 0.1 / sigma0^2 ,
                          tau0 = (1e-3) ^ 2,
                          zeta_gamma = 2.5 ^ 2,
                          trace = T,
                          plot = F,
                          true_beta = NULL,
                          true_eta = NULL,
                          true_sigma = NULL,
                          #tol = 1e-3,
                          tol.beta = 1e-3,
                          tol.obj = 1e-4) {
  t0 = proc.time()
  time.prepare = 0
  time.gamma = 0
  time.sigma = 0
  time.eta = 0
  time.beta = 0
  time.diag = 0
  temp.time = proc.time()
  obj = Inf

  beta = beta0
  if (is.null(beta)) {
    #num_voxel
    beta = matrix(0, num_voxel, num_covs)
  }
  alpha = matrix(0, num_lambda, num_obs)
  eta = matrix(0, num_voxel, num_obs)
  gamma = matrix(0, num_voxel, num_obs)

  sigma = sigma # get the sigma from voxel wise estimate, each voxel can have different sigma
  if (is.null(sigma)) {
    sigma = 1
  }
  sigma0 = sigma0

  A_alpha = t(psi_mat / sigma) %*% (psi_mat / sigma) + diag(1 / lambda)
  A_alpha = solve(A_alpha, t(psi_mat / sigma ^ 2))

  x = rnorm(10000000)
  x = x[x ^ 2 <= zeta_gamma]
  scaling_sigma = median(abs(x - median(x))) #use median absolute deviation instead of sd, more stable



  ## adjacency matrix
  #adj <- coordinates_to_adjacency_matrix(voxels)
  adj_matrix = adj$adj_matrix
  neig_id <- as.matrix(summary(adj_matrix))
  neig_id <- neig_id[neig_id[, "i"] < neig_id[, "j"], ]

  adj_x <- adj$adj_x
  adj_y <- adj$adj_y
  adj_z <- adj$adj_z
  idx <- as.matrix(summary(adj_x))
  idx <- idx[idx[, "i"] < idx[, "j"], ]

  idy <- as.matrix(summary(adj_y))
  idy <- idy[idy[, "i"] < idy[, "j"], ]

  idz <- as.matrix(summary(adj_z))
  idz <- idz[idz[, "i"] < idz[, "j"], ]

  num_neig <- dim(idx)[1] + dim(idy)[1] + dim(idz)[1]

  #smooth
  ss1 = Matrix(0, nrow = num_covs * dim(idx)[1], ncol = 0, sparse = T)
  tA <- Matrix(0, nrow = dim(idx)[1], ncol = num_voxel, sparse = T)
  tA[cbind(1:dim(idx)[1], idx[,"i"])] <- sqrt(zeta) #sqrt(num_obs * zeta)
  tA[cbind(1:dim(idx)[1], idx[,"j"])] <- -sqrt(zeta) #-sqrt(num_obs * zeta)
  for (q in 1:num_covs) {
    ttA <- Matrix(0, nrow = num_covs * dim(idx)[1], ncol =  num_voxel, sparse = T)
    ttA[((q-1)*dim(idx)[1]) + (1: dim(idx)[1]),] <- tA
    ss1 <- cbind(ss1, ttA)
  }

  ss2 = Matrix(0, nrow = num_covs * dim(idy)[1], ncol = 0, sparse = T)
  tA <- Matrix(0, nrow = dim(idy)[1], ncol = num_voxel, sparse = T)
  tA[cbind(1:dim(idy)[1], idy[,"i"])] <- sqrt(zeta) #sqrt(num_obs * zeta)
  tA[cbind(1:dim(idy)[1], idy[,"j"])] <- -sqrt(zeta) #-sqrt(num_obs * zeta)
  for (q in 1:num_covs) {
    ttA <- Matrix(0, nrow = num_covs * dim(idy)[1], ncol =  num_voxel, sparse = T)
    ttA[(q-1)*dim(idy)[1] + (1: dim(idy)[1]),] <- tA
    ss2 <- cbind(ss2, ttA)
  }

  ss3 = Matrix(0, nrow = num_covs * dim(idz)[1], ncol = 0, sparse = T)
  tA <- Matrix(0, nrow = dim(idz)[1], ncol = num_voxel, sparse = T)
  tA[cbind(1:dim(idz)[1], idz[,"i"])] <- sqrt(zeta) #sqrt(num_obs * zeta)
  tA[cbind(1:dim(idz)[1], idz[,"j"])] <- -sqrt(zeta) #-sqrt(num_obs * zeta)
  for (q in 1:num_covs) {
    ttA <- Matrix(0, nrow = num_covs * dim(idz)[1], ncol =  num_voxel, sparse = T)
    ttA[(q-1)*dim(idz)[1] + (1: dim(idz)[1]),] <- tA
    ss3 <- cbind(ss3, ttA)
  }

  beta_A_smooth = rbind(ss1, ss2, ss3)



  beta_A_sparse = Diagonal(n = num_voxel * num_covs, x = sqrt(zeta0 / tau0)) #x = sqrt(num_obs * zeta0 / tau0)
  beta_A_residual = as(kronecker(-cov_mat, Diagonal(n = num_voxel, x = 1/sigma/sqrt(num_obs))), "dgCMatrix")  # as(kronecker(-cov_mat, Diagonal(n = num_voxel, x = 1/sigma)), "dgCMatrix") / num_obs
  A = rbind(beta_A_residual, beta_A_smooth, beta_A_sparse)
  lam_residual = rep(zeta_gamma / sqrt(num_obs), num_voxel * num_obs) #no change rep(zeta_gamma, num_voxel * num_obs)
  if (length(tau) == 1) {
    lam_smooth = rep(zeta * tau, num_covs * num_neig) #rep(num_obs * zeta * tau, num_covs * num_neig)
  } else {
    lam_smooth = rep(zeta*tau, each = num_neig)
  }

  lam_sparse = rep(zeta0, num_voxel * num_covs) #rep(num_obs * zeta0, num_voxel * num_covs)
  lam = c(lam_residual, lam_smooth, lam_sparse)
  b_smooth = rep(0, num_covs * num_neig)
  b_sparse = rep(0, num_covs * num_voxel)
  b_residual = (img - eta) / sigma /sqrt(num_obs) #(img - eta) / sigma
  b = c(b_residual, b_smooth, b_sparse)

  time.prepare = proc.time() - temp.time

  for (it in 1:max_it) {
    temp.time = proc.time()

    sigma.old = sigma
    beta.old = beta
    obj.old = obj

    residual = img - beta %*% t(cov_mat) - eta
    gamma = residual
    outlier_detected = gamma ^ 2 > zeta_gamma * sigma ^ 2
    gamma[gamma ^ 2 <= zeta_gamma * sigma ^ 2] =  0 #detect outliers report outliers after the last iterations
    #if not outlier, gamma = 0
    time.gamma = time.gamma + proc.time() - temp.time
    temp.time = proc.time()

    if (update_sigma) {
      sigma = apply(residual - gamma, 1, function(x) { #calculate only non-outliers
        y = x[x != 0] ## some rows are all 0 then y will be NA
        return(median(abs(y - median(y))) / scaling_sigma)
      })  ##some sigma is NA then all becomes NA after the next line
      sigma[is.na(sigma)] = mean(sigma, na.rm = T)
      sigma = 0.5 * sigma + 0.5 * mean(sigma) ##sigma: vector of length of voxels




      beta_A_residual = as(kronecker(-cov_mat, Diagonal(n = num_voxel, x = 1/sigma/sqrt(num_obs))), "dgCMatrix") #as(kronecker(-cov_mat, Diagonal(n = num_voxel, x = 1/sigma)), "dgCMatrix")
      A = rbind(beta_A_residual, beta_A_smooth, beta_A_sparse)

      A_alpha = t(psi_mat / sigma) %*% (psi_mat / sigma) + diag(1 / lambda)
      A_alpha = solve(A_alpha, t(psi_mat / sigma ^ 2))

      if (trace) {
        print(norm.2(sigma - sigma.old) / norm.2(sigma.old))
      }

      if (norm.2(sigma - sigma.old) / norm.2(sigma.old) < 0.05) {update_sigma = F}
    }

    time.sigma = time.sigma + proc.time() - temp.time
    temp.time = proc.time()

    eta_est = img - gamma - beta %*% t(cov_mat)
    alpha = A_alpha %*% eta_est
    eta = psi_mat %*% alpha
    b_residual = (img - eta) / sigma /sqrt(num_obs) #(img - eta) / sigma
    b = c(b_residual, b_smooth, b_sparse)


    time.eta = time.eta + proc.time() - temp.time
    temp.time = proc.time()

    beta[1:length(beta)] = quadratic_cpp(
      A = A,
      b = b,
      x0 = c(beta),
      lambda = lam,
      max_it = 1,
      trace = F
    )$par

    time.beta = time.beta + proc.time() - temp.time
    temp.time = proc.time()

    if (trace) {
      obj = obj.value(
        img = img,
        cov_mat = cov_mat,
        psi_mat = psi_mat,
        lambda = lambda,
        zeta = zeta,
        tau = tau,
        zeta0 = zeta0,
        tau0 = tau0,
        zeta_gamma = zeta_gamma,
        beta = beta,
        alpha = alpha,
        sigma = sigma,
        num_obs = num_obs,
        num_voxel = num_voxel,
        neig_id = neig_id
      )
      diagnostic = c(it,
                     mean(sigma),
                     sd(sigma),
                     max(abs(beta - beta.old)),
                     obj,
                     obj - obj.old
      )
      names(diagnostic) = c(
        "iteration",
        "mean(sigma)",
        "sd(sigma)",
        "delta_beta",
        "obj_value",
        "delta_obj_value"
      )
      print(diagnostic)
      if (obj > obj.old && !isTRUE(all.equal(sigma, sigma.old))) {
        warning("objective values increased!")
      }
    }

    beta <- ifelse(abs(beta) > 1e-4, beta, 0)  ####import for calculation of AIC/BIC sum(beta !=0)

    if (plot) {
      plot1 = multi.figs.levelplot(
        cbind(true_beta, beta),
        voxels[, 1],
        voxels[, 2],
        titles = c(
          paste("beta", 1:num_covs, "(v)", sep = ""),
          paste("beta_est", 1:num_covs, "(v)", sep = "")
        ),
        layout = c(num_covs, 2)
      )
      num_eta_show = 8
      plot2 = multi.figs.levelplot(
        cbind(true_eta[, 1:num_eta_show], eta[, 1:num_eta_show]),
        voxels[, 1],
        voxels[, 2],
        titles = c(
          paste("eta", 1:num_eta_show, "(v)", sep = ""),
          paste("eta_est", 1:num_eta_show, "(v)", sep = "")
        ),
        layout = c(num_eta_show, 2)
      )
      plot3 = multi.figs.levelplot(
        cbind(true_sigma, sigma),
        voxels[, 1],
        voxels[, 2],
        titles = c("sigma", "sigma_est"),
        layout = c(2, 1)
      )
      grid.arrange(plot1,
                   plot2,
                   plot3,
                   nrow = 2,
                   ncol = 2)
    }

    time.diag = time.diag + proc.time() - temp.time

    if ((max(abs(beta - beta.old)) < tol.beta)) {
      break
    } else if (abs(obj - obj.old) != Inf) {
      if ((abs(obj - obj.old)/ abs(obj.old + 1e-8) < tol.obj))
        break
    }
  }
  if (trace) {
    print("time.prepare")
    print(time.prepare)
    print("time.gamma")
    print(time.gamma)
    print("time.sigma")
    print(time.sigma)
    print("time.eta")
    print(time.eta)
    print("time.beta")
    print(time.beta)
    print("time.diag")
    print(time.diag)
  }
  total_time = proc.time() - t0


  likelihood = num_obs*sum(log(sigma))
  AIC = 2*num_obs*sum(log(sigma)) + sum(beta !=0) * 2
  BIC = 2*num_obs*sum(log(sigma)) + sum(beta !=0) * log(num_obs * num_voxel)
  adjAIC = 2*num_obs*sum(log(sigma)) + calculate_passtau(beta, adj$adj_matrix, tau)*2
  adjBIC = 2*num_obs*sum(log(sigma)) + calculate_passtau(beta, adj$adj_matrix, tau) * log(num_obs * num_voxel)

  if (!is.null(true_beta)) {
    estimation_mm = metrics(true_beta, beta)
  }

  
  if (is.null(true_beta)) {
    return(list(
      beta = beta,
      outlier_detected = outlier_detected,
      alpha = alpha,
      sigma = sigma,
      eta = eta,
      total_time = total_time[3],
      AIC = AIC,
      BIC = BIC,
      likelihood = likelihood,
      adjAIC = adjAIC,
      adjBIC = adjBIC
    ))
  } else {
    return(list(
      beta = beta,
      outlier_detected = outlier_detected,
      alpha = alpha,
      sigma = sigma,
      eta = eta,
      total_time = total_time[3],
      AIC = AIC,
      BIC = BIC,
      likelihood = likelihood,
      adjAIC = adjAIC,
      adjBIC = adjBIC,
      RMSEsig = mean(estimation_mm$RMSEsig),
      FNRsig = mean(estimation_mm$fnr_sig),
      FDR = mean(estimation_mm$FDR)
    ))
  }


}





#' Performance metrics for estimated voxel-wise regression coefficients
#'
#' Computes signal-region and non-signal-region accuracy metrics for coefficient
#' estimation, including RMSE within signal/non-signal regions, true/false
#' positive/negative rates, and false discovery rate (FDR). Elements of the
#' estimated coefficients with absolute value below a small threshold are set to
#' zero prior to evaluation.
#'
#' @param true_beta Numeric matrix (or coercible to matrix) of true coefficients,
#'   typically of dimension \eqn{V \times p}, where \eqn{V} is the number of voxels
#'   and \eqn{p} is the number of covariates.
#' @param est_beta Numeric matrix (or coercible to matrix) of estimated
#'   coefficients with the same dimension as \code{true_beta}.
#'
#' @return A list with entries (each a numeric vector of length \eqn{p}):
#' \describe{
#'   \item{RMSEsig}{RMSE restricted to voxels in the true signal region (\code{true_beta != 0}).}
#'   \item{RMSEnonsig}{RMSE restricted to voxels in the true non-signal region (\code{true_beta == 0}).}
#'   \item{tpr_sig}{True positive rate in the signal region (1 - FNR).}
#'   \item{fnr_sig}{False negative rate in the signal region.}
#'   \item{fpr_nonsig}{False positive rate in the non-signal region.}
#'   \item{FDR}{False discovery rate, \eqn{FP/(FP+TP)}.}
#' }
#'
#' @details
#' The function thresholds \code{est_beta} at \code{1e-5} in absolute value before
#' computing metrics. Metrics are computed separately for each coefficient
#' column (covariate).
#'
#' Note that if a coefficient column has no true signals (or no true non-signals),
#' some rates may be \code{NaN} due to division by zero.
#'
#' @importFrom stats na.omit
#' @export
metrics <- function(true_beta, est_beta) {
  true_beta <- as.matrix(true_beta)
  est_beta <- as.matrix(est_beta)
  est_beta <- ifelse(abs(est_beta) > 1e-5,est_beta, 0)

  num_beta <- ncol(true_beta)
  signal_region = true_beta != 0

  fnr_sig <- apply((est_beta == 0) * signal_region, 2, sum) / apply(signal_region, 2, sum)
  tpr_sig <- 1 - fnr_sig
  fpr_nonsig <- apply((est_beta != 0) * !signal_region, 2, sum) / apply(!signal_region, 2, sum)

  RMSEsig <- sqrt(apply((est_beta * signal_region - true_beta * signal_region) ^ 2, 2, sum) / apply(signal_region, 2, sum))
  RMSEnonsig <- sqrt(apply((est_beta * !signal_region - true_beta * !signal_region) ^ 2, 2,sum) / apply(!signal_region, 2, sum))

  fp = apply((true_beta == 0) & (est_beta != 0), 2, sum)  #sum((true_beta == 0) & (est_beta != 0))
  tp = apply((true_beta != 0) & (est_beta != 0), 2, sum) #sum((true_beta != 0) & (est_beta != 0))
  FDR = fp / (fp + tp)#apply((est_beta != 0) * !signal_region, 2, sum) / apply((est_beta != 0), 2, sum)

  return(
    list(
      RMSEsig = RMSEsig,
      RMSEnonsig = RMSEnonsig,
      tpr_sig = tpr_sig,
      fnr_sig = fnr_sig,
      fpr_nonsig = fpr_nonsig,
      FDR = FDR
    )
  )
}


