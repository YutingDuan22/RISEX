library(Matrix)

#minimized sum of truncated quadratic functions in 1D
#each summand is weight[i] * min((x - mu[i])^2, lambda[i])
#if specified using coef, then each summand is min(coef[i,1] * x^2 + coef[i,2] * x + coef[i,3], lambda[i])

quadratic.1D = function(coef = NULL,
                        mu = NULL,
                        lambda = Inf,
                        weights = 1) {
    if (!is.null(coef)) {
        a = coef[, 1]
        b = coef[, 2]
        c = coef[, 3]
        weights = a
        mu = -b / (2 * a)
        offset = c - mu ^ 2 * a
        lambda = (lambda - offset) / a
        offset = sum(offset)
    } else {
        offset = 0
    }
    m = length(mu)
    if (length(lambda) == 1)
        lambda = rep(lambda, m)
    if (length(weights) == 1)
        weights = rep(weights, m)
    untruncated = (lambda == Inf)
    truncated = (lambda != Inf & lambda > 0)
    mu.l = mu[truncated] - sqrt(lambda[truncated])
    mu.u = mu[truncated] + sqrt(lambda[truncated])
    end.pts = c(mu.l, mu.u)
    pts.idx = rep(which(truncated), 2)
    s = sort(end.pts, index.return = T)
    pts.idx = pts.idx[s$ix]
    included = untruncated
    w.sum = sum(weights[included])
    w.mu.sum = sum(weights[included] * mu[included])
    w.mu2.sum = sum(weights[included] * mu[included] ^ 2)
    w.lambda.sum = sum(weights[!untruncated] * lambda[!untruncated])
    par = ifelse(w.sum > 0, w.mu.sum / w.sum, 0)
    value = w.mu2.sum - w.sum * par ^ 2 + w.lambda.sum
    for (p in pts.idx) {
        included[p] = !included[p]
        sign = included[p] * 2 - 1
        w.sum = w.sum + sign * weights[p]
        w.mu.sum = w.mu.sum + sign * weights[p] * mu[p]
        w.mu2.sum = w.mu2.sum + sign * weights[p] * mu[p] ^ 2
        w.lambda.sum = w.lambda.sum - sign * weights[p] * lambda[p]
        if (w.sum > 0) {
            x = w.mu.sum / w.sum
            f = w.mu2.sum - w.sum * x ^ 2 + w.lambda.sum
            if (f < value) {
                value = f
                par = x
            }
        }
    }
    return(list(par = par, value = value + offset))
}

obj.value.quadratic.1D <-
    function(x,
             coef = NULL,
             mu = NULL,
             lambda = Inf,
             weights = 1) {
        if (!is.null(coef)) {
            a = coef[, 1]
            b = coef[, 2]
            c = coef[, 3]
            return(sum(pmin(a * x ^ 2 + b * x + c, lambda)))
        } else {
            return(sum(weights * pmin((x - mu) ^ 2, lambda)))
        }
    }

test.quadratic.1D = function() {
    m = 100
    ngrid = 1000
    n.Inf = 5
    mu = rnorm(m)
    x = seq(-3, 3, length.out = ngrid)
    y = rep(0, ngrid)
    lambda = runif(m,-0.05, 0.1)
    weights = runif(m)
    lambda[1:n.Inf] = Inf
    weights[1:n.Inf] = 0.01
    for (i in 1:ngrid) {
        y[i] = obj.value.quadratic.1D(x[i],
                                      mu = mu,
                                      lambda = lambda,
                                      weights = weights)
    }
    print(system.time(sol <-
                          quadratic.1D(
                              mu = mu,
                              lambda = lambda,
                              weights = weights
                          )))
    #system.time(for (i in 1:(100000/m)) {sol = quadratic.1D(mu = mu, lambda = lambda, weights = weights)})
    print(system.time(
        sol1 <-
            quadratic_1D_mu(
                mu = mu,
                lambda = lambda,
                weights = weights
            )
    ))
    #system.time(for (i in 1:(100000/m)) {sol1 = quadratic_1D_mu(mu = mu, lambda = lambda, weights = weights)})
    if (!all.equal(sol, sol1))
        stop("error")
    print(sol)
    par(mar = c(4, 4, 1, 1) + 0.1)
    plot(x,
         y,
         type = 'l',
         xlab = "d",
         ylab = "G'(d)")
    abline(v = sol$par,
           col = 4,
           lty = 2)

    a = runif(m)
    a[1:n.Inf] = 0.01
    b = rnorm(m) * 2 * a
    c = rnorm(m)
    lambda = runif(m,-0.05, 0.1) * a + c - b ^ 2 / 4 / a
    lambda[1:n.Inf] = Inf
    for (i in 1:ngrid) {
        y[i] = obj.value.quadratic.1D(x[i], coef = cbind(a, b, c), lambda = lambda)
    }
    print(system.time(sol <-
                          quadratic.1D(
                              coef = cbind(a, b, c), lambda = lambda
                          )))
    #system.time(for (i in 1:(100000/m)) {sol = quadratic.1D(coef = cbind(a, b, c), lambda = lambda)})
    print(system.time(sol1 <-
                          quadratic_1D_coef(
                              coef = cbind(a, b, c), lambda = lambda
                          )))
    #system.time(for (i in 1:(100000/m)) {sol1 = quadratic_1D_coef(coef = cbind(a, b, c), lambda = lambda)})
    if (!all.equal(sol, sol1))
        stop("error")
    print(sol)
    par(mar = c(4, 4, 1, 1) + 0.1)
    plot(x,
         y,
         type = 'l',
         xlab = "d",
         ylab = "G'(d)")
    abline(v = sol$par,
           col = 4,
           lty = 2)
}

#objective function for sum of truncated quadratic functions in high-D
#each summand is min((A_i x + b_i)^2, lambda[i])
#matrix A can be either dense or sparse
obj.value.quadratic = function(x, A, b = rep(0, nrow(A)), lambda = rep(Inf, ncol(A))) {
    return(sum(pmin((A %*% x + b)^2, lambda)))
}

#minimized sum of truncated quadratic functions in high-D
#each summand is min((A_i x + b_i)^2, lambda[i])
#matrix A can be either dense or sparse
quadratic = function(A, b = rep(0, nrow(A)), x0 = rep(0, ncol(A)), lambda = rep(Inf, nrow(A)), max_it = 1000, tol = 1e-3, trace = F) {
    n = nrow(A)
    p = ncol(A)
    x = x0
    Axb = A %*% x + b
    obj = Inf

    temp_A = lapply(1:p, function(i){
        A_i = A[,i,drop=F]
        index = which(A_i != 0)
        coef = A_i[index]
        return(cbind(index, coef))
    })

    for (it in 1:max_it) {
        x.old = x
        obj.old = obj

        for (i in 1:p) {
            index = temp_A[[i]][,1]
            coef = temp_A[[i]][,2]
            temp_Axb = Axb[index] - coef * x[i]
            x[i] = quadratic_1D_mu(mu = - temp_Axb / coef, lambda = lambda[index] / coef ^ 2, weights = coef ^ 2)$par
            Axb[index] = temp_Axb + coef * x[i]
        }

        if (trace) {
            obj = obj.value.quadratic(x, A, b, lambda)
            diagnostic = c(it,
                           max(abs(x - x.old)),
                           obj,
                           obj - obj.old)
            names(diagnostic) = c(
                "iteration",
                "delta_x",
                "obj_value",
                "delta_obj_value"
            )
            print(diagnostic)
            if (obj > obj.old) {
                warning("objective values increased!")
            }
        }

        if (max(abs(x - x.old)) < tol) {
            break
        }
    }

    obj = obj.value.quadratic(x, A, b, lambda)
    return(list(par = x, value = obj))
}


library(nloptr)

test.quadratic = function() {
    m = 100
    ngrid = 100
    n.Inf = 5
    A = matrix(rnorm(m*2), m, 2)
    A[1:n.Inf, ] = A[1:n.Inf, ] / 10
    b = rnorm(m)
    x = seq(-3, 3, length.out = ngrid)
    y = seq(-3, 3, length.out = ngrid)
    z = matrix(0, ngrid, ngrid)
    lambda = runif(m, -1, 1)
    lambda[1:n.Inf] = Inf
    for (i in 1:ngrid) {
        for (j in 1:ngrid) {
            z[i, j] = obj.value.quadratic(x = c(x[i], y[j]), A = A, b = b, lambda = lambda)
        }
    }
    print(system.time(sol <- quadratic(A = A, b = b, lambda = lambda, tol = 1e-6)))
    print(system.time(sol1 <- quadratic_cpp(A = A, b = b, lambda = lambda, tol = 1e-6)))
    print(system.time(sol2 <- quadratic_cpp(A = as(A, "dgCMatrix"), b = b, lambda = lambda, tol = 1e-6)))
    print(system.time(sol3 <- direct(obj.value.quadratic, lower = c(-3, -3), upper = c(3, 3), control = list(xtol_rel = 1e-6, maxeval = 10000, stopval = -Inf, ftol_rel = 0.0, ftol_abs = 0.0, check_derivatives = F), A = A, b = b, lambda = lambda)))
    A.list = NULL
    b.list = NULL
    for (i in 1:m) {
        A.list = c(A.list, list(2*A[i,]%*%t(A[i,])))
        b.list = c(b.list, list(2*b[i]*A[i,]))
    }
    print(system.time(sol4 <- quadratic.DC(A.list = A.list, b.list = b.list, c = b^2, lambda = lambda, tol = 1e-6)))
    print("CD solution")
    print(sol)
    if (!all.equal(sol, sol1)) stop("error!")
    if (!all.equal(sol, sol2)) stop("error!")
    if (max(abs(sol$par - sol3$par)) > 1e-4 || max(abs(sol$value - sol3$value)) > 1e-4) {
        print("DIRECT solution")
        print(sol3)
    }
    if (max(abs(sol$par - sol4$par)) > 1e-4 || max(abs(sol$value - sol4$value)) > 1e-4) {
        print("DC solution")
        print(sol4)
    }

    require(grDevices)
    filled.contour3(x = x, y = y, z = z, color = terrain.colors, xlim = c(-3,3), ylim = c(-3,3), asp = 1)
    points(x = sol$par[1], y = sol$par[2], cex = 0.5, pch = 3)
    points(x = sol3$par[1], y = sol3$par[2], cex = 0.5, pch = 5)
    points(x = sol4$par[1], y = sol4$par[2], cex = 0.5, pch = 7)
}

#testing
#test.quadratic.1D()
#test.quadratic()
