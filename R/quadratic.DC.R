#DC algorithm for sum of truncated quadratic problems
#there are n summands, each summand is min(1/2x^TA_ix+x^Tb_i+c_i, lambda_i)
quadratic.DC = function(A.list, b.list, c = 0, lambda = 0, x0 = 0, maxit = 1000, tol = 1e-6) {
    p = nrow(A.list[[1]])
    if (length(x0) == 1) {
        x0 = rep(x0, p)
    }
    x = x0
    obj = Inf
    n = length(A.list)
    if (length(c) == 1) {
        c = rep(c, n)
    }
    if (length(lambda) == 1) {
        lambda = rep(lambda, n)
    }
    A = matrix(0, p, p)
    b = rep(0, p)
    for (i in 1:n) {
        A = A + A.list[[i]]
        b = b + b.list[[i]]
    }
    for (it in 1:maxit) {
        x.old = x
        obj.old = obj
        f = sapply(1:n, function(i){return(1/2*t(x)%*%A.list[[i]]%*%x+t(x)%*%b.list[[i]]+c[i])})
        obj = sum(pmin(f, lambda))
        if (obj > obj.old + tol) {
            stop("error")
        }
        truncated = (f > lambda)
        b1 = b
        for (i in which(truncated)) {
            b1 = b1 - A.list[[i]]%*%x - b.list[[i]]
        }
        x = -solve(A, b1)
        if (max(abs(x - x.old)) < tol) {
            break
        }
    }
    return(list(par = x, value = obj))
}
