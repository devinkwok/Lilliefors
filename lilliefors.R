# for testing
source("include.R")

# get mean
mu = function(vec) {
    # n = length(vec)
    # m = sum(vec) / n
    # assert_equal(m, mean(vec))
    # return(m)
    return(mean(vec))
}

# get sd
std_dev = function(vec) {
    # m = mu(vec)
    # n = length(vec)
    # s2 = sum( (vec - m)^2 ) / (n - 1)
    # s = sqrt(s2)
    # assert_equal(s, sd(vec))
    # return(s)
    return(sd(vec))
}

normalize = function(vec) {
    m = mu(vec)
    s = std_dev(vec)
    z = (vec - m) / s
    return(z)
}

double_up = function(vec) {
    paired = c(rbind(vec, vec))
    return(paired)
}

ecdf_x = function(vec) {
    margin = (max(vec) - min(vec)) * 100
    x = sort(vec)
    x = double_up(x)
    x = c(-1 * margin, x, margin)
    return(x)
}

ecdf_y = function(vec) {
    y = 0:length(vec) / length(vec)
    y = double_up(y)
    return(y)
}

ks_max = function(vec) {
    v = sort(vec)
    x = ecdf_x(v)
    y = ecdf_y(v)
    x = x[2:length(x)-1]
    y = y[2:length(y)-1]
    d = abs(pnorm(x) - y)
    t = max(d)
    i = which.max(d)
    return(c(x[i], y[i]))
}

ks_statistic = function(vec) {
    pt = ks_max(vec)
    t = abs(pnorm(pt[1]) - pt[2])
    
    # v = sort(vec)
    # p = 0:length(v) / length(v)
    # low = abs(pnorm(v) - p[1:(length(p)-1)])
    # hi = abs(pnorm(v) - p[2:length(p)])
    # t1 = max(low, hi)
    # assert_equal(t, t1)
    return(t)
}

lilliefors_t = function(vec) {
    z = normalize(vec)
    t = ks_statistic(z)
    return(t)
}

empirical_quantile = function(vec, p) {
    EPSILON = 0.00000001
    v = sort(vec)
    if (p == 0) {
        return(v[1])
    }
    if (p == 1) {
        return(v[n])
    }
    n = length(v)
    i = as.integer(ceiling(n*p))
    if (i - n*p < EPSILON) {
        q = mean(v[i:(i+1)])
        # q1 = (v[i] + v[i+1]) / 2
        # assert_equal(q, q1)
    }
    else {
        q = v[i]
    }
    return(q)
}

lilliefors_crit = function(n, rfun, reps, alphas) {
    x = rfun(n * reps)
    x = split(x, 1:reps)
    t = sapply(x, lilliefors_t)
    t_crit = sapply(alphas, function(a) empirical_quantile(t, a))
    return(t_crit)
}

lilliefors_table = function(ns, rfun, reps, alphas) {
    table = sapply(ns, function(n) lilliefors_crit(n, rfun, reps, alphas))
    table = t(table)
    colnames(table) = alphas
    rownames(table) = ns
    return(table)
}

save_table = function(ns, rfun, reps, alphas, filename) {
    table = lilliefors_table(ns, rfun, reps, alphas)
    write.table(table, file=filename, sep=" ", header=TRUE)
}

plot_ks = function(vec, pfun) {
    MARGIN = 0.8
    PTS = 200

    z = sort(normalize(vec))
    pt = ks_max(z)
    xd = pt[1]
    yd = pt[2]
    d = pnorm(xd) - yd

    a = min(z)
    b = max(z)
    r = b - a
    extra = r * MARGIN
    x = seq(a - extra, b + extra, length.out=PTS)
    y = pfun(x)

    plot(x, y, type="l", lwd=6, col="black", xlab="Z score", ylab="Cumulative probability", main="Calculating Kolmogorov-Smirnov statistic")
    lines(ecdf_x(z), ecdf_y(z), lwd=1, cex=0, col="red")
    segments(xd,yd, y1=(yd + d), col="black", lty="dotted", lwd=3)
    legend("bottomright", lwd=c(6,2, 3), lty=c("solid", "solid", "dotted"), col=c("black", "red", "black"), legend=c("Hypothesized CDF", "Empirical CDF", "K-S test statistic"))
}

# p less than 0.100, sample size 5: 100
nllss_approx = function(n, dmax) {
    if (n > 100) {
        dmax = dmax * (n / 100) ^ 0.49
        n = 100
    }
    p = exp(-7.01256 * dmax^2 * (n + 2.78019)
        + 2.99587 * dmax * sqrt(n + 2.78019)
        - 0.122119 + .974598 / sqrt(n) + 1.67997 / n)
    return(p)
}

nllss_approx_table = function(table) {
    ns = as.integer(rownames(table))
    analytic = NULL
    for (i in 1:nrow(table)) {
        n = ns[i]
        row = table[i,]
        approx = sapply(row, function(dmax) nllss_approx(n, dmax))
        analytic = rbind(analytic, approx)
    }
    rownames(analytic) = rownames(table)
    return(analytic)
}

# recipes

# # to make a plot of K-S statistic
#
# rfun = rnorm
# pfun = pnorm
# n = 4
# plot_ks(rfun(n), pfun)

# # to get crit values
#
# n = 4
# rfun = rnorm
# reps = 1000000
# alphas = c(0.8, 0.85, 0.9, 0.95, 0.99)
# lilliefors_crit(4, rfun, reps, alphas)
# [1] 0.3029507 0.3215354 0.3454923 0.3753628 0.4133755

# # to get the table
#
# ns = c(4:20, 25, 30)
# rfun = rnorm
# reps = 1000000
# alphas = c(0.8, 0.85, 0.9, 0.95, 0.99)
# filename = paste(file.path("tables", "lilliefors-table"), toString(reps), ".csv", sep="")
# save_table(ns, rfun, reps, alphas, filename)

# # to compare approximations against table
#
# table = read.table(file.path("tables", "lilliefors-table-1000000.csv"), sep=" ", header=TRUE)
# nllss_approx_table(table)