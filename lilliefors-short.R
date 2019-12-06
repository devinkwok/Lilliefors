# this code was included in report (version 1.0)
normalize = function(vec) {
  m = mean(vec)
  s = sd(vec)
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
  if (p == 0) { return(v[1]) }
  if (p == 1) { return(v[n]) }
  n = length(v)
  i = as.integer(ceiling(n*p))
  if (i - n*p < EPSILON) { q = mean(v[i:(i+1)]) }
  else { q = v[i] }
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
