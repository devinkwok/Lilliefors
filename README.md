Lilliefors test
===============
A quick implementation of Lilliefors test in R, which can be applied to find critical values of any probability distribution function (version 1.0). See references below for more information:

- Lilliefors, H. (1969). On the Kolmogorov-Smirnov Test for the Exponential Distribution with Mean Unknown. Journal of the American Statistical Association, 64(325), 387-389. doi:10.2307/2283748
- Dallal, G., & Wilkinson, L. (1986). An Analytic Approximation to the Distribution of Lilliefors's Test Statistic for Normality. The American Statistician, 40(4), 294-296. doi:10.2307/2684607


Recipes
-------

Get critical values (one-sided) for Lilliefors test:
```
n = 4
rfun = rnorm
reps = 1000000
alphas = c(0.8, 0.85, 0.9, 0.95, 0.99)
lilliefors_crit(4, rfun, reps, alphas)
[1] 0.3029507 0.3215354 0.3454923 0.3753628 0.4133755
```

Get table of critical values:
```
ns = c(4:20, 25, 30)
rfun = rnorm
reps = 1000000
alphas = c(0.8, 0.85, 0.9, 0.95, 0.99)
filename = paste(file.path("tables", "lilliefors-table"), toString(reps), ".csv", sep="")
save_table(ns, rfun, reps, alphas, filename)
```

Compare Dallal & Wilkinson p-value approximation formula against a table:
```
table = read.table(file.path("tables", "lilliefors-table-1000000.csv"), sep=" ", header=TRUE)
nllss_approx_table(table)
```

Plot hypothesized and empirical CDF, and Kolmogorov-Smirnov test statistic:
```
rfun = rnorm
pfun = pnorm
n = 4
plot_ks(rfun(n), pfun)
```
