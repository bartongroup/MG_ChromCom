library(ggplot2)
set.seed(666)

# generate data
x <- seq(100) / 100
y <- 0.5 * x + rnorm(100, sd = 0.03) + 0.2

# function to fit
# in my real case the function is non-analytic
f <- function(x, a, b) {
  a * x + b
}

# error function to minimise: RMS
errfun <- function(par, x, y) {
  a <- par[1]
  b <- par[2]
  err <- sqrt(sum((f(x, a, b) - y)^2))
}

# use optim to fit the model to the data
par <- c(1, 0)
res <- optim(par, errfun, gr=NULL, x, y)

# best-fitting parameters
best_a <- res$par[1]
best_b <- res$par[2]
best_a
best_b



boot_a <- boot_b <- numeric(1000)
for(i in 1:1000){
  j <- sample(100,replace=TRUE)
  x.boot <- x[j]; y.boot <- y[j]
  par <- c(1, 0)
  res <- optim(par, errfun, gr=NULL, x.boot, y.boot)

  # best-fitting parameters
  boot_a[i] <- res$par[1]
  boot_b[i] <- res$par[2]
}

quantile(boot_a, c(0.025, 0.975))
quantile(boot_b, c(0.025, 0.975))

# compare with linear fit
df <- data.frame(
  x = x,
  y = y,
  fy = f(x, best_a, best_b)
)

lf <- lm(y ~ x, data=df)
confint(lf)

# plot the result
ggplot(df, aes(x, y)) +
  geom_point() +
  geom_line(aes(x, fy))
