n_person <- 3000
s <- rnorm(n_person)
lambda <- c(0.7, 0.7, 0.4)
sigma_res <- diag(1 - diag(lambda %*% t(lambda)))
res <- MASS::mvrnorm(n_person, mu = rep(0, 3), Sigma = sigma_res)
w <- outer(s, lambda, "*") + res
########## estimate loadings #########
r_w <- cor(w)
lambda_hat <- rep(NA, 3)
lambda_hat[c(1, 2)] <- sqrt(r_w[1, 2])
lambda_hat[3] <- 0.5 * (r_w[3, 1] + r_w[3, 2]) / lambda[1]
lambda_hat <- t(t(lambda_hat))
#######################################
########## estimate error covariance matrix###
sigma_w <- cov(w)
sigma_res <- sigma_w - lambda_hat %*% t(lambda_hat)
blup_coef <- solve(sigma_res) %*% lambda_hat / (t(lambda_hat) %*% solve
  (sigma_res) %*% lambda_hat)[1]
blp_coef <- solve(sigma_res) %*% lambda_hat / (1 + (t(lambda_hat) %*% solve
  (sigma_res) %*% lambda_hat)[1])
######################################
s_hat <- w %*% blup_coef
plot(s_hat ~ s)
abline(a = 0, b = 1)
######################################
s_hat <- w %*% blp_coef
points(s_hat ~ s, col = "red")
abline(a = 0, b = 1)