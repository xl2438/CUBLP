source("blup.R")
genearte_data <- function(lambda, n_person) {
  ############# generate data #########
  s <- rnorm(n_person)
  sigma_res <- diag(1 - diag(lambda %*% t(lambda)))
  res <- MASS::mvrnorm(n_person, mu = rep(0, length(lambda)), Sigma = sigma_res)
  w <- outer(s, lambda, "*") + res
  return(w)
}
##############################################

lambda <- matrix(c(0.8, 0.8, 0.6, 0.3,
                   0.5, 0.5, 0.6, 0.3), nrow = 2, byrow = TRUE)
n_person <- seq(from = 500, to = 3000, by = 500)
n_rep <- 1000
mat_lambda_hat <- array(NA, dim = c(4, dim(lambda)[1], length(n_person), n_rep))
for (i in seq_len(dim(lambda)[1])) {
  for (j in seq_len(length(n_person))) {
    for (k in seq_len(n_rep)) {
      w <- genearte_data(lambda[i, ], n_person[j])
      fit <- compute_coef(cor(w), 2, 2)
      mat_lambda_hat[, i, j, k] <- fit$lambda_hat
    }
  }
}
##################################################
########## plotting results #################
lambda_hat_mean <- apply(mat_lambda_hat, c(1, 2, 3), mean)
png("bias_equal_discrim_high.png", width = 6, height = 6, units = "in", res =
  800)
plot(lambda_hat_mean[1, 1, ] - 0.8 ~ n_person, type = "l", xlab = "sample size",
 ylab = "bias", ylim = c(-0.005, 0.005), lty = 2)
for (i in c(3:4)) {
  lines(n_person, lambda_hat_mean[i, 1, ] - lambda[1, i], lty = i)
}
abline(h = 0)
legend("topright", legend = c(expression(hat(lambda[Y])), expression(hat(lambda
  [X[1]])), expression(hat(lambda[X[2]]))), lty = c(2:4))
dev.off()
#############
png("bias_equal_discrim_low.png", width = 6, height = 6, units = "in", res =
  800)
plot(lambda_hat_mean[1, 2, ] - 0.5 ~ n_person, type = "l", xlab = "sample size",
 ylab = "bias", ylim = c(-0.005, 0.005), lty = 2)
for (i in c(3:4)) {
  lines(n_person, lambda_hat_mean[i, 2, ] - lambda[2, i], lty = i)
}
abline(h = 0)
legend("topright", legend = c(expression(hat(lambda[Y])), expression(hat(lambda
  [X[1]])), expression(hat(lambda[X[2]]))), lty = c(2:4))
dev.off()
#############################
png("var_equal_discrim_high.png", width = 6, height = 6, units = "in", res =
  800)
lambda_hat_var <- apply(mat_lambda_hat, c(1, 2, 3), var)
plot(lambda_hat_var[1, 1, ] ~ n_person, type = "l", xlab = "sample size", ylab
  = "variance", ylim = c(0, 0.005), lty = 2)
lines(lambda_hat_var[3, 1, ] ~ n_person, lty = 3)
lines(lambda_hat_var[4, 1, ] ~ n_person, lty = 4)
legend("topright", legend = c(expression(hat(lambda[Y])), expression(hat(lambda
  [X[1]])), expression(hat(lambda[X[2]]))), lty = c(2:4))
dev.off()
######################
png("var_equal_discrim_low.png", width = 6, height = 6, units = "in", res =
  800)
lambda_hat_var <- apply(mat_lambda_hat, c(1, 2, 3), var)
plot(lambda_hat_var[1, 2, ] ~ n_person, type = "l", xlab = "sample size", ylab
  = "variance", ylim = c(0, 0.005), lty = 2)
lines(lambda_hat_var[3, 2, ] ~ n_person, lty = 3)
lines(lambda_hat_var[4, 2, ] ~ n_person, lty = 4)
legend("topright", legend = c(expression(hat(lambda[Y])), expression(hat(lambda
  [X[1]])), expression(hat(lambda[X[2]]))), lty = c(2:4))
dev.off()
#########################################################
####################### unequal loadings ###############
lambda <- matrix(c(0.9, 0.7, 0.6, 0.3,
                   0.6, 0.4, 0.6, 0.3), nrow = 2, byrow = TRUE)
n_person <- seq(from = 500, to = 3000, by = 500)
n_rep <- 1000
mat_lambda_hat <- array(NA, dim = c(4, dim(lambda)[1], length(n_person), n_rep))
for (i in seq_len(dim(lambda)[1])) {
  for (j in seq_len(length(n_person))) {
    for (k in seq_len(n_rep)) {
      w <- genearte_data(lambda[i, ], n_person[j])
      fit <- compute_coef_uneq_loadings(cor(w), 2, 2)
      mat_lambda_hat[, i, j, k] <- fit$lambda_hat
    }
  }
}
########## plotting results #################
lambda_hat_mean <- apply(mat_lambda_hat, c(1, 2, 3), mean)
png("bias_uneq_discrim_high.png", width = 6, height = 6, units = "in", res =
  800)
plot(lambda_hat_mean[1, 1, ] - lambda[1, 1] ~ n_person, type = "l", xlab =
 "sample size", ylab = "bias", ylim = c(-0.005, 0.005), lty = 2)
for (i in c(2:4)) {
  lines(n_person, lambda_hat_mean[i, 1, ] - lambda[1, i], lty = i + 1)
}
abline(h = 0)
my_lengend <- c(expression(hat(lambda[Y[1]])), expression(hat(lambda[Y[2]])),
 expression(hat(lambda[X[1]])), expression(hat(lambda[X[2]])))
legend("topright", legend = my_lengend, lty = c(2:5))
dev.off()
#############
png("bias_uneq_discrim_low.png", width = 6, height = 6, units = "in", res =
  800)
plot(lambda_hat_mean[1, 2, ] - lambda[2, 1] ~ n_person, type = "l", xlab =
 "sample size", ylab = "bias", ylim = c(-0.005, 0.005), lty = 2)
for (i in c(2:4)) {
  lines(n_person, lambda_hat_mean[i, 2, ] - lambda[2, i], lty = i + 1)
}
abline(h = 0)
legend("topright", legend = my_lengend, lty = c(2:5))
dev.off()
#############################
png("var_uneq_discrim_high.png", width = 6, height = 6, units = "in", res =
  800)
lambda_hat_var <- apply(mat_lambda_hat, c(1, 2, 3), var)
plot(lambda_hat_var[1, 1, ] ~ n_person, type = "l", xlab = "sample size", ylab
  = "variance", ylim = c(0, 0.005), lty = 2)
for (i in 2:4) {
  lines(lambda_hat_var[i, 1, ] ~ n_person, lty = i + 1)
}
legend("topright", legend = my_lengend, lty = c(2:5))
dev.off()
######################
png("var_uneq_discrim_low.png", width = 6, height = 6, units = "in", res =
  800)
lambda_hat_var <- apply(mat_lambda_hat, c(1, 2, 3), var)
plot(lambda_hat_var[1, 2, ] ~ n_person, type = "l", xlab = "sample size", ylab
  = "variance", ylim = c(0, 0.005), lty = 2)
for (i in 2:4) {
  lines(lambda_hat_var[i, 2, ] ~ n_person, lty = i + 1)
}
legend("topright", legend = my_lengend, lty = c(2:5))
dev.off()
#####################################################
##################prediction bias####################
lambda <- matrix(c(0.8, 0.8, 0.6, 0.3,
                   0.5, 0.5, 0.6, 0.3), nrow = 2, byrow = TRUE)
n_person <- c(1000, 3000)
w <- genearte_data(lambda[1, ], n_person[1])
fit <- compute_coef(cor(w), 2, 2)
n_rep <- 1000
s <- seq(from = -3, to = 3, by = 0.1)
sigma_res <- diag(1 - diag(lambda[1, ] %*% t(lambda[1, ])))
res <- MASS::mvrnorm(n_rep, mu = rep(0, length(lambda[1, ])), Sigma = sigma_res)
w <- outer(s, lambda[1, ], "*") + res
