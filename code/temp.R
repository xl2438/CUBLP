source("blup.R")
genearte_data <- function(lambda, n_person) {
  ############# generate data #########
  s <- rnorm(n_person)
  sigma_res <- diag(1 - diag(lambda %*% t(lambda)))
  res <- MASS::mvrnorm(length(s), mu = rep(0, length(lambda)), Sigma =
   sigma_res)
  w <- outer(s, lambda, "*") + res
  return(list(w, s))
}
##############################################

lambda <- matrix(c(0.7, 0.7, 0.4,
                   0.5, 0.5, 0.5), nrow = 2, byrow = TRUE)
temp <- genearte_data(lambda[1, ], 3000)
w <- temp[[1]]
s <- temp[[2]]
fit <- compute_coef(cor(w), 2, 1)
w_hat_blp <- w %*% fit$blp_coef
w_hat_blup <- w %*% fit$blup_coef
png("prediction.png")
plot(w_hat_blp ~ s, xlim = c(-4, 4), ylim = c(-4, 4), xlab = "true score", ylab
  = "predicted score")
points(w_hat_blup ~ s, col = 2)
abline(a = 0, b = 1)
legend("bottomright", legend = c("BLP", "CUBLP"), pch = 1, col =  c(1, 2))
dev.off()

s <- rep(seq(from = -3, to = 3, by = 0.01), 100)
s <- sort(s)
lambda <- c(0.7, 0.7, 0.4)
sigma_res <- diag(1 - diag(lambda %*% t(lambda)))
res <- MASS::mvrnorm(length(s), mu = rep(0, length(lambda)), Sigma =
 sigma_res)
w <- outer(s, lambda, "*") + res
s_hat_blp <- w %*% fit$blp_coef
s_hat_blup <- w %*% fit$blup_coef
s_mat_blp <- matrix(s_hat_blp, nrow = 100)
s_mat_blup <- matrix(s_hat_blup, nrow = 100)
bias_blp <- apply(s_mat_blp, 2, mean) - seq(from = -3, to = 3, by = 0.01)
bias_blup <- apply(s_mat_blup, 2, mean) - seq(from = -3, to = 3, by = 0.01)
s_temp <- seq(from = -3, to = 3, by = 0.01)
png("prediction_bias.png")
plot(bias_blp ~ s_temp, xlab = "true score", ylab = "prediction bias")
points(bias_blup ~ s_temp, col = 2)
legend("topright", legend = c("BLP", "CUBLP"), pch = 1, col =  c(1, 2))
dev.off()