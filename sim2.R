source("blup.R")
genearte_data <- function(lambda, s) {
  ############# generate data #########
  n_person <- length(s)
  sigma_res <- diag(1 - diag(lambda %*% t(lambda)))
  res <- MASS::mvrnorm(n_person, mu = rep(0, length(lambda)), Sigma = sigma_res)
  w <- outer(s, lambda, "*") + res
  return(w)
}
##############################################

lambda <- matrix(c(0.8, 0.8, 0.6, 0.3,
                   0.5, 0.5, 0.6, 0.3), nrow = 2, byrow = TRUE)
n_person <- c(500, 3000)
n_rep <- 1000

mat_bias_blup <- array(NA, dim = c(2, n_person[2], ncol = n_rep))
mat_bias_blp <- array(NA, dim = c(2, n_person[2], ncol = n_rep))

s <- rnorm(n_person[2])
for (i in seq_len(2)) {
  for (j in seq_len(n_rep)) {
    w <- genearte_data(lambda[i, ], s)
    fit <- compute_coef(cor(w), 2, 2)
    mat_bias_blup[i, , j] <- w %*% fit$blup_coef - s
    mat_bias_blp[i, , j] <- w %*% fit$blp_coef - s
  }
}

mean_bias_blup <- apply(mat_bias_blup, c(1, 2), mean)
mean_bias_blp <- apply(mat_bias_blp, c(1, 2), mean)

plot(mean_bias_blup[1, ] ~ s, ylab = "bias", ylim = c(-2, 2))
points(mean_bias_blp[1, ] ~ s, col = "red")