compute_coef <- function(r_w, length_y, length_x, var_s) {
  ######### estimate loadings ##############
  lambda_hat_y <- rep(NA, length_y)
  lambda_hat_x <- rep(NA, length_x)
  r_y <- r_w[1:length_y, 1:length_y]
  diag(r_y) <- 0
  lambda_hat_y[1:length_y] <- sqrt(sum(r_y) / (length_y * (length_y - 1)))
  r_xy <- r_w[-c(1:length_y), c(1:length_y)]
  r_xy <- matrix(r_xy, nrow = length_x, ncol = length_y)
  lambda_hat_x[1:length_x] <- apply(r_xy, 1, sum) / (length_y * lambda_hat_y[1])
  lambda_hat <- c(lambda_hat_y, lambda_hat_x)
  ########## estimate error covariance matrix###
  r_res <- r_w - lambda_hat %*% t(lambda_hat)
  blp_coef <- coef <- solve(r_res) %*% lambda_hat / (1/var_s + (t(lambda_hat) %*%
 solve(r_res) %*% lambda_hat)[1]) #BLP coefficients
  blup_coef <- coef <- solve(r_res) %*% lambda_hat / (t(lambda_hat) %*% solve
      (r_res) %*% lambda_hat)[1] #BLUP
  return(list(blup_coef = blup_coef, blp_coef = blp_coef, lambda_hat =
    lambda_hat))
}
compute_coef_uneq_loadings <- function(r_w, length_y, length_x, var_s) {
  lambda_hat_y <- rep(NA, length_y)
  lambda_hat_x <- rep(NA, length_x)
  r_y <- r_w[1:length_y, 1:length_y]
  diag(r_y) <- 0
  r_xy <- r_w[-c(1:length_y), c(1:length_y)]
  r_xy <- matrix(r_xy, nrow = length_x, ncol = length_y)
  lambda_hat_y <- rnorm(length(lambda_hat_y), 0, 0.3)
  lambda_hat_x <- rnorm(length(lambda_hat_x), 0, 0.3)
  lambda_hat_last <- c(lambda_hat_y, lambda_hat_x)
  repeat {
    for (i in seq_len(length(lambda_hat_y))) {
      lambda_hat_y[i] <- (lambda_hat_y[-i] %*% r_y[i, -i] + lambda_hat_x %*%
       r_xy[, i]) / (sum(lambda_hat_y[-i]^2) + sum(lambda_hat_x^2))
    }
    lambda_hat_y <- abs(lambda_hat_y)
    lambda_hat_x <- as.vector(r_xy %*% lambda_hat_y / sum(lambda_hat_y^2))
    lambda_hat <- c(lambda_hat_y, lambda_hat_x)
    if (sum((lambda_hat - lambda_hat_last)^2) < 1e-20) {
      break
    }
    lambda_hat_last <- lambda_hat
  }
  r_res <- r_w - lambda_hat %*% t(lambda_hat)
  blp_coef <- solve(r_res) %*% lambda_hat / (1/var_s + (t(lambda_hat) %*% solve
      (r_res) %*% lambda_hat)[1]) #BLP coefficients
  blup_coef <- solve(r_res) %*% lambda_hat / (t(lambda_hat) %*% solve
      (r_res) %*% lambda_hat)[1] #BLUP
  return(list(blup_coef = blup_coef, blp_coef = blp_coef, lambda_hat =
   lambda_hat))
}