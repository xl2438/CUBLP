source("blup.R")
dat <- read.csv("LR.csv")
head(dat)
dat <- dat[, -1]
dim(dat)
dat <- na.omit(dat)
dim(dat)
dat[, c(1, 2)] <- (dat[, c(1, 2)] - colMeans(dat[, c(1, 2)])) /
apply(dat[, c(1, 2)], 2, sd)
r <- cor(dat, use = "complete.obs")
lambda_hat <- rep(NA, 6)
lambda_hat[c(1, 2)] <- sqrt(r[1, 2])
for (i in 3:6) {
  lambda_hat[i] <- 0.5 * sum(r[i, c(1, 2)]) / lambda_hat[1]
}
temp <- compute_coef(r, 2, 4, 1)
my_tab <- cbind(temp$blup_coef, temp$blp_coef, temp$lambda_hat)
colnames(my_tab) <- c("CUBLP", "BLP", "lambda_hat")
stargazer::stargazer(my_tab)
lambda_hat <- temp$lambda_hat
names(lambda_hat) <- c("rating 1", "rating 2", "accuracy", "pronunciation",
 "fluency", "rhythm")
stargazer::stargazer(lambda_hat)

s_hat <- as.matrix(dat) %*% blup_coef
s_hat_blp <- as.matrix(dat) %*% blp_coef
s_hat_alt <- apply(dat[, c(1, 2)], 1, mean)
ind_low <- which(s_hat_blp < -0.5)
ind_high <- which(s_hat_blp > 0.5)
png("lower_S.png")
plot(s_hat[ind_low] ~ s_hat_blp[ind_low], xlab = "BLP", ylab = "CUBLP")
abline(a = 0, b = 1)
dev.off()
png("higher_S.png")
plot(s_hat[ind_high] ~ s_hat_blp[ind_high], xlab = "BLP", ylab = "CUBLP")
abline(a = 0, b = 1)
dev.off()
#####################################
t(temp$lambda_hat) %*% solve(r) %*% temp$lambda_hat
########### using only features ########
sigma_x <- r[-c(1, 2), -c(1, 2)]
lambda_x_hat <- lambda_hat[-c(1, 2)]
sigma_x_res <- sigma_x - lambda_x_hat %*% t(lambda_x_hat)
blup_x_coef <- solve(sigma_x_res) %*% lambda_x_hat / (t(lambda_x_hat) %*% solve
  (sigma_x_res) %*% lambda_x_hat)[1]
blp_x_coef <- solve(sigma_x_res) %*% lambda_x_hat / (1 + (t(lambda_x_hat) %*%
 solve(sigma_x_res) %*% lambda_x_hat)[1])
blup_x_coef
s_hat_x <- as.matrix(dat[, -c(1, 2)]) %*% blup_x_coef
# plot((s_hat - s_hat_alt) ~ s_hat_alt, ylim = c(-4, 4))
# abline(h = 0)
# points((s_hat_x - s_hat_alt) ~ s_hat_alt, col = "red")
# abline(h = 0)
cor(s_hat_alt, s_hat, method = "spearman")
cor(s_hat_alt, s_hat_x, method = "spearman")
