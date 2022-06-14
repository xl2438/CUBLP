s <- seq(from = -3, to = 3, by = 0.01)
png("conditional_bias.png", width = 5, height = 5, units = "in", res = 400)
plot(s * (0.2 - 1) ~ s, type = "l", xlab = "S", ylab = paste
  ("conditional", "bias"))
lines(s * (0.5 - 1) ~ s, lty = 2)
lines(s * (0.8 - 1) ~ s, lty = 3)
legend("topright", legend = c(paste0("ratio", "=", "0.2"), paste0("ratio", "=",
  "0.5"), paste0("ratio", "=", "0.8")), lty = c(1, 2, 3))
dev.off()