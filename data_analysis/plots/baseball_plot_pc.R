
mc_upper <- min(100*(0.95 + 1.96*sqrt(0.05*0.95/200)),100)
mc_lower <- max(100*(0.95 - 1.96*sqrt(0.05*0.95/200)),0)

pdf("baseball_PC.pdf")
plot(c(1:m_real), post_covS2[,1], 
     type = "n",
     ylim = c(min(c(95, mc_lower, post_covS1[,1:2], post_covS1[,1]))-1, 
              101), xlim = c(0, m_real),
     xaxt = 'n', yaxt = "n",
     xlab = "", ylab = "", cex.axis = 2, cex.lab = 2)
axis(1, at = 1:m_real, labels = baseball$player, cex.axis = 1.1, las = 2)
axis(2, at = seq(80,100, by = 5), labels = TRUE, cex.axis = 2)
abline(h = mc_lower, col = "green", lwd = 2)
abline(h = mc_upper, col = "green", lwd = 2)
polygon(
  x = c(-3, m_real+3, m_real+3, -3),
  y = c(mc_lower, mc_lower, mc_upper, mc_upper),
  col = adjustcolor("green", alpha.f = 0.1),
  border = NA
)
abline(h = 95, col = "black", lwd = 2)
points(c(1:m_real), post_covS2[,1], pch = 19, col = "darkgreen")
points(c(1:m_real), post_covS1[,1], pch = 19, col = "blue")
points(c(1:m_real), post_covS1[,2], pch = 19, col = "red")
legend(
  x = "bottomleft",
  pch = 19,
  col = c( "darkgreen", "blue", "red", "green"),
  legend = c(
    expression(I[i]^YL * " (M3)"),
    expression(I[i]^N * " (M4)"),
    expression(I[i]^YL * " (M4)"),
    "MC error"
  ), 
  cex = 1.5
)
dev.off()


# post_covS2 <- post_cov
# post_covS1 <- post_cov
