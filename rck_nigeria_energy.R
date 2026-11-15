###############################################################################
# Ramsey-Cass-Koopmans model with an energy input, calibrated to Nigeria
# Experiment: the May 2023 petrol subsidy removal as a permanent, unanticipated
#             rise in the relative price of energy
# Author: Theresa Ogbeche
#
# This mirrors the Mathematica RCK.nb from class (old steady state -> new
# steady state -> linearize -> eigenvalues -> saddle path -> time paths),
# rewritten in R and with Nigerian numbers instead of textbook ones.
###############################################################################

rm(list = ls())
dir.create("figures", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)

# ---------------------------------------------------------------------------
# 1. MODEL
# ---------------------------------------------------------------------------
# Gross output per worker:   Y = A * k^alpha * E^gamma
# Energy E is bought at relative price p, so firms pick E to maximize Y - p*E.
# First-order condition:     gamma * Y / E = p
# Plugging that back in gives net output (GDP / value added):
#   f(k) = B(p) * k^a,  with
#   a    = alpha / (1 - gamma)
#   B(p) = (1 - gamma) * A^(1/(1-gamma)) * (gamma/p)^(gamma/(1-gamma))
# So the model collapses to the same y = B k^a form we used in class, and an
# energy price shock works through B. Higher p -> lower B.
#
# Households: u(c) = c^(1-theta)/(1-theta), discount rate rho
#   Euler:   cdot = (c/theta) * (f'(k) - delta - rho)
#   Capital: kdot = f(k) - c - delta*k

# ---------------------------------------------------------------------------
# 2. CALIBRATION (Nigeria)
# ---------------------------------------------------------------------------
# Penn World Table 10.01, Nigeria, 2010-2019 averages:
#   labsh (labour share of GDP) = 0.4782  -> capital share of GDP a = 0.5218
#   delta (depreciation rate)   = 0.0375
#   irr   (real return on capital) = 0.0939 -> used as rho, since in steady
#                                     state f'(k) - delta = rho
# Set use_pwt <- TRUE to recompute these yourself from the pwt10 package.
use_pwt <- FALSE
if (use_pwt) {
  if (!requireNamespace("pwt10", quietly = TRUE)) install.packages("pwt10")
  pw  <- pwt10::pwt10.01
  ng  <- subset(pw, isocode == "NGA" & year >= 2010 & year <= 2019)
  a_share <- 1 - mean(ng$labsh, na.rm = TRUE)
  delta   <- mean(ng$delta, na.rm = TRUE)
  rho     <- mean(ng$irr,   na.rm = TRUE)
} else {
  a_share <- 1 - 0.4782
  delta   <- 0.0375
  rho     <- 0.0939
}

theta <- 2       # Havranek, Horvath, Irsova & Rusnak (2015, J. Int. Econ.
                 # 96:100-118): mean EIS across 2,735 estimates is 0.5, so
                 # theta = 1/0.5 = 2. They find lower EIS in poorer countries,
                 # so theta = 3 is reported as a robustness check.

# Energy = petrol, so the cost share and the price shock refer to the same good.
# gamma = annual petrol spending / nominal GDP, 2024 (first full post-reform year)
#   consumption: 51.8 million litres/day in 2024 (NMDPRA, as reported)
#   price: N1,189.12/litre, December 2024 (NBS PMS Price Watch), when pump
#          prices were fully market-determined
#   nominal GDP 2024: N372.82 trillion (NBS rebased GDP, 2019 base year)
litres_2024 <- 51.8e6 * 365
gdp_2024    <- 372.82e12
gamma <- litres_2024 * 1189.12 / gdp_2024          # = 0.0603
gamma_low <- litres_2024 * (1037.66 / 1.3832) / gdp_2024
# gamma_low uses the June 2024 price (N750.19, backed out of NBS's June 2025
# figure of N1,037.66, +38.32% y/y) = 0.0380, a lower bound.
A     <- 1       # Units normalization, not a parameter: it scales every level
                 # and cancels out of every % change, eigenvalue and half-life.

# Energy price shock (NBS PMS Price Watch, August 2022 vs August 2023):
#   petrol: N189.46 -> N626.70 per litre (+230.78%)
#   headline CPI inflation, Aug 2023 y/y: 25.80%
# Real (relative) petrol price ratio:
p0 <- 1
p1 <- (626.70 / 189.46) / 1.2580
# p0 = 1 is also a units normalization; only the ratio p1/p0 matters.

alpha <- a_share * (1 - gamma)   # implied capital elasticity in gross output

# ---------------------------------------------------------------------------
# 3. STEADY STATES
# ---------------------------------------------------------------------------
B_of_p <- function(p, gamma, A = 1) {
  (1 - gamma) * A^(1 / (1 - gamma)) * (gamma / p)^(gamma / (1 - gamma))
}

steady_state <- function(B, a, delta, rho) {
  k <- (a * B / (rho + delta))^(1 / (1 - a))   # from f'(k) = rho + delta
  c <- B * k^a - delta * k                      # from kdot = 0
  list(k = k, c = c, y = B * k^a)
}

B0 <- B_of_p(p0, gamma, A)
B1 <- B_of_p(p1, gamma, A)
ss0 <- steady_state(B0, a_share, delta, rho)   # before subsidy removal
ss1 <- steady_state(B1, a_share, delta, rho)   # after

cat("Relative energy price ratio p1/p0:", round(p1, 3), "\n")
cat("Petrol cost share gamma:", round(gamma, 4), "(lower bound", round(gamma_low, 4), ")\n\n")
cat("Initial steady state:  k =", round(ss0$k, 4), " c =", round(ss0$c, 4),
    " y =", round(ss0$y, 4), "\n")
cat("New steady state:      k =", round(ss1$k, 4), " c =", round(ss1$c, 4),
    " y =", round(ss1$y, 4), "\n")
cat(sprintf("Long-run change: k %.2f%%, c %.2f%%, y %.2f%%\n\n",
            100 * (ss1$k / ss0$k - 1), 100 * (ss1$c / ss0$c - 1),
            100 * (ss1$y / ss0$y - 1)))

# ---------------------------------------------------------------------------
# 4. LINEARIZE AROUND THE NEW STEADY STATE (same matrix as in class)
# ---------------------------------------------------------------------------
# [cdot]   [a11 a12] [c - c*]
# [kdot] = [a21 a22] [k - k*]
keq <- ss1$k; ceq <- ss1$c
fpp <- B1 * a_share * (a_share - 1) * keq^(a_share - 2)   # f''(k*)
fp  <- B1 * a_share * keq^(a_share - 1)                   # f'(k*)

a11 <- 0
a12 <- (ceq / theta) * fpp
a21 <- -1
a22 <- fp - delta          # equals rho at the steady state
MM  <- matrix(c(a11, a12, a21, a22), nrow = 2, byrow = TRUE)
cat("Jacobian matrix:\n"); print(round(MM, 6))

ev  <- eigen(MM)
mu  <- sort(Re(ev$values))
sse  <- mu[1]   # negative (stable) root
lse1 <- mu[2]   # positive (unstable) root
cat("\nEigenvalues:", round(mu, 6), "\n")
cat("Characteristic polynomial: mu^2 -", round(a22, 4), "mu +",
    round(-a12, 6), "\n")
cat("Saddle point? ", sse < 0 & lse1 > 0, "\n")
cat("Speed of convergence:", round(-sse, 4),
    " -> half-life:", round(log(2) / -sse, 1), "years\n\n")

# ---------------------------------------------------------------------------
# 5. SADDLE PATH AND TIME PATHS
# ---------------------------------------------------------------------------
# k starts at the old steady state (capital can't jump); c jumps onto the
# stable arm. From row 1: mu * (c - c*) = a12 * (k - k*).
keqb <- ss0$k
A1   <- (keqb - keq) / ((sse - a11) / a12)    # same constant as in class
ct <- function(t) ceq + A1 * exp(sse * t)
kt <- function(t) keq + ((sse - a11) / a12) * A1 * exp(sse * t)

cat("c(t) =", round(ceq, 5), "+", round(A1, 5), "* exp(", round(sse, 5), "t )\n")
cat("k(t) =", round(keq, 5), "+", round((sse - a11) / a12 * A1, 5),
    "* exp(", round(sse, 5), "t )\n")
cat(sprintf("Consumption jumps on impact from %.4f to %.4f (%.2f%%)\n\n",
            ss0$c, ct(0), 100 * (ct(0) / ss0$c - 1)))

tt <- seq(0, 60, by = 0.25)
paths <- data.frame(t = tt, k = kt(tt), c = ct(tt), y = B1 * kt(tt)^a_share)
paths$k_pct <- 100 * (paths$k / ss0$k - 1)
paths$c_pct <- 100 * (paths$c / ss0$c - 1)
paths$y_pct <- 100 * (paths$y / ss0$y - 1)
write.csv(paths, "results/rck_nigeria_paths.csv", row.names = FALSE)

# ---------------------------------------------------------------------------
# 6. FIGURES
# ---------------------------------------------------------------------------
png("figures/fig1_time_paths.png", width = 1400, height = 500, res = 150)
par(mfrow = c(1, 3), mar = c(4, 4, 2, 1))
plot(tt, paths$c, type = "l", lwd = 2, xlab = "t (years)", ylab = "c",
     main = "Consumption")
abline(h = c(ss0$c, ceq), lty = c(3, 2), col = "grey40")
plot(tt, paths$k, type = "l", lwd = 2, xlab = "t (years)", ylab = "k",
     main = "Capital")
abline(h = c(ss0$k, keq), lty = c(3, 2), col = "grey40")
plot(tt, paths$y, type = "l", lwd = 2, xlab = "t (years)", ylab = "y",
     main = "Output (GDP)")
abline(h = c(ss0$y, ss1$y), lty = c(3, 2), col = "grey40")
dev.off()

png("figures/fig2_phase_diagram.png", width = 900, height = 750, res = 150)
kgrid <- seq(0.5 * keq, 1.3 * ss0$k, length.out = 300)
plot(kgrid, B1 * kgrid^a_share - delta * kgrid, type = "l", lwd = 2,
     xlab = "k", ylab = "c", main = "Phase diagram: subsidy removal",
     ylim = range(c(0.8 * ceq, 1.1 * ss0$c)))
lines(kgrid, B0 * kgrid^a_share - delta * kgrid, lty = 3)
abline(v = keq, lwd = 2); abline(v = ss0$k, lty = 3)
lines(paths$k, paths$c, col = "firebrick", lwd = 2)
points(c(ss0$k, keq, ss0$k), c(ss0$c, ceq, ct(0)), pch = 19,
       col = c("grey40", "black", "firebrick"))
legend("topleft", bty = "n", cex = 0.8,
       legend = c("kdot = 0 (new)", "kdot = 0 (old)", "cdot = 0 (new)",
                  "saddle path"),
       lty = c(1, 3, 1, 1), col = c(1, 1, 1, "firebrick"))
dev.off()

png("figures/fig3_percent_deviations.png", width = 1000, height = 650, res = 150)
plot(tt, paths$k_pct, type = "l", lwd = 2, ylim = range(paths[, 5:7]),
     xlab = "t (years)", ylab = "% change from pre-reform steady state",
     main = "Transition after subsidy removal")
lines(tt, paths$c_pct, lwd = 2, lty = 2)
lines(tt, paths$y_pct, lwd = 2, lty = 4)
abline(h = 0, col = "grey60")
legend("bottomright", bty = "n", legend = c("k", "c", "y"),
       lty = c(1, 2, 4), lwd = 2)
dev.off()

# ---------------------------------------------------------------------------
# 7. ROBUSTNESS: June vs December 2024 petrol price for gamma; EIS 0.5 vs 0.33
# ---------------------------------------------------------------------------
sens <- expand.grid(gamma = round(c(gamma_low, gamma), 4), theta = c(2, 3))
sens$k_change_pct <- NA; sens$c_impact_pct <- NA; sens$half_life <- NA
for (i in seq_len(nrow(sens))) {
  g  <- sens$gamma[i]; th <- sens$theta[i]
  b0 <- B_of_p(p0, g); b1 <- B_of_p(p1, g)
  s0 <- steady_state(b0, a_share, delta, rho)
  s1 <- steady_state(b1, a_share, delta, rho)
  m  <- matrix(c(0, (s1$c / th) * b1 * a_share * (a_share - 1) * s1$k^(a_share - 2),
                 -1, rho), 2, byrow = TRUE)
  neg <- min(Re(eigen(m)$values))
  c0  <- s1$c + (m[1, 2] / neg) * (s0$k - s1$k)   # c on saddle path at k0
  sens$k_change_pct[i] <- round(100 * (s1$k / s0$k - 1), 2)
  sens$c_impact_pct[i] <- round(100 * (c0 / s0$c - 1), 2)
  sens$half_life[i]    <- round(log(2) / -neg, 1)
}
cat("Sensitivity analysis:\n"); print(sens)
write.csv(sens, "results/rck_nigeria_sensitivity.csv", row.names = FALSE)
