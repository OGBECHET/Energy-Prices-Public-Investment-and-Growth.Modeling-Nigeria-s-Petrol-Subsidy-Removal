###############################################################################
# Petrol subsidy removal in Nigeria: where the money goes matters
# Ramsey-Cass-Koopmans model with energy, public capital and a government
# Author: Theresa Ogbeche
#
# Extends rck_nigeria_energy.R. Key change: before 2023 firms paid a subsidized
# pump price while the economy paid the full import cost. Removing the subsidy
# raises what firms pay AND frees government money. The question is what the
# government does with that money:
#   Scenario 1 "Leak":      savings lost (waste, external debt service)
#   Scenario 2 "Transfers": savings paid to households as lump-sum transfers
#   Scenario 3 "Invest":    savings spent on public capital (roads, power)
###############################################################################

rm(list = ls())
dir.create("figures", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)

# ---------------------------------------------------------------------------
# 1. MODEL (per worker, labour normalized to 1)
# ---------------------------------------------------------------------------
# Gross output: Y = A k^alpha kg^eta E^gamma,  kg = public capital
# Firms buy petrol E at the price they face, pf (pump price).
# Their choice of E gives value added (household income):
#   V = B(pf) * k^a * kg^h,  B(p) = (1-gamma) A^(1/(1-gamma)) (gamma/p)^(gamma/(1-gamma))
#   a = capital share of GDP, h = output elasticity of public capital
# Households:  cdot  = (c/theta) * (a*V/k - delta - rho)
#              kdot  = V - c - delta*k + T
# Government:  kgdot = eps * Ig - delta_g * kg     (eps = investment efficiency)
#   Before reform: pays subsidy S0 = (pw - p0)*E0 and baseline investment Ig0.
#   After reform:  S = 0; the freed S0 goes to leak / transfers T / extra Ig.

# ---------------------------------------------------------------------------
# 2. CALIBRATION, with sources
# ---------------------------------------------------------------------------
a       <- 1 - 0.4782  # PWT 10.01 Nigeria 2010-19: 1 - labour share
delta   <- 0.0375      # PWT 10.01 Nigeria 2010-19: depreciation rate
rho     <- 0.0939      # PWT 10.01 Nigeria 2010-19: real return on capital (irr)
theta   <- 2           # Havranek et al. (2015, JIE): mean EIS = 0.5
h       <- 0.128       # Ezzahid & Rafik (2024, Economies 12(5):118): long-run
                       # output elasticity of public capital, panel of 48
                       # African countries incl. Nigeria, 1960-2015 (IMF ICSD
                       # + PWT). Close to Bom & Ligthart's (2014) 0.122.
eps     <- 0.5         # Gupta et al. (2014): about half of public investment
                       # spending becomes capital in low/middle-income countries.
                       # Note: the baseline level cancels out of every % result;
                       # what matters is eps_new below.
eps_new <- 1           # efficiency of the freed money relative to existing
                       # public spending (1 = same; 0.5 = half; section 7)
delta_g <- delta       # ASSUMPTION: public capital depreciates like private
                       # capital (sensitivity in section 7)
sg      <- 0.030 / 1.361  # IMF Nigeria Article IV (2022): public investment
                          # ~3% of old-base GDP. NBS: rebased 2022 GDP is 36.1%
                          # higher, so ~2.2% of rebased GDP (same base as below)

# Full-year 2022 petrol figures (all NNPC / NBS):
#   subsidy bill 2022:      N4.39 trillion (NNPC)
#   supply:                 68 million litres/day (NNPC, Jan-Aug 2022 imports)
#   official pump price:    N165/litre (NNPC: N462 cost - N297 subsidy, Q2 2022)
#   nominal GDP 2022:       N274.23 trillion (NBS rebased)
litres_2022   <- 68e6 * 365
sub_per_litre <- 4.39e12 / litres_2022            # ~N177/litre, full-year avg
p_pump        <- 165
p0 <- 1
pw <- (p_pump + sub_per_litre) / p_pump           # full cost / pump price ~2.07
subsidy_share <- 4.39 / 274.23                    # 1.60% of GDP
gamma <- subsidy_share / (pw / p0 - 1)            # pump-price petrol share ~1.5%
# Check: gamma * GDP / pump price reproduces the 68m litres/day exactly.

A <- 1   # units normalization (cancels out of all reported % changes)

B <- function(p) (1 - gamma) * A^(1 / (1 - gamma)) * (gamma / p)^(gamma / (1 - gamma))
V_of <- function(k, kg, p) B(p) * k^a * kg^h

# ---------------------------------------------------------------------------
# 3. STEADY STATES
# ---------------------------------------------------------------------------
# For a given level of public investment Ig, kg = eps*Ig/delta_g, and k solves
# a*V/k = rho + delta.
k_given_kg <- function(kg, p) (a * B(p) * kg^h / (rho + delta))^(1 / (1 - a))

# Pre-reform: Ig0 = sg * V0 (fixed point, since V depends on kg)
kg0 <- 1
for (i in 1:500) {
  k0  <- k_given_kg(kg0, p0)
  V0  <- V_of(k0, kg0, p0)
  kg0 <- eps * sg * V0 / delta_g
}
k0  <- k_given_kg(kg0, p0); V0 <- V_of(k0, kg0, p0)
Ig0 <- sg * V0
c0  <- V0 - delta * k0
Y0  <- V0 / (1 - gamma)
Phi <- (pw - p0) * gamma * Y0 / p0   # freed subsidy money, held fixed in naira terms

cat(sprintf("gamma = %.4f, pw/p0 = %.3f, subsidy = %.2f%% of GDP\n\n",
            gamma, pw, 100 * Phi / V0))

scenarios <- list(
  Leak      = c(g = 0, T = 0),
  Transfers = c(g = 0, T = 1),
  Invest    = c(g = 1, T = 0)
)

post_ss <- function(sc, en = eps_new, dg = delta_g) {
  Ig  <- Ig0 + en * sc["g"] * Phi      # effective spending
  kg  <- eps * Ig / dg
  k   <- k_given_kg(kg, pw)
  V   <- V_of(k, kg, pw)
  T   <- sc["T"] * Phi
  list(k = unname(k), kg = unname(kg), c = unname(V - delta * k + T),
       V = unname(V), T = unname(T), Ig = unname(Ig))
}

# ---------------------------------------------------------------------------
# 4. DYNAMICS: linearize the 3-variable system (c jumps; k, kg predetermined)
# ---------------------------------------------------------------------------
rhs <- function(x, T, Ig, dg = delta_g) {
  c <- x[1]; k <- x[2]; kg <- x[3]
  V <- V_of(k, kg, pw)
  c(c / theta * (a * V / k - delta - rho),
    V - c - delta * k + T,
    eps * Ig - dg * kg)
}

solve_path <- function(sc, en = eps_new, dg = delta_g, Tmax = 400, dt = 0.05) {
  ss <- post_ss(sc, en, dg)
  xs <- c(ss$c, ss$k, ss$kg)
  J  <- matrix(0, 3, 3)                     # numerical Jacobian
  for (j in 1:3) {
    hstep <- 1e-6 * max(1, abs(xs[j]))
    up <- xs; dn <- xs; up[j] <- up[j] + hstep; dn[j] <- dn[j] - hstep
    J[, j] <- (rhs(up, ss$T, ss$Ig, dg) - rhs(dn, ss$T, ss$Ig, dg)) / (2 * hstep)
  }
  e   <- eigen(J)
  st  <- which(Re(e$values) < 0)
  if (length(st) != 2) stop("Not a saddle: need exactly 2 stable roots")
  Vs  <- e$vectors[, st]
  # k0 and kg0 are given from the pre-reform steady state
  kg_start <- eps * Ig0 / dg
  bcoef <- solve(Vs[2:3, ], c(k0 - ss$k, kg_start - ss$kg))
  tt <- seq(0, Tmax, by = dt)
  X  <- sapply(tt, function(t) Re(Vs %*% (bcoef * exp(e$values[st] * t)))) + xs
  list(t = tt, c = X[1, ], k = X[2, ], kg = X[3, ], ss = ss,
       eig = e$values, J = J, V = V_of(X[2, ], X[3, ], pw))
}

# Welfare: consumption-equivalent gain vs staying in the pre-reform steady
# state forever (what % of c0 would households pay / need to be indifferent)
welfare_gain <- function(path) {
  u  <- path$c^(1 - theta) / (1 - theta)
  W1 <- sum(exp(-rho * path$t) * u) * (path$t[2] - path$t[1])
  W0 <- c0^(1 - theta) / ((1 - theta) * rho)
  100 * ((W1 / W0)^(1 / (1 - theta)) - 1)
}

# ---------------------------------------------------------------------------
# 5. RESULTS
# ---------------------------------------------------------------------------
res <- lapply(scenarios, solve_path)
tab <- data.frame(
  scenario        = names(res),
  k_longrun_pct   = sapply(res, function(r) 100 * (r$ss$k  / k0  - 1)),
  kg_longrun_pct  = sapply(res, function(r) 100 * (r$ss$kg / kg0 - 1)),
  V_longrun_pct   = sapply(res, function(r) 100 * (r$ss$V  / V0  - 1)),
  c_impact_pct    = sapply(res, function(r) 100 * (r$c[1]  / c0  - 1)),
  c_longrun_pct   = sapply(res, function(r) 100 * (r$ss$c  / c0  - 1)),
  # private-capital root (the other stable root is just -delta_g, from kg)
  half_life_yrs   = sapply(res, function(r) { st <- Re(r$eig[Re(r$eig) < 0])
                      log(2) / -st[which.max(abs(st + delta_g))] }),
  welfare_pct     = sapply(res, welfare_gain)
)
tab[, -1] <- round(tab[, -1], 2)
cat("Main results (% change from pre-reform steady state):\n")
print(tab, row.names = FALSE)
write.csv(tab, "results/fiscal_scenarios.csv", row.names = FALSE)

# Keep the baseline pre-reform values (section 7 overwrites c0, k0, V0)
base <- list(c0 = c0, k0 = k0, V0 = V0, kg0 = kg0)

# ---------------------------------------------------------------------------
# 6. FIGURES
# ---------------------------------------------------------------------------
cols <- c(Leak = "firebrick", Transfers = "darkorange", Invest = "steelblue")
show <- res[[1]]$t <= 80
png("figures/fig4_scenarios.png", width = 1500, height = 550, res = 150)
par(mfrow = c(1, 3), mar = c(4, 4.2, 2, 1))
for (v in c("c", "k", "V")) {
  ref <- switch(v, c = c0, k = k0, V = V0)
  ys  <- lapply(res, function(r) 100 * (r[[v]][show] / ref - 1))
  plot(NA, xlim = c(0, 80), ylim = range(unlist(ys), 0),
       xlab = "Years after reform", ylab = "% change",
       main = switch(v, c = "Consumption", k = "Private capital", V = "Output (GDP)"))
  abline(h = 0, col = "grey70")
  for (s in names(res)) lines(res[[s]]$t[show], ys[[s]], col = cols[s], lwd = 2)
  if (v == "c") legend("topleft", names(res), col = cols, lwd = 2, bty = "n")
  if (v != "c") mtext("Leak and Transfers overlap", side = 3, line = 0.1, cex = 0.6)
}
dev.off()

# ---------------------------------------------------------------------------
# 7. ROBUSTNESS: investment efficiency and public-capital depreciation
# ---------------------------------------------------------------------------
rob <- expand.grid(eps_new = c(0.5, 1.0), delta_g = c(0.025, 0.0375, 0.05))
rob$welfare_Invest <- NA; rob$welfare_Transfers <- NA; rob$welfare_Leak <- NA
dg_save <- delta_g
for (i in seq_len(nrow(rob))) {
  en <- rob$eps_new[i]; delta_g <- rob$delta_g[i]
  # recompute the pre-reform steady state under this delta_g
  kg0 <- 1
  for (j in 1:500) { k0 <- k_given_kg(kg0, p0); V0 <- V_of(k0, kg0, p0); kg0 <- eps * sg * V0 / delta_g }
  k0 <- k_given_kg(kg0, p0); V0 <- V_of(k0, kg0, p0); Ig0 <- sg * V0
  c0 <- V0 - delta * k0; Y0 <- V0 / (1 - gamma); Phi <- (pw - p0) * gamma * Y0 / p0
  rob$welfare_Invest[i]    <- round(welfare_gain(solve_path(scenarios$Invest,    en, delta_g)), 2)
  rob$welfare_Transfers[i] <- round(welfare_gain(solve_path(scenarios$Transfers, en, delta_g)), 2)
  rob$welfare_Leak[i]      <- round(welfare_gain(solve_path(scenarios$Leak,      en, delta_g)), 2)
}
delta_g <- dg_save
cat("\nRobustness (welfare gain, % of consumption):\n"); print(rob, row.names = FALSE)
write.csv(rob, "results/fiscal_robustness.csv", row.names = FALSE)
