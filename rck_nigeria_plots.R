###############################################################################
# Draw the graphs on screen (RStudio "Plots" pane)
# Keep this file in the same folder as rck_nigeria_fiscal.R
# Use the left/right arrows in the Plots pane to flip between the graphs.
###############################################################################

# 1. Find and run the model --------------------------------------------------
if (!file.exists("rck_nigeria_fiscal.R")) {
  message("Can't find rck_nigeria_fiscal.R in this folder. ",
          "A window will open: pick rck_nigeria_fiscal.R")
  f <- file.choose()
  setwd(dirname(f))
}
source("rck_nigeria_fiscal.R")

if (!exists("base")) {
  stop("You are using an old copy of rck_nigeria_fiscal.R. ",
       "Download the new one from the chat and replace the old file.")
}

# 2. Settings -----------------------------------------------------------------
cols  <- c(Leak = "firebrick", Transfers = "darkorange", Invest = "steelblue")
years <- res$Leak$t <= 80                # show the first 80 years
pct   <- function(x, ref) 100 * (x / ref - 1)
par(mar = c(4, 4, 2, 1), cex.main = 0.9) # small margins so it fits the pane

draw_one <- function(var, ref, title) {
  ys <- lapply(res, function(r) pct(r[[var]][years], ref))
  plot(NA, xlim = c(0, 80), ylim = range(unlist(ys), 0),
       xlab = "Years after subsidy removal", ylab = "% change", main = title)
  abline(h = 0, col = "grey60", lty = 2)
  for (s in names(res)) lines(res[[s]]$t[years], ys[[s]], col = cols[s], lwd = 3)
  legend("topleft", legend = names(res), col = cols, lwd = 3, bty = "n", cex = 0.8)
}

# 3. Graphs -------------------------------------------------------------------
draw_one("c", base$c0, "Consumption")
draw_one("k", base$k0, "Private capital")
draw_one("V", base$V0, "Output (GDP)")

w <- tab$welfare_pct
names(w) <- tab$scenario
bp <- barplot(w, col = cols[names(w)], ylim = range(c(w, 0)) * 1.3,
              ylab = "Welfare gain (% of consumption)", main = "Welfare by scenario")
abline(h = 0)
text(bp, w, labels = sprintf("%.2f%%", w), pos = ifelse(w > 0, 3, 1))
