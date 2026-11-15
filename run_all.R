###############################################################################
# Run everything: reproduces every table and figure in the proposal.
# Set the working directory to this folder first, then source this file.
###############################################################################
source("rck_nigeria_energy.R")            # benchmark: class RCK model, price shock only
source("rck_nigeria_fiscal.R")            # main model: energy + public capital + government
source("rck_nigeria_proposal_figures.R")  # figures used in the proposal PDF
cat("\nDone. Figures are in figures/ and tables are in results/.\n")
