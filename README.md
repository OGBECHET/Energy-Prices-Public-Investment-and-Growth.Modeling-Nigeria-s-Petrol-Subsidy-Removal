# Petrol Subsidy Removal in Nigeria: Does It Matter Where the Money Goes?

A Ramsey-Cass-Koopmans growth model, calibrated to Nigerian data, that asks whether the 2023 petrol subsidy removal leaves Nigerians better off, and how much the answer depends on what the government does with the money it saves.

**Author:** Theresa Ogbeche, PhD candidate in Econometrics and Quantitative Economics, Southern Illinois University Carbondale

## The question

Nigeria ended its petrol subsidy in May 2023. In 2022 the subsidy cost ₦4.39 trillion, about 1.6% of GDP, and pump prices more than tripled once it was gone. Most of the debate has been about the price shock. This project looks at the other half of the reform: the freed fiscal space.

I compare three uses of the savings:

| Scenario | What happens to the savings |
| --- | --- |
| **Leak** | Lost to waste or external debt service |
| **Transfers** | Returned to households as cash |
| **Invest** | Spent on public capital such as roads and power |

## Main results

| Scenario | Long-run GDP | Consumption on impact | Welfare gain |
| --- | ---: | ---: | ---: |
| Leak | −2.28% | −0.78% | −1.51% |
| Transfers | −2.28% | +1.14% | +0.39% |
| Invest | +13.29% | +1.06% | +2.03% |

Welfare gain is the permanent change in consumption that leaves households equally well off.

The reform pays off only if the money is used. Investing it beats cash transfers only when the new spending is about as efficient as existing public spending. If only half of it becomes real capital, transfers do better.

![Transition paths](figures/fig_paths.png)

![Welfare by scenario](figures/fig_welfare.png)

## The model

- **Firms** use private capital, public capital and petrol. Removing the subsidy raises the pump price firms pay to the full cost of petrol.
- **Households** save and consume as in the standard Ramsey model with CRRA utility.
- **Government** builds public capital from public investment, and only part of each naira spent becomes capital.

The model is solved by linearizing around the new steady state and finding the saddle path from the eigenvalues of the Jacobian, following the same steps as the classroom Mathematica exercise.

## Calibration

| Parameter | Value | Source |
| --- | ---: | --- |
| Capital share of GDP | 0.522 | Penn World Table 10.01, Nigeria, 2010–2019 |
| Depreciation rate | 0.0375 | Penn World Table 10.01, Nigeria, 2010–2019 |
| Discount rate | 0.094 | Penn World Table 10.01, real return on capital |
| Full cost / pump price of petrol | 2.07 | NNPC subsidy bill and supply, 2022 |
| Petrol cost share | 0.015 | NNPC and NBS rebased GDP, 2022 |
| Public investment / GDP | 0.022 | IMF Article IV, converted to rebased GDP |
| Output elasticity of public capital | 0.128 | Ezzahid & Rafik (2024), 48 African countries |
| Inverse elasticity of intertemporal substitution | 2 | Havranek et al. (2015) |
| Public investment efficiency | 0.5 | Gupta et al. (2014) |
| Public capital depreciation | 0.0375 | Assumed; 0.025–0.05 tested |

## Files

| File | What it does |
| --- | --- |
| `run_all.R` | Runs everything below in order |
| `rck_nigeria_energy.R` | Benchmark: the classroom RCK model with petrol, price shock only |
| `rck_nigeria_fiscal.R` | Main model with public capital and the three scenarios |
| `rck_nigeria_proposal_figures.R` | Figures used in the proposal |
| `rck_nigeria_plots.R` | Draws the main graphs on screen |
| `figures/` | All saved figures |
| `results/` | Results tables as CSV files |

## How to run it

1. Install [R](https://cran.r-project.org) (and optionally [RStudio](https://posit.co/download/rstudio-desktop/)). No extra packages are needed.
2. Download this repository (green **Code** button, then **Download ZIP**) and unzip it.
3. Set the working directory to the unzipped folder. In RStudio: **Session → Set Working Directory → Choose Directory**. In R.app on a Mac: **Misc → Change Working Directory**.
4. Run:

```r
source("run_all.R")
```

## References

- Bom, P. R. D., & Ligthart, J. E. (2014). What have we learned from three decades of research on the productivity of public capital? *Journal of Economic Surveys*, 28(5), 889–916.
- Ezzahid, E., & Rafik, H. (2024). Analyzing the impact of public capital on private capital productivity in a panel of African nations. *Economies*, 12(5), 118.
- Feenstra, R. C., Inklaar, R., & Timmer, M. P. (2015). The next generation of the Penn World Table. *American Economic Review*, 105(10), 3150–3182.
- Gupta, S., Kangur, A., Papageorgiou, C., & Wane, A. (2014). Efficiency-adjusted public capital and growth. *World Development*, 57, 164–178.
- Havranek, T., Horvath, R., Irsova, Z., & Rusnak, M. (2015). Cross-country heterogeneity in intertemporal substitution. *Journal of International Economics*, 96(1), 100–118.

## Status

Work in progress: term paper for Macroeconomic Theory III. Planned extensions include CES production, a comparison with Nigeria's actual 2023–2025 data, a nonlinear solution, and an open-economy version.
