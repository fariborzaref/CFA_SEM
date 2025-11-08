# Structural Equation Modeling (SEM) — Dr. Fariborz Aref
# Dataset: 27 OECD countries, 2002–2022
# Goal: Estimate a latent Social Inequality factor and compare structures across
#       the Great Recession (2008–2009) and the Pandemic (2020–2022).

# 0) Packages
required <- c("lavaan","semTools","tidyverse","data.table","psych","semPlot","janitor")
to_install <- setdiff(required, rownames(installed.packages()))
if (length(to_install)) install.packages(to_install, repos = "https://cloud.r-project.org")
invisible(lapply(required, library, character.only = TRUE))

# 1) Data
dat <- data.table::fread("oecd_inequality_2002_2022.csv") |>
  janitor::clean_names()

# expected columns:
# country, year, gini_income, health_ineq, labor_ineq, gdp_pc, unemployment, openness
dat <- dat |>
  dplyr::mutate(
    crisis_period = dplyr::case_when(
      year %in% 2008:2009 ~ "Great_Recession",
      year %in% 2020:2022 ~ "Pandemic",
      TRUE ~ "Normal"
    ),
    crisis_period = factor(crisis_period, levels = c("Great_Recession","Pandemic")),
    log_gdp  = log(pmax(gdp_pc, 1)),
    log_open = log(pmax(openness, 0) + 1)
  ) |>
  dplyr::filter(!is.na(crisis_period))

# 2) Model
model_ineq <- '
  # Measurement
  SocialIneq =~ gini_income + health_ineq + labor_ineq

  # Structural
  SocialIneq ~ log_gdp + unemployment + log_open

  # Covariances (exogenous)
  log_gdp ~~ unemployment + log_open
  unemployment ~~ log_open
'

# 3) Multi-group SEM
fit_multi <- lavaan::sem(
  model_ineq,
  data      = dat,
  group     = "crisis_period",
  estimator = "MLR",
  missing   = "fiml",
  std.lv    = TRUE
)

print(summary(fit_multi, fit.measures = TRUE, standardized = TRUE))

# 4) Measurement invariance (configural → metric → scalar)
inv <- semTools::measurementInvariance(
  model = model_ineq,
  data  = dat,
  group = "crisis_period",
  estimator = "MLR"
)
print(inv)

# 5) Fit indices (compact)
fit_idx <- lavaan::fitMeasures(
  fit_multi, c("cfi","tli","rmsea","srmr","aic","bic","chisq","df","pvalue")
)
print(round(fit_idx, 3))

# 6) Standardized solutions (loadings and paths only)
std_sol <- semTools::standardizedSolution(fit_multi) |>
  dplyr::filter(op %in% c("=~","~")) |>
  dplyr::arrange(group, op, lhs, rhs)

print(std_sol, n = nrow(std_sol))

# 7) Cross-crisis comparison table for key structural paths
param_diff <- lavaan::parameterEstimates(fit_multi, standardized = TRUE) |>
  dplyr::filter(op == "~", lhs == "SocialIneq", rhs %in% c("log_gdp","unemployment","log_open")) |>
  dplyr::select(group, rhs, est, se, z, pvalue, std.all) |>
  dplyr::arrange(group, rhs)

cat("\nCross-crisis structural paths to SocialIneq (standardized shown in std.all):\n")
print(param_diff)

# 8) Visualization (minimal, academic sizing; no heavy lines)
# semPlot uses base graphics. Keep labels readable but not oversized.
png("CFA_SEM/fig_sem_paths.png", width = 1600, height = 1000, res = 180)
semPlot::semPaths(
  object = fit_multi,
  what = "std",
  layout = "tree",
  whatLabels = "std",
  title = TRUE,
  nCharNodes = 0,
  sizeMan = 5.5,   # modest
  sizeLat = 6.5,   # modest
  edge.label.cex = 0.75,
  mar = c(4,4,4,4),
  curvePivot = TRUE,
  groups = "crisis_period",
  pastel = FALSE,  # keep it neutral
  residuals = FALSE
)
dev.off()

# 9) Notes printed to console (concise)
cat("\nInterpretive notes:\n")
cat("• Loadings reflect each dimension’s contribution to SocialIneq.\n")
cat("• Compare standardized paths across crises: GDP, unemployment, openness.\n")
cat("• Expect stronger labor and health contributions in the Pandemic window.\n")

# 10) Save tidy outputs
if (!dir.exists("CFA_SEM")) dir.create("CFA_SEM", recursive = TRUE)
data.table::fwrite(as.data.frame(std_sol),  "CFA_SEM/sem_standardized_loadings_paths.csv")
data.table::fwrite(as.data.frame(param_diff),"CFA_SEM/sem_structural_paths_by_crisis.csv")
data.table::fwrite(
  data.frame(metric = names(fit_idx), value = as.numeric(fit_idx)),
  "CFA_SEM/sem_fit_indices.csv"
)
saveRDS(
  list(fit_multi = fit_multi, std_solution = std_sol, fit_indices = fit_idx, param_diff = param_diff),
  file = "CFA_SEM/sem_oecd_crisis_2008_2022.rds"
)

cat("\n✓ Completed. Tables in CFA_SEM/, figure saved to CFA_SEM/fig_sem_paths.png\n")

