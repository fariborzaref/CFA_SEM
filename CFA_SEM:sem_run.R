############################################################
# STRUCTURAL EQUATION MODELING (SEM)
# Author: Dr. Fariborz Aref
# Data: 27 OECD countries, 2002 to 2022
# Goal: Estimate a latent Social Inequality factor and compare structures
#       across the Great Recession and the Pandemic windows.
############################################################

# 0) Packages
required <- c("lavaan","semTools","tidyverse","data.table","psych","semPlot","janitor")
to_install <- setdiff(required, rownames(installed.packages()))
if (length(to_install)) install.packages(to_install, repos = "https://cloud.r-project.org")
suppressPackageStartupMessages(invisible(lapply(required, library, character.only = TRUE)))

set.seed(2025)
options(stringsAsFactors = FALSE)

# 1) Data
csv_path <- "oecd_inequality_2002_2022.csv"   # place next to sem_run.R or set a full path
if (!file.exists(csv_path)) stop("File not found: ", csv_path)

dat <- data.table::fread(csv_path) |>
  janitor::clean_names()

need <- c("country","year","gini_income","health_ineq","labor_ineq","gdp_pc","unemployment","openness")
miss <- setdiff(need, names(dat))
if (length(miss)) stop("Missing columns: ", paste(miss, collapse = ", "))

dat <- dat |>
  dplyr::mutate(
    crisis_period = dplyr::case_when(
      year %in% 2008:2009 ~ "Great_Recession",
      year %in% 2020:2022 ~ "Pandemic",
      TRUE ~ NA_character_
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

  # Covariances among exogenous predictors
  log_gdp ~~ unemployment + log_open
  unemployment ~~ log_open
'

# 3) Multi-group SEM by crisis period
fit_multi <- lavaan::sem(
  model     = model_ineq,
  data      = dat,
  group     = "crisis_period",
  estimator = "MLR",
  missing   = "fiml",
  std.lv    = TRUE
)

cat("\n=== MULTI GROUP SEM SUMMARY (standardized) ===\n")
print(summary(fit_multi, fit.measures = TRUE, standardized = TRUE))

# 4) Measurement invariance
inv <- semTools::measurementInvariance(
  model     = model_ineq,
  data      = dat,
  group     = "crisis_period",
  estimator = "MLR"
)
cat("\n=== MEASUREMENT INVARIANCE (configural, metric, scalar) ===\n")
print(inv)

# 5) Fit indices
fit_idx <- lavaan::fitMeasures(
  fit_multi, c("cfi","tli","rmsea","srmr","aic","bic","chisq","df","pvalue")
)
cat("\n=== FIT INDICES ===\n")
print(round(fit_idx, 3))

# 6) Standardized loadings and paths
std_sol <- semTools::standardizedSolution(fit_multi) |>
  dplyr::filter(op %in% c("=~","~")) |>
  dplyr::arrange(group, op, lhs, rhs)

cat("\n=== STANDARDIZED SOLUTION (loadings and paths) ===\n")
print(std_sol, n = nrow(std_sol))

# 7) Cross crisis structural paths to SocialIneq
param_diff <- lavaan::parameterEstimates(fit_multi, standardized = TRUE) |>
  dplyr::filter(op == "~", lhs == "SocialIneq", rhs %in% c("log_gdp","unemployment","log_open")) |>
  dplyr::select(group, rhs, est, se, z, pvalue, std.all) |>
  dplyr::arrange(group, rhs)

cat("\n=== CROSS CRISIS STRUCTURAL PATHS TO SocialIneq (std.all) ===\n")
print(param_diff)

# 8) Path diagram
if (!dir.exists("CFA_SEM/figs")) dir.create("CFA_SEM/figs", recursive = TRUE)
png("CFA_SEM/figs/sem_paths.png", width = 1600, height = 1000, res = 180)
semPlot::semPaths(
  object         = fit_multi,
  what           = "std",
  layout         = "tree",
  whatLabels     = "std",
  title          = TRUE,
  nCharNodes     = 0,
  sizeMan        = 5.5,
  sizeLat        = 6.5,
  edge.label.cex = 0.75,
  mar            = c(4,4,4,4),
  curvePivot     = TRUE,
  groups         = "crisis_period",
  pastel         = FALSE,
  residuals      = FALSE
)
dev.off()

# 9) Save tidy outputs
if (!dir.exists("CFA_SEM/out")) dir.create("CFA_SEM/out", recursive = TRUE)

data.table::fwrite(as.data.frame(std_sol),   "CFA_SEM/out/sem_standardized_loadings_paths.csv")
data.table::fwrite(as.data.frame(param_diff),"CFA_SEM/out/sem_structural_paths_by_crisis.csv")
data.table::fwrite(
  data.frame(metric = names(fit_idx), value = as.numeric(fit_idx)),
  "CFA_SEM/out/sem_fit_indices.csv"
)

saveRDS(
  list(fit_multi = fit_multi, std_solution = std_sol, fit_indices = fit_idx, param_diff = param_diff),
  file = "CFA_SEM/out/sem_oecd_crisis_2008_2022.rds"
)

# 10) Console notes
cat("\nInterpretive notes:\n")
cat("• Loadings show contributions of gini, health, and labor to SocialIneq.\n")
cat("• Compare standardized paths across crises for GDP, unemployment, and openness.\n")
cat("• Expect stronger labor and health links during the Pandemic window.\n")
cat("\n✓ Completed. Tables in CFA_SEM/out and figure in CFA_SEM/figs/sem_paths.png\n")
