# Structural Equation Modeling for OECD Inequality

**Author:** Dr. Fariborz Aref  
**Discipline:** Quantitative Sociology and Inequality Research  
**License:** MIT  

### Purpose
Estimate a latent **Social Inequality** construct and compare measurement and structural relations across two crisis windows. Uses multi group SEM with robust estimation and full information maximum likelihood for missing data.

### Structure
├── sem_run.R


### Methods
- Latent factor **SocialIneq** measured by gini_income, health_ineq, labor_ineq  
- Structural regressions on `log_gdp`, `unemployment`, `log_open`  
- Multi group comparison for **Great_Recession** and **Pandemic**  
- Measurement invariance sequence: configural, metric, scalar  
- Robust estimator **MLR**, missing handled via **FIML**

### Quick example
```r
library(lavaan)
model <- '
  SocialIneq =~ gini_income + health_ineq + labor_ineq
  SocialIneq ~ log_gdp + unemployment + log_open
  log_gdp ~~ unemployment + log_open
  unemployment ~~ log_open
'
fit <- sem(model, data = dat, group = "crisis_period", estimator = "MLR", missing = "fiml", std.lv = TRUE)
summary(fit, fit.measures = TRUE, standardized = TRUE)
