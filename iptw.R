
# Inverse Probability of Treatment Weighting (IPTW) is a method used 
# in observational studies to estimate causal effects by creating 
# a pseudo-population in which the treatment assignment is independent of measured confounders. 
# The idea is to weight each individual's contribution to the analysis by the 
# inverse of the probability of receiving the treatment they actually received, given their covariates.
# After weighting, covariates should be balanced between treated and untreated groups, mimicking randomization,
# allowing for unbiased estimation of treatment effects.
# Core assumption: Consistency (treatment well defined), exchangeability (no unmeasured confounding), positivity (non-zero probability of treatment for all covariate patterns).
# correct PS model specification is crucial for IPTW to yield unbiased estimates.
# Estimands: Average Treatment Effect (ATE) - the average effect of treatment in the entire population; Average Treatment Effect on the Treated (ATT) - the average effect of treatment among those who received the treatment.


# Step by step IPTW in R 
# Esimate the propensity score
ps_model <- glm(
  treat ~ age + sex + education + baseline_score,
  data = df,
  family = binomial()
)

df$ps <- predict(ps_model, type = "response")

# compute stabilized IPTW weights
p_treat <- mean(df$treat == 1)
df$w_iptw <- ifelse(df$treat == 1, p_treat / df$ps, (1 - p_treat) / (1 - df$ps))

# check positivity and extreme weights
summary(df$w_iptw)
quantile(df$w_iptw, probs = c(0.01, 0.99))
         
# optional trimming
df$w_iptw_trimming <- pmin(df$w_iptw, quantile(df$w_iptw, 0.99))

# Assess covariate balance 
library(cobalt)
bal.tab(
  treat ~ age + sex + education + baseline_score,
  data = df,
  weights = df$w_iptw,
  method = "weighting"
)
# Target: standardized mean differences (SMD) < 0.1 for all covariates

# Estimate treatment effect using weighted outcome model
# Continuous outcome
library(survey)
design <- svydesign(ids = ~1, data = df, weights = ~w_iptw)
outcome_model <- svyglm(outcome ~ treat, design = design)

# binary outcome
binary_outcome_model <- svyglm(outcome ~ treat, design = design, family = binomial())
# The coefficient on treat estimates the ATE
# Robust (sandwich) standard errors account for weighting

# common pitfalls 
# Extreme weights (Fix: Use stabilized weigths, trim)
# Pool balance (Fix: Improve PS model, e.g., add nonlinear terms, interactions)
# Violating positivity (Fix: Restrict sample)
# Using unweighted SEs (Fix: Use survey package)
# Matching + IPTW (Fix: Avoid unless justified)

# When NOT to use IPTW
# Very small samples
# Severe lack of overlap
# Strong unmeasured confounding

# Alternatives
# Propensity score matching
# Overlap weights
# Doubly robust estimators (IPTW + outcome model)
# Targeted Maximum likelihood (TMLE)





