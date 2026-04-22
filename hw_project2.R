


# Carry out an IPTW causal analysis ----
# For this assignment we will use data from Lalonde(1986), that aimed to evaluate the 
# impact of National Supported Work (NSW) Demonstration, which is a labor training program,
# on post-intervention income levels. Interest is in estimating the causal effect of this training
# program on income.

# load packages ----
library(tidyverse)
library(tableone)
library(Matching)
library(ipw)
library(survey)
library(MatchIt)
library(sandwich) 

# load the lalonde data which is in the MatchIt package ----
data(lalonde)
# The data have n=614 subjects, and 10 variables
# age in years, educ years of schooling, black indicator variable for blacks
# hispan indicator variable for Hispanics
# married indicator variable for married, nodegree indicator variable for high school degree
# re74 real earnings in 1974, re75 real earnings in 1975, re78 real earnings in 1978
# treat indicator variable for participation in the NSW program (1 = participated, 0 = did not participate)
# The outcome is re78 (real earnings in 1978) and the treatment is treat (participation in the NSW program).
# The potential confounding variables are age, educ, black, hispan, married, nodegree, re74, and re75.

# Fit a propensity score model. Use a logistic regression model, where the outcome is treatment. Include 
# the 8 confounding variables in the model as predictors, with no interaction terms or non-linear terms (such as squared terms).
# Obtain the propensity score for each subject. Next, obtain the inverse probability of treatment weights for each subject.
lalonde$black <- as.numeric(lalonde$race == "black")
lalonde$hispan <- as.numeric(lalonde$race == "hispan")
lalonde$pscore <- glm(treat ~ age + educ + black + hispan + married + nodegree + re74 + re75, 
                      data = lalonde, family = binomial)$fitted.values

# create weigths 
weight <- ifelse(lalonde$treat==1, 1/lalonde$pscore, 1/(1-lalonde$pscore))
summary(weight)  # Q1

# compute stabilized IPTW weights
# p_treat <- mean(lalonde$treat == 1)
# lalonde$w_iptw <- ifelse(lalonde$treat == 1, p_treat / lalonde$pscore, (1 - p_treat) / (1 - lalonde$pscore))

# apply weights to data
weighteddata <- svydesign(ids = ~1, data = lalonde, weights = ~weight)

# weighted table 1
weightedtable <- svyCreateTableOne(vars = c("age", "educ", "black", "hispan", "married", "nodegree", "re74", "re75"), 
                          strata = "treat", data = weighteddata, test = FALSE)

print(weightedtable, smd = TRUE) # Q2


# Using IPTW find the estimate and 95% confidence interval for the average causal effect. 
# This can be obtained from svyglm
outcome_model <- svyglm(re78 ~ treat, design = weighteddata)
betaiptw <- coef(outcome_model)["treat"]  # Q3) point estimate 224.6763

SE <- sqrt(vcov(outcome_model)["treat", "treat"])
# SE2 <- sqrt(diag(vcovHC(outcome_model, type = "HC0")))  # robust SE ?? two large??
CI_lower <- betaiptw - 1.96 * SE
CI_upper <- betaiptw + 1.96 * SE
CI <- c(CI_lower, CI_upper)  # Q4) 95% CI (-1559.353, 2008.706)


# Now truncate the weights at the 1st and 99th percentiles. This can be done with the trunc = 0.01 option in svyglm (does not work)
# using IPTW with the truncated weights, find the estimate and 95% confidence interval for the average causal effect.
truncated_outcome_model <- svyglm(re78 ~ treat, design = weighteddata, trunc = 0.01)
beta_truncated <- coef(truncated_outcome_model)["treat"]  # Q5) point estimate  224.6763
SE_truncated <- sqrt(vcov(truncated_outcome_model)["treat", "treat"])
CI_truncated_lower <- beta_truncated - 1.96 * SE_truncated
CI_truncated_upper <- beta_truncated + 1.96 * SE_truncated
CI_truncated <- c(CI_truncated_lower, CI_truncated_upper)  # Q6) 95% CI (-1559.353, 2008.706)


# use ipw package instead ----
ipw_out <- ipwpoint(exposure = treat, family = "binomial", link = "logit", 
                    denominator = ~ age + educ + black + hispan + married + nodegree + re74 + re75, 
                    data = lalonde, trunc = 0.01)

# summary of trunc weight 
summary(ipw_out$weights.trunc)
