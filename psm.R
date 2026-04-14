# Purpose: Data Analysis project - analyze data in R using propensity score matching for A Crash Course in Causality: Inferring Causal Effects from Observational Data

# Propensity score matching is used in observational studies to reduce confounding when comparing a treated 
# group to a control group.
# Each subject has a propensity score: the probability of receiving treatment given observed covariates.
# You then match treated and untreated subjects with similar propensity scores.
# After matching, the distribution of covariates should be balanced between groups, mimicking a randomized experiment (on observed variables).
# Key assumption: No unmeasured confoudning (strong ignorability)
# SMD < 0.1 (good balance)


# load packages ----
library(tidyverse)
library(tableone)
library(MatchIt)

# load data ----
load(url("http://www.ats.ucla.edu/stat/data/psm.RData"))

View(psm)

# optional step - create new dataset with only variables of interest, convert character to numeric

# create table 1 ----
# Describes the overall sample before matching
vars <- c("age", "educ", "black", "hispan", "married", "nodegree", "re74", "re75")
table1 <- CreateTableOne(vars = vars, data = psm, test = FALSE)
print(table1)
# estimate propensity scores ----
# the goal is prediction of treatment, not causal interpretation of coefficients
# Include variables related to treatment and/or outcome
psm$pscore <- glm(treat ~ age + educ + black + hispan + married + nodegree + re74 + re75, 
                  data = psm, family = binomial)$fitted.values

# Plot of propensity score, pre-matching
# check overlap, we want treated and control groups to have overlapping PS distributions
ggplot(psm, aes(x = pscore, fill = factor(treat))) +
  geom_density(alpha = 0.5) +
  labs(title = "Propensity Score Distribution Before Matching", x = "Propensity Score", fill = "Treatment") +
  theme_minimal()

# perform matching ----
# re-estimates the propensity score internally and matches treated to control subjects using nearest-neighbor matching
# matches on the logit of the propensity score (common best practice)
matchit_out <- matchit(treat ~ age + educ + black + hispan + married + nodegree + re74 + re75, 
                       data = psm, method = "nearest", distance = "logit")

summary(matchit_out)
# number of treated and matched controls, balance statistics before and after matching and SMDs for each covariate
plot(matchit_out, type = "jitter")
plot(matchit_out, type = "qq", interactive = FALSE)
plot(matchit_out, type = "hist", interactive = FALSE)

# propensity score plot 
matched_data <- match.data(matchit_out)
# Keeps only matched observations and adds weights and subclass info
ggplot(matched_data, aes(x = pscore, fill = factor(treat))) +
  geom_density(alpha = 0.5) +
  labs(title = "Propensity Score Distribution After Matching", x = "Propensity Score", fill = "Treatment") +
  theme_minimal()


# Matching without caliper ----
# match on logit(propensity score) without a caliper 
# do greedy matching on logit(ps)
psmatch <- Match(Tr = psm$treat, X = log(psm$pscore / (1 - psm$pscore)), M = 1, replace = FALSE)
summary(psmatch)
matched <- psm[unlist(psmatch[c("index.treated", "index.control")]), ]
xvars <- c("age", "female", "meanbp1","aps")
matchedtab1 <- CreateTableOne(vars = xvars, strata = "treat", data = matched, test = FALSE)
print(matchedtab1, smd = TRUE)


# PS matching with caliper ----
# re-do matching using a caliper 
# restricts matches to be within a maximum distance, typically defined as 0.2*SD of the logit(PS)
# prevents poor matches, improve balance, and may discard some treated subjects
# better internal validity and smaller matched sample
psmatch <- Match(Tr = psm$treat, M=1, X = logit(psm$pscore), replace = FALSE, caliper = .2)
matched <- psm[unlist(psmatch[c("index.treated", "index.control")]), ]
xvars <- c("age", "female", "meanbp1","aps")
matchedtab1 <- CreateTableOne(vars = xvars, strata = "treat", data = matched, test = FALSE)
print(matchedtab1, smd = TRUE)



