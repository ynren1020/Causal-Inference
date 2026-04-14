

# PSM data analysis project - hw1 



library(tableone)
library(Matching)
library(MatchIt)
library(ggplot2)

# load the lalonde data which is in the MatchIt package ----
data(lalonde)

# Find the standardized differences for all of the confounding variables (pre-matching).
# What is the standardized difference for married (to nearest hundredth) # 0.72
table1 <- CreateTableOne(vars = c("age", "educ", "race", "married", "nodegree", "re74", "re75"), 
                         strata = "treat", data = lalonde, test = FALSE)
print(table1, smd = TRUE)


# What is the raw (unadjusted) mean of real earnings in 1978 for treated subjects
# minus the mean of real earnings in 1978 for untreated subjects?  # -635
mean(lalonde$re78[lalonde$treat == 1]) - mean(lalonde$re78[lalonde$treat == 0])

# Fit a propensity score model. Use a logistic regression model, where the outcome is treatment.
# Include the 8 confounding variables in the model as predictors, with no interaction terms or 
# non-linear terms (such as squared terms). Obtain the propensity score for each subject from the fitted model.
lalonde$pscore <- glm(treat ~ age + educ + race + married + nodegree + re74 + re75, 
                  data = lalonde, family = binomial)$fitted.values
summary(lalonde$pscore)  # min 0.00908 max 0.85315

# Plot of propensity score, pre-matching
ggplot(lalonde, aes(x = pscore, fill = factor(treat))) +
  geom_density(alpha = 0.5) +
  labs(title = "Propensity Score Distribution Before Matching", x = "Propensity Score", fill = "Treatment") +
  theme_minimal()


# Now carry out propensity score matching using the Match function. 
# Before using the Match function, first do:
set.seed(931139)
# Setting the seed will ensure that you end up with a matched data set that is the same as the one used to create 
# the solutions. Use options to specify pair matching, without replacement, no caliper.
# Match on the propensity score itself, not logit of the propensity score. 
# Obtain the standardized differences for the matched data.
# Question 4 - What is the standardized difference for married (to nearest hundredth) after matching? # 0.027
match_out <- Match(Y = lalonde$re78, Tr = lalonde$treat, X = lalonde$pscore, M = 1, replace = FALSE)
matched_data <- lalonde[c(match_out$index.treated, match_out$index.control), ]
table2 <- CreateTableOne(vars = c("age", "educ", "race", "married", "nodegree", "re74", "re75"), 
                         strata = "treat", data = matched_data, test = FALSE)
print(table2, smd = TRUE)


# For the propensity score matched data: Question 5 which variable has the largest standardized difference? nodegree


# Re-do the matching, but use a caliper this time, set the caliper = 0.1 in the options in the Match function.
# Again, before running the Match function, set the seed;
set.seed(931139)
# What is the standardized difference for married (to nearest hundredth) after matching with a caliper? # 0.002
match_out_caliper <- Match(Y = lalonde$re78, Tr = lalonde$treat, X = lalonde$pscore, M = 1, replace = FALSE, caliper = 0.1)
matched_data_caliper <- lalonde[c(match_out_caliper$index.treated, match_out_caliper$index.control), ]

# Question 6 How many matched pairs are there?
length(match_out_caliper$index.treated)  # 111 pairs


# Use the matched data set (from propensity score matching with a caliper=0.1) to carry out the outcome analysis
# Question 7: For the matched data, what is the mean of real earnings in 1978 for treated subjects minus the 
# mean of real earnings in 1978 for untreated subjects?   # 1246.81
mean(matched_data_caliper$re78[matched_data_caliper$treat == 1]) - mean(matched_data_caliper$re78[matched_data_caliper$treat == 0])  # 1794.5

# carry out a paired t-test for the effect of treatment on earnings. What are the values of the 95% confidence interval?
t.test(matched_data_caliper$re78[matched_data_caliper$treat == 1], matched_data_caliper$re78[matched_data_caliper$treat == 0], paired = TRUE)

# best practice 
library(dplyr)
paired-df <- matched_data_caliper |>
  arrange(pair_id, treat) |> # ensure consistent ordering
  group_by(pair_id) |>
  summarise(
    re78_treat = re78[treat ==1],
    re78_ctrl  = re78[treat ==0],
    .groups = "drop"
  )
t.test(paired_df$re78_treat, paired_df$re78_ctrl, paired = TRUE)








