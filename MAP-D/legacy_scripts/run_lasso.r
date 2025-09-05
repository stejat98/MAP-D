## lasso script


library(glmnet)
library(dplyr)
library(tidyr)
library(purrr)
library(caret)

# Define outcomes and covariates
outcomes <- c("BMI", "TRIG_HDL_RATIO", "HDL", "LDL", "systolic_BP", "diastolic_BP", "HbA1c")
baseline_vars <- adjustments
proteomic_vars <- validated_proteins_unique_all_step1_step_2

set.seed(123)

# Create stratified train/test split
split_df <- data %>%
  group_by(GlycemicStatus) %>%
  mutate(row_id = row_number()) %>%
  mutate(n = n()) %>%
  mutate(train_flag = row_id <= floor(0.7 * n)) %>%
  ungroup()

train_df <- split_df %>% filter(train_flag) %>% select(-row_id, -n, -train_flag)
test_df  <- split_df %>% filter(!train_flag) %>% select(-row_id, -n, -train_flag)

# LASSO (no bootstrap)
run_lasso <- function(y_var, predictors, train_df, test_df) {
  x_train <- as.matrix(train_df[, predictors])
  y_train <- train_df[[y_var]]
  x_test  <- as.matrix(test_df[, predictors])
  y_test  <- test_df[[y_var]]
  
  complete_train <- complete.cases(x_train, y_train)
  x_train <- x_train[complete_train, ]
  y_train <- y_train[complete_train]
  
  complete_test <- complete.cases(x_test, y_test)
  x_test <- x_test[complete_test, ]
  y_test <- y_test[complete_test]
  
  if (length(unique(y_train)) < 2 || nrow(x_train) < 10) return(NA)
  
  cv_fit <- tryCatch(cv.glmnet(x_train, y_train, alpha = 1, standardize = TRUE),
                     error = function(e) return(NULL))
  if (is.null(cv_fit)) return(NA)
  
  pred <- predict(cv_fit, s = cv_fit$lambda.min, newx = x_test)
  if (any(is.na(pred)) || any(is.na(y_test)) || length(pred) != length(y_test)) return(NA)
  
  cor(y_test, pred)^2
}

# Compute R2
results <- expand_grid(
  Outcome = outcomes,
  ModelType = c("BaselineOnly", "ProteomicOnly", "BaselinePlusProteomic")
) %>%
  mutate(
    Predictors = case_when(
      ModelType == "BaselineOnly" ~ list(baseline_vars),
      ModelType == "ProteomicOnly" ~ list(proteomic_vars),
      ModelType == "BaselinePlusProteomic" ~ list(c(baseline_vars, proteomic_vars))
    )
  ) %>%
  mutate(
    R2 = map2_dbl(Outcome, Predictors, ~ run_lasso(.x, .y, train_df, test_df))
  )

# Ensure model type ordering
results <- results %>%
  mutate(ModelType = factor(ModelType, levels = c("BaselineOnly", "ProteomicOnly", "BaselinePlusProteomic")))


saveRDS(results, "/n/groups/patel/sivateja/UKB/PEWAS_results/lasso_results_R2.RDS")


library(ggplot2)
library(dplyr)


saveRDS(results, "/n/groups/patel/sivateja/UKB/PEWAS_results/lasso_results_R2_bootstrap_CIs.RDS")

results <- readRDS("/n/groups/patel/sivateja/UKB/PEWAS_results/lasso_results_R2.RDS")



pdf("/n/groups/patel/sivateja/UKB/PEWAS_results/forest_plot_R2.pdf")

# Order outcomes by R2 from the BaselinePlusProteomic model
ordered_outcomes <- results %>%
  filter(ModelType == "BaselinePlusProteomic") %>%
  arrange(R2) %>%
  pull(Outcome)

# Update factor levels
results <- results %>%
  mutate(Outcome = factor(Outcome, levels = ordered_outcomes))

# Forest plot with CIs
ggplot(results, aes(x = R2, y = Outcome, color = ModelType)) +
  geom_point(position = position_dodge(width = 0.6), size = 3) +
  geom_errorbarh(aes(xmin = CI_lower, xmax = CI_upper), height = 0.2,
                 position = position_dodge(width = 0.6)) +
  scale_x_continuous(limits = c(0, 1)) +
  labs(
    x = "R2 (Test Set with 95% CI)",
    y = "Clinical Risk Factor",
    color = "Model Type"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "right",
    axis.title.y = element_blank()
  )

dev.off()



# Order outcomes by R2 from the BaselinePlusProteomic model
ordered_outcomes <- results %>%
  filter(ModelType == "BaselinePlusProteomic") %>%
  arrange(R2) %>%
  pull(Outcome)

# Update factor levels
results <- results %>%
  mutate(Outcome = factor(Outcome, levels = ordered_outcomes))

# Forest plot with CIs
ggplot(results, aes(x = R2, y = Outcome, color = ModelType)) +
  geom_point(position = position_dodge(width = 0.6), size = 3) +
  scale_x_continuous(limits = c(0, 1)) +
  labs(
    x = "R2 (Test)",
    y = "Clinical Risk Factor",
    color = "Model Type"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "right",
    axis.title.y = element_blank()
  )








