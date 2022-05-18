# -----------------------------------------------------------------------------
# In this script, I do some exploratory data analysis on the plasma samples
# from the Biocrates data. This will follow the same sequence as with the 
# lavage samples. 
# Sarah Samorodnitsky
# -----------------------------------------------------------------------------

# Setting the working directory
currentwd <- "/Users/sarahsamorodnitsky/Dropbox/HIV-COPD/Plasma_Biocrates/"

# Loading in the cleaned metabolite and patient data (loads in both)
load(paste0(currentwd, "BiocratesPlasmaProcessed.rda"))

# -----------------------------------------------------------------------------
# Creating histograms of the scaled and centered data and also the log-
# transformed versions. 
# -----------------------------------------------------------------------------

# Creating a new dataset to store the scaled data in 
plasma_scaled <- plasma_processed[,3:ncol(plasma_processed)]
plasma_log_scaled <- plasma_processed[,3:ncol(plasma_processed)]

# Scaling all the rows because we want each metabolite to be centered, not each subject. 
plasma_scaled <- t(apply(plasma_scaled, 1, function(x) scale(x, center = TRUE, scale = TRUE)))
plasma_log_scaled <- t(apply(plasma_log_scaled, 1, function(x) scale(log(x+1), center = TRUE, scale = TRUE)))

# Adding the column names back
colnames(plasma_scaled) <- colnames(plasma_processed)[-c(1:2)]
colnames(plasma_log_scaled) <- colnames(plasma_processed)[-c(1:2)]

# -----------------------------------------------------------------------------
# Creating a PCA plot of the metabolite data.
# -----------------------------------------------------------------------------

# First, running an SVD on the metabolite data. 
# Using the logged and scaled data

# Transpose the matrix so that the rows are the patients
# Then the UD part of the SVD corresponds to patient-level variation. 
plasma_svd <- svd(t(plasma_log_scaled))

# -----------------------------------------------------------------------------
# Univariate testing: using paired t-tests to compare metabolite levels
# between cases and controls in each matched pair.
# Using a false-discovery rate correction after repeated testing. 
# -----------------------------------------------------------------------------

# Creating a dataframe to store the results
case_control_t.test_res <- data.frame(test.stat = numeric(nrow(plasma_log_scaled)),
                                      p.value = numeric(nrow(plasma_log_scaled)))

for (i in 1:nrow(case_control_t.test_res)) { # iterating through each metabolite
  # extracting the case and control ids
  case_ids <- subject_processed_plasma$id[subject_processed_plasma$ccstat == 1]
  control_ids <- subject_processed_plasma$id[subject_processed_plasma$ccstat == 2]
  
  # subsetting the metabolite i levels for cases and controls
  case_metab_i <- plasma_log_scaled[i, which(colnames(plasma_log_scaled) %in% case_ids)]
  control_metab_i <- plasma_log_scaled[i, which(colnames(plasma_log_scaled) %in% control_ids)]
  
  # performing the paired t-test
  res_i <- t.test(case_metab_i, control_metab_i, paired = TRUE)
  
  # storing the results
  case_control_t.test_res$test.stat[i] <- res_i$statistic
  case_control_t.test_res$p.value[i] <- res_i$p.value
}

# Performing a false discovery rate correction
p.values <- case_control_t.test_res$p.value
names(p.values) <- plasma_processed$Metabolite

case_control_t.test_fdr_adjust <- p.adjust(p.values, method = "fdr")
names(case_control_t.test_fdr_adjust) <- plasma_processed$Metabolite


# -----------------------------------------------------------------------------
# Performing DWD on the metabolite data to compare cases and controls. 
# This was the standard DWD that I initially started with. 
# -----------------------------------------------------------------------------

library("SparseM")
library("DWDLargeR")

# Creating a vector of case-control labels, either -1 or 1
case_control_labels <- ifelse(subject_processed_plasma$ccstat == 1, 1, -1) # change cases to 1, controls to -1

# Changing the form of the data (not sure why Adam did this)
plasma_log_scaled.DWD <- SparseM::`[.matrix.csr`(plasma_log_scaled)

# # Adjusting the penalty parameter
penalty.plasma <- DWDLargeR::penaltyParameter(X = plasma_log_scaled.DWD,
                                              y = case_control_labels,
                                              expon = 1)

# Cross-validated DWD:
# Iterating through each case-control pair
# ind_of_pairs <- seq(from = 1, to = length(case_control_labels), by = 2)
# cv.scores <- lapply(ind_of_pairs, function(i) {
#   # For each case-control pair, i
#   w.vec = genDWD(X = plasma_log_scaled.DWD[,-c(i, i+1)],
#                  y = case_control_labels[-c(i, i+1)],
#                  expon = 1,
#                  C = penalty.plasma,
#                  scaleFea = FALSE)$w
#   cv.scores = w.vec%*%plasma_log_scaled.DWD[,c(i, i+1)]
#   cv.scores
# })
# 
# # Doing a t-test on the scores (why is this not paired?)
# # Putting the scores in a matrix where the first column corresponds to the score
# # for the case and the second for the control. Makes it easier to feed into the
# # t.test function.
# cv.scores.mat <- matrix(unlist(cv.scores), ncol = 2, byrow = T)
# colnames(cv.scores.mat) <- c("case", "control")
# case_control_dwd_t.test <- t.test(x = cv.scores.mat[,1], y = cv.scores.mat[,2],
#                                   paired = TRUE,
#                                   alternative = "greater")
# 
# # Saving the results
# results_to_save <- list(cv.scores = cv.scores, t.test.res = case_control_dwd_t.test)
# save(results_to_save, file = paste0(currentwd, "DWD_On_Original_Data_Results.rda"))

# -----------------------------------------------------------------------------
# Fitting the DWD model on the full training dataset to get the weights of 
# each metabolite
# -----------------------------------------------------------------------------

full_training_res_plasma <- genDWD(X = plasma_log_scaled.DWD,
                                  y = case_control_labels,
                                  expon = 1,
                                  C = penalty.plasma,
                                  scaleFea = FALSE)

# -----------------------------------------------------------------------------
# Implementing a permutation testing-based approach to DWD.
# -----------------------------------------------------------------------------

library("diproperm")
library("SparseM")
library("DWDLargeR")
library("gtools")

# # Creating a vector of case-control labels, either -1 or 1
# case_control_labels <- ifelse(subject_processed_plasma$ccstat == 1, 1, -1) # change cases to 1, controls to -1
# 
# # Changing the form of the data (not sure why Adam did this)
# plasma_log_scaled.DWD <- SparseM::`[.matrix.csr`(plasma_log_scaled)
# 
# # Adjusting the penalty parameter
# penalty.plasma <- DWDLargeR::penaltyParameter(X = plasma_log_scaled.DWD,
#                                               y = case_control_labels,
#                                               expon = 1)
# 
# # Changing the case-control vector to a list so that it is easier to permute
# ind_of_pairs <- seq(from = 1, to = length(case_control_labels), by = 2)
# case_control_labels.list <- lapply(ind_of_pairs, function(pair) {
#   case_control_labels[c(pair, pair+1)]
# })
# 
# # Changing the form of the data (not sure why Adam did this)
# plasma_log_scaled.DWD <- SparseM::`[.matrix.csr`(plasma_log_scaled)
# 
# # Permutation testing approach to DWD:
# # For 500 replications:
# #     Randomly permute who is a case and who is a control within each pair
# #     For each pair:
# #         Do DWD, leaving one pair out each time.
# 
# # Number of permutations to do
# npermute <- 500
# 
# # Test statistics
# DWD_t.stats_permute_plasma <- c()
# 
# for (permute_iter in 1:npermute) {
#   # Randomly permute which subject in each pair is a case and which is a control
#   case_control_labels.permute.list <- lapply(case_control_labels.list, function(pair) permute(pair))
#   case_control_labels.permute <- unlist(case_control_labels.permute.list)
# 
#   # LOOCV DWD on pairs
#   cv.scores.permute <- lapply(ind_of_pairs, function(i) {
#     # For each case-control pair, i
#     w.vec = genDWD(X = plasma_log_scaled.DWD[,-c(i, i+1)],
#                    y = case_control_labels.permute[-c(i, i+1)],
#                    expon = 1,
#                    C = penalty.plasma,
#                    scaleFea = FALSE)$w
#     cv.scores = w.vec%*%plasma_log_scaled.DWD[,c(i, i+1)]
#     cv.scores
#   })
# 
#   # Storing the scores in a matrix to do a paired t-test
#   cv.scores.permute.mat <- matrix(unlist(cv.scores.permute), ncol = 2, byrow = T)
#   colnames(cv.scores.permute.mat) <- c("case", "control")
# 
#   # Doing a t-test on the scores (why is this not paired?)
#   # Unlist the cv.scores so that each case and control entry within a pair is adjacent.
#   case_control_dwd_t.test.permute <- t.test(x = cv.scores.permute.mat[,1],
#                                             y = cv.scores.permute.mat[,2],
#                                             paired = TRUE,
#                                             alternative = "greater")
# 
#    # Storing the results
#    DWD_t.stats_permute_plasma[permute_iter] <- case_control_dwd_t.test.permute$statistic
# }
# 
# save(DWD_t.stats_permute_plasma, file = paste0(currentwd, "Permutation_Testing_DWD_Test_Stats.rda"))
# 


# -----------------------------------------------------------------------------
# Summing across the metabolites for each metabolite family and comparing cases 
# and controls. This gives us a sense of whether the cumulative levels of each
# metabolite differ between cases and controls. 
# The idea here is to sum within each subject, so retain the subjects and the
# metabolite families, and then compare between cases and controls. 
# -----------------------------------------------------------------------------

# First, adding in the metabolite information back
plasma_log_scaled_info <- cbind(plasma_processed[,1:2], plasma_log_scaled)

# Splitting up the data into cases and controls
plasma_log_scaled_cases <- plasma_log_scaled_info[,colnames(plasma_log_scaled_info) %in% subject_processed_plasma[subject_processed_plasma$ccstat == 1,]$id,]
plasma_log_scaled_cases <- cbind(plasma_processed[,1:2], plasma_log_scaled_cases)
plasma_log_scaled_controls <- plasma_log_scaled_info[,colnames(plasma_log_scaled_info) %in% subject_processed_plasma[subject_processed_plasma$ccstat == 2,]$id,]
plasma_log_scaled_controls <- cbind(plasma_processed[,1:2], plasma_log_scaled_controls)

# Which are the families of metabolits?
metabolite_families <- unique(plasma_processed$Class)
n_families <- length(metabolite_families)

# Adding within the families but retaining the subjects
plasma_log_scaled_cases_cumulative <- matrix(nrow = n_families,
                                             ncol = ncol(plasma_log_scaled_cases) - 2)

plasma_log_scaled_controls_cumulative <- matrix(nrow = n_families,
                                                ncol = ncol(plasma_log_scaled_cases) - 2)

for (i in 1:n_families) {
  metabolite_data_cases <- plasma_log_scaled_cases[plasma_log_scaled_cases$Class == metabolite_families[i], -c(1:2)]
  metabolite_data_controls <- plasma_log_scaled_controls[plasma_log_scaled_controls$Class == metabolite_families[i], -c(1:2)]
  
  plasma_log_scaled_cases_cumulative[i,] <- apply(metabolite_data_cases, 2, sum)
  plasma_log_scaled_controls_cumulative[i,] <- apply(metabolite_data_controls, 2, sum)
}

# Doing a paired t-test for each family
metabolite_casecontrol_cumulative <- data.frame(Family = metabolite_families,
                                                Test.Stat = numeric(n_families),
                                                P.value = numeric(n_families))

for (i in 1:n_families) {
  res_metab_i <- t.test(plasma_log_scaled_cases_cumulative[i,], 
                        plasma_log_scaled_controls_cumulative[i,], 
                        paired = TRUE)
  metabolite_casecontrol_cumulative$Test.Stat[i] <- res_metab_i$statistic
  metabolite_casecontrol_cumulative$P.value[i] <- res_metab_i$p.value
}

# -----------------------------------------------------------------------------
# Considering the FEV1_percent_predicted (FEV1pp) column as an outcome.
#    * Computing the correlations between FEV1pp and each metabolite
#    * Running a LASSO model with FEV1pp as the outcome with all metabolites
# -----------------------------------------------------------------------------

# First, checking the order still matches between the clinical and metabolite data --
# (since its been awhile) 
all(colnames(plasma_log_scaled) == subject_processed_plasma$id)

# Computing the correlations between each metabolite (each row in the lavage data)
# with FEV1_percent_predicted
n_metabolite <- nrow(plasma_log_scaled)
corr_fev1pp_metabolites <- data.frame(Metabolite = character(n_metabolite),
                                      Correlation = numeric(n_metabolite),
                                      P.Value = numeric(n_metabolite),
                                      Q.Value = numeric(n_metabolite))

for (i in 1:n_metabolite) {
  # Save the current metabolite
  metab_i <- plasma_processed$Metabolite[i]
  
  # Calculate correlation
  corr_i <- cor.test(plasma_log_scaled[i,],
                     subject_processed_plasma$FEV1_percent_predicted)
  
  # Save results
  corr_fev1pp_metabolites$Metabolite[i] <- metab_i
  corr_fev1pp_metabolites$Correlation[i] <- corr_i$estimate
  corr_fev1pp_metabolites$P.Value[i] <- corr_i$p.value
}

# Applying an FDR correction
corr_fev1pp_metabolites$Q.Value <- p.adjust(corr_fev1pp_metabolites$P.Value, method = "fdr")

# Sorting the rows in order of Q-value
corr_fev1pp_metabolites <- corr_fev1pp_metabolites[order(corr_fev1pp_metabolites$Q.Value),]

# Adding indices for the rows
rownames(corr_fev1pp_metabolites) <- 1:n_metabolite

# How many were below each threshold?
sum(corr_fev1pp_metabolites$Q.Value < 0.05)
sum(corr_fev1pp_metabolites$Q.Value < 0.1)
sum(corr_fev1pp_metabolites$Q.Value < 0.2)

# -----------------------------------------------------------------------------
# Considering a LASSO regression model on the data -- 
# -----------------------------------------------------------------------------

# Loading in the package
library(glmnet)

# Training the model, finding the best lambda value
set.seed(1)

# Save the outcome
outcome <- subject_processed_plasma$FEV1_percent_predicted

# Full training data
X <- t(plasma_log_scaled)
Y <- outcome

# Fit the model
lasso_cv <- cv.glmnet(X, Y, family = "gaussian", alpha = 1)
lambda.min <- lasso_cv$lambda.min

# Fit the model again on the test data with the min lambda
lasso_final <- glmnet(X, Y, family = "gaussian", alpha = 1, lambda = lambda.min)

# Which coefficients are non-zero?
lasso_beta <- matrix(lasso_final$beta)
rownames(lasso_beta) <- plasma_processed$Metabolite
non_zero_betas <- lasso_beta[lasso_beta != 0, ]


# Consider cross validation --
case_control_labels <- ifelse(subject_processed_plasma$ccstat == 1, 1, -1)
ind_of_pairs <- seq(from = 1, to = length(case_control_labels), by = 2)

cv.lasso <- lapply(ind_of_pairs, function(i) {
  # For each case-control pair, i
  # find optimal penalty
  cv.fit <- cv.glmnet(x = X[-c(i, i+1),],
                      y = Y[-c(i, i+1)],
                      family = "gaussian",,
                      alpha = 1)
  
  # save optimal penalty
  lambda.min <- cv.fit$lambda.min
  
  # refit model with penalty
  fit <- glmnet(x = X[-c(i, i+1),],
                y = Y[-c(i, i+1)],
                family = "gaussian", 
                alpha = 1,
                lambda = lambda.min)
  # calculate scores on test pair
  pred.out <- matrix(fit$a0 + X[c(i, i+1),] %*% fit$beta, ncol = 2)
  names(pred.out) <- rownames(X[c(i, i+1),])
  
  # Saving the betas
  list(pred.out = pred.out, beta.out = fit$beta)
})

# Saving the predicted outcome on the held-out pair
cv.lasso.list <- c(sapply(cv.lasso, function(iter) iter$pred.out))

cor(cv.lasso.list, outcome)

# Saving the betas from each cross-validated run
cv.beta.list <- sapply(cv.lasso, function(iter) t(iter$beta.out))
cv.beta.matrix <- do.call(rbind, cv.beta.list)
colnames(cv.beta.matrix) <- plasma_processed$Metabolite

# Colmeans on the betas
cv.average.beta <- colMeans(cv.beta.matrix)
cv.average.beta[order(abs(cv.average.beta), decreasing = TRUE)]

# -----------------------------------------------------------------------------
# Correlation with DLCO values
# -----------------------------------------------------------------------------

# Loading in the Pittsburgh samples only
pittsburgh_dlco <- read_xlsx("/Users/sarahsamorodnitsky/Dropbox/HIV-COPD/Lavage_Biocrates/Pittsburgh_Samples_DLCO_pp.xlsx")

# Convert the DLCOpp column to numeric - this will induce NAs because there are missing values
pittsburgh_dlco$DLCO_precent_predicted <-
  as.numeric(pittsburgh_dlco$DLCO_precent_predicted)

# Save the Pittsburgh patient IDs with available DLCOpp values
pittsburgh_ids <- pittsburgh_dlco$id[!is.na(pittsburgh_dlco$DLCO_precent_predicted)]

# Subset the Biocrates data to contain just these Pittsburgh patient IDs
plasma_log_scaled_pittsburgh <- plasma_log_scaled[,colnames(plasma_log_scaled) %in% pittsburgh_ids]

# Refine the list of Pittsburgh patient IDs since one Pittsburgh individual was not in the Somascan data
pittsburgh_ids <- pittsburgh_ids[pittsburgh_ids %in% colnames(plasma_log_scaled_pittsburgh)]

# Subset the DLCOpp data to contain just the subjects for whom we had Biocrates data
pittsburgh_dlco_avail <- pittsburgh_dlco %>% dplyr::filter(id %in% pittsburgh_ids)

# Ensure the order of subjects matches with the order of DLCO values
all(pittsburgh_dlco_avail$id == colnames(plasma_log_scaled_pittsburgh)) # TRUE!

# Initialize datafame with columns for protein, test stat, p-value, q-value from 
# Pearson correlation test
dlco_corr <- data.frame(Metabolite = character(nrow(plasma_log_scaled)),
                        TestStat = numeric(nrow(plasma_log_scaled)),
                        PValue = numeric(nrow(plasma_log_scaled)),
                        QValue = numeric(nrow(plasma_log_scaled)))

# Iterate through proteins, test for correlation with DLCO, save in dataframe
for (i in 1:nrow(plasma_log_scaled)) {
  # Calculate Pearson correlation test with DLCO
  res <- cor.test(plasma_log_scaled_pittsburgh[i,],
                  pittsburgh_dlco_avail$DLCO_precent_predicted)
  
  # Save the results
  dlco_corr$Metabolite[i] <- plasma_processed$Metabolite[i]
  dlco_corr$TestStat[i] <- res$statistic
  dlco_corr$PValue[i] <- res$p.value
}

# Calculate multiple comparisons adjustment
dlco_corr$QValue <- p.adjust(dlco_corr$PValue, method = "fdr")

# Order by q-value
dlco_corr <- dlco_corr[order(dlco_corr$QValue, decreasing = FALSE),]

