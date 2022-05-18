# -----------------------------------------------------------------------------
# In this script, I do some exploratory data analysis on the lavage samples
# from the Biocrates data. 
# Sarah Samorodnitsky
# -----------------------------------------------------------------------------

library(tidyverse)

# Setting the working directory
currentwd <- "/Users/sarahsamorodnitsky/Dropbox/HIV-COPD/Lavage_Biocrates/"

# Loading in the cleaned metabolite and patient data (loads in both)
load(paste0(currentwd, "BiocratesLavageProcessed.rda"))

# -----------------------------------------------------------------------------
# Creating histograms of the scaled and centered data and also the log-
# transformed versions. 
# -----------------------------------------------------------------------------

# Creating a new dataset to store the scaled data in 
lavage_scaled <- lavage_processed[,3:ncol(lavage_processed)]
lavage_log_scaled <- lavage_processed[,3:ncol(lavage_processed)]

# Scaling all the rows because we want each metabolite to be centered, not each subject. 
lavage_scaled <- t(apply(lavage_scaled, 1, function(x) scale(x, center = TRUE, scale = TRUE)))
lavage_log_scaled <- t(apply(lavage_log_scaled, 1, function(x) scale(log(x+1), center = TRUE, scale = TRUE)))

# Adding the column names back
colnames(lavage_scaled) <- colnames(lavage_processed)[-c(1:2)]
colnames(lavage_log_scaled) <- colnames(lavage_processed)[-c(1:2)]

# -----------------------------------------------------------------------------
# Creating a PCA plot of the metabolite data.
# -----------------------------------------------------------------------------

# First, running an SVD on the metabolite data. 
# Using the logged and scaled data

# Transpose the matrix so that the rows are the patients
# Then the UD part of the SVD corresponds to patient-level variation. 
lavage_svd <- svd(t(lavage_log_scaled))

# -----------------------------------------------------------------------------
# Demographics table
# -----------------------------------------------------------------------------

library(table1)

# Creating a copy of the clinical data to add in factor levels
subject_processed_v2 <- subject_processed

# Relabeling the levels -- 

# Ethnicity/race
subject_processed_v2$ethnicity <- as.factor(subject_processed_v2$ethnicity)
levels(subject_processed_v2$ethnicity) <- c("Black, Non-Hispanic",
                                             "White, Hispanic/Latino",
                                             "Asian/Pacific Islander")

# Sex
subject_processed_v2$sex <- as.factor(subject_processed_v2$sex)
levels(subject_processed_v2$sex) <- c("Male", "Female")

# Smoker
subject_processed_v2$smoker <- as.factor(subject_processed_v2$smoker)
levels(subject_processed_v2$smoker) <- c("Yes", "No")

# ART
subject_processed_v2$art <- as.factor(subject_processed_v2$art)
levels(subject_processed_v2$art) <- c("Yes", "No")

# Case:Control (cases = 1, controls = 2)
subject_processed_v2$ccstat <- as.factor(subject_processed_v2$ccstat)
levels(subject_processed_v2$ccstat) <- c("Case", "Control")

# Renaming the columns
colnames(subject_processed_v2)[colnames(subject_processed_v2) == "art"] <- "ART"
colnames(subject_processed_v2)[colnames(subject_processed_v2) == "age"] <- "Age"
colnames(subject_processed_v2)[colnames(subject_processed_v2) == "smoker"] <- "Smoker"
colnames(subject_processed_v2)[colnames(subject_processed_v2) == "ethnicity"] <- "Ethnicity"
colnames(subject_processed_v2)[colnames(subject_processed_v2) == "sex"] <- "Sex"
colnames(subject_processed_v2)[colnames(subject_processed_v2) == "FEV1_percent_predicted"] <- "FEV1 Percent Predicted"

# -----------------------------------------------------------------------------
# Univariate testing: using paired t-tests to compare metabolite levels
# between cases and controls in each matched pair.
# Using a false-discovery rate correction after repeated testing. 
# -----------------------------------------------------------------------------

# Creating a dataframe to store the results
case_control_t.test_res <- data.frame(test.stat = numeric(nrow(lavage_log_scaled)),
                                      p.value = numeric(nrow(lavage_log_scaled)))

for (i in 1:nrow(case_control_t.test_res)) { # iterating through each metabolite
  # the current metabolite
  metabolite_i <- lavage_processed$Metabolite[i]
  
  # extracting the case and control ids
  case_ids <- subject_processed$id[subject_processed$ccstat == 1]
  control_ids <- subject_processed$id[subject_processed$ccstat == 2]
  
  # subsetting the metabolite i levels for cases and controls
  case_metab_i <- lavage_log_scaled[i, which(colnames(lavage_log_scaled) %in% case_ids)]
  control_metab_i <- lavage_log_scaled[i, which(colnames(lavage_log_scaled) %in% control_ids)]
  
  # performing the paired t-test
  res_i <- t.test(case_metab_i, control_metab_i, paired = TRUE)
  
  # storing the results
  case_control_t.test_res$test.stat[i] <- res_i$statistic
  case_control_t.test_res$p.value[i] <- res_i$p.value
}

# Performing a false discovery rate correction
# Note that the order in which the q-values are returned is the same
# order in which the p-values are given. 
p.values <- case_control_t.test_res$p.value
names(p.values) <- lavage_processed$Metabolite
case_control_t.test_fdr_adjust <- p.adjust(p.values, method = "fdr")

# -----------------------------------------------------------------------------
# Performing DWD on the metabolite data to compare cases and controls. 
# This was the standard DWD that I initially started with. 
# -----------------------------------------------------------------------------

library("SparseM")
library("DWDLargeR")

# # Creating a vector of case-control labels, either -1 or 1
case_control_labels <- ifelse(subject_processed$ccstat == 1, 1, -1) # change cases to 1, controls to -1
# 
# Changing the form of the data (not sure why Adam did this)
lavage_log_scaled.DWD <- SparseM::`[.matrix.csr`(lavage_log_scaled)
# 
# # Adjusting the penalty parameter
penalty.lavage <- DWDLargeR::penaltyParameter(X = lavage_log_scaled.DWD,
                                              y = case_control_labels,
                                              expon = 1)
# 
# # Cross-validated DWD:
# # Iterating through each case-control pair
# ind_of_pairs <- seq(from = 1, to = length(case_control_labels), by = 2)
# cv.scores <- lapply(ind_of_pairs, function(i) {
#   # For each case-control pair, i
#   w.vec = genDWD(X = lavage_log_scaled.DWD[,-c(i, i+1)],
#                  y = case_control_labels[-c(i, i+1)],
#                  expon = 1,
#                  C = penalty.lavage,
#                  scaleFea = FALSE)$w
#   cv.scores = w.vec%*%lavage_log_scaled.DWD[,c(i, i+1)]
#   cv.scores
# })
# 
# # Doing a t-test on the scores 
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

full_training_res <- genDWD(X = lavage_log_scaled.DWD,
                           y = case_control_labels,
                           expon = 1,
                           C = penalty.lavage,
                           scaleFea = FALSE)

# -----------------------------------------------------------------------------
# Implementing a permutation testing-based approach to DWD.
# -----------------------------------------------------------------------------

library("diproperm")
library("SparseM")
library("DWDLargeR")
library("gtools")
# 
# # Creating a vector of case-control labels, either -1 or 1
# case_control_labels <- ifelse(subject_processed$ccstat == 1, 1, -1) # change cases to 1, controls to -1
# 
# # Changing the form of the data (not sure why Adam did this)
# lavage_log_scaled.DWD <- SparseM::`[.matrix.csr`(lavage_log_scaled)
# 
# # Adjusting the penalty parameter
# penalty.lavage <- DWDLargeR::penaltyParameter(X = lavage_log_scaled.DWD,
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
# lavage_log_scaled.DWD <- SparseM::`[.matrix.csr`(lavage_log_scaled)
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
# DWD_t.stats_permute <- c()
# 
# for (permute_iter in 1:npermute) {
#   # Randomly permute which subject in each pair is a case and which is a control
#   case_control_labels.permute.list <- lapply(case_control_labels.list, function(pair) permute(pair))
#   case_control_labels.permute <- unlist(case_control_labels.permute.list)
# 
#   # LOOCV DWD on pairs
#   cv.scores.permute <- lapply(ind_of_pairs, function(i) {
#     # For each case-control pair, i
#     w.vec = genDWD(X = lavage_log_scaled.DWD[,-c(i, i+1)],
#                    y = case_control_labels.permute[-c(i, i+1)],
#                    expon = 1,
#                    C = penalty.lavage,
#                    scaleFea = FALSE)$w
#     cv.scores = w.vec%*%lavage_log_scaled.DWD[,c(i, i+1)]
#     cv.scores
#   })
# 
#   # Storing the scores in a matrix to do a paired t-test
#   cv.scores.permute.mat <- matrix(unlist(cv.scores.permute), ncol = 2, byrow = T)
#   colnames(cv.scores.permute.mat) <- c("case", "control")
# 
#   # Doing a t-test on the scores
#   # Unlist the cv.scores so that each case and control entry within a pair is adjacent.
#   case_control_dwd_t.test.permute <- t.test(x = cv.scores.permute.mat[,1],
#                                             y = cv.scores.permute.mat[,2],
#                                             paired = TRUE,
#                                             alternative = "greater")
# 
#   # Storing the results
#   DWD_t.stats_permute[permute_iter] <- case_control_dwd_t.test.permute$statistic
# }
# 
# save(DWD_t.stats_permute, file = paste0(currentwd, "Permutation_Testing_DWD_Test_Stats.rda"))

# -----------------------------------------------------------------------------
# Summing across the metabolites for each metabolite family and comparing cases 
# and controls. This gives us a sense of whether the cumulative levels of each
# metabolite differ between cases and controls. 
# The idea here is to sum within each subject, so retain the subjects and the
# metabolite families, and then compare between cases and controls. 
# -----------------------------------------------------------------------------

# First, adding in the metabolite information back
lavage_log_scaled_info <- cbind(lavage_processed[,1:2], lavage_log_scaled)

# Splitting up the data into cases and controls
lavage_log_scaled_cases <- lavage_log_scaled_info[,colnames(lavage_log_scaled_info) %in% subject_processed[subject_processed$ccstat == 1,]$id,]
lavage_log_scaled_cases <- cbind(lavage_processed[,1:2], lavage_log_scaled_cases)
lavage_log_scaled_controls <- lavage_log_scaled_info[,colnames(lavage_log_scaled_info) %in% subject_processed[subject_processed$ccstat == 2,]$id,]
lavage_log_scaled_controls <- cbind(lavage_processed[,1:2], lavage_log_scaled_controls)

# Which are the families of metabolits?
metabolite_families <- unique(lavage_processed$Class)
n_families <- length(metabolite_families)

# Adding within the families but retaining the subjects
lavage_log_scaled_cases_cumulative <- matrix(nrow = n_families,
                                             ncol = ncol(lavage_log_scaled_cases) - 2)

lavage_log_scaled_controls_cumulative <- matrix(nrow = n_families,
                                                ncol = ncol(lavage_log_scaled_cases) - 2)

for (i in 1:n_families) {
  metabolite_data_cases <- lavage_log_scaled_cases[lavage_log_scaled_cases$Class == metabolite_families[i], -c(1:2)]
  metabolite_data_controls <- lavage_log_scaled_controls[lavage_log_scaled_controls$Class == metabolite_families[i], -c(1:2)]
  
  lavage_log_scaled_cases_cumulative[i,] <- apply(metabolite_data_cases, 2, sum)
  lavage_log_scaled_controls_cumulative[i,] <- apply(metabolite_data_controls, 2, sum)
}

# Doing a paired t-test for each family
metabolite_casecontrol_cumulative <- data.frame(Family = metabolite_families,
                                                Test.Stat = numeric(n_families),
                                                P.value = numeric(n_families))

for (i in 1:n_families) {
  res_metab_i <- t.test(lavage_log_scaled_cases_cumulative[i,], 
                        lavage_log_scaled_controls_cumulative[i,], 
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
all(colnames(lavage_log_scaled) == subject_processed$id)

# Computing the correlations between each metabolite (each row in the lavage data)
# with FEV1_percent_predicted
n_metabolite <- nrow(lavage_log_scaled)
corr_fev1pp_metabolites <- data.frame(Metabolite = character(n_metabolite),
                                      Correlation = numeric(n_metabolite),
                                      P.Value = numeric(n_metabolite),
                                      Q.Value = numeric(n_metabolite))

for (i in 1:n_metabolite) {
  # Save the current metabolite
  metab_i <- lavage_processed$Metabolite[i]
    
  # Calculate correlation
  corr_i <- cor.test(lavage_log_scaled[i,],
                      subject_processed$FEV1_percent_predicted)
  
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
sum(corr_fev1pp_metabolites$Q.Value < 0.01)
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
outcome <- subject_processed$FEV1_percent_predicted

# Full training data
X <- t(lavage_log_scaled)
Y <- outcome

# Fit the model
lasso_cv <- cv.glmnet(X, Y, family = "gaussian", alpha = 1)
lambda.min <- lasso_cv$lambda.min

# Fit the model again on the test data with the min lambda
lasso_final <- glmnet(X, Y, family = "gaussian", alpha = 1, lambda = lambda.min)

# Which coefficients are non-zero?
lasso_beta <- matrix(lasso_final$beta)
rownames(lasso_beta) <- lavage_processed$Metabolite
non_zero_betas <- lasso_beta[lasso_beta != 0, ]


# Consider cross validation --
case_control_labels <- ifelse(subject_processed$ccstat == 1, 1, -1)
ind_of_pairs <- seq(from = 1, to = length(case_control_labels), by = 2)

cv.lasso <- lapply(ind_of_pairs, function(i) {
  # For each case-control pair, i
  # find optimal penalty
  cv.fit <- cv.glmnet(x = X[-c(i, i+1),],
                   y = Y[-c(i, i+1)],
                   family = "gaussian",
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

# Saving the betas from each cross-validated run
cv.beta.list <- sapply(cv.lasso, function(iter) t(iter$beta.out))
cv.beta.matrix <- do.call(rbind, cv.beta.list)
colnames(cv.beta.matrix) <- lavage_processed$Metabolite

# Calculate the average coefficient for each metabolite
cv.average.beta <- colMeans(cv.beta.matrix)

# Calculate the proportion of cross validation iterations each metabolite was selected
cv.prop.nonzero <- apply(cv.beta.matrix, 2, function(metabolite) sum(metabolite != 0)/length(metabolite))

# -----------------------------------------------------------------------------
# Pathway analysis for metabolites combined with proteins
# -----------------------------------------------------------------------------

# Loading in the original Biocrates data to obtain the ChEBI IDs
lavage_raw_data <- read_excel(paste0(currentwd, "Pitt+Van Lavage Biocrates Data.xlsx"))
metabolite_ids_raw_data <- lavage_raw_data[17:nrow(lavage_raw_data), 1:3]
metabolite_ids_raw_data[1,1] <- "Metabolite Name"
colnames(metabolite_ids_raw_data) <- metabolite_ids_raw_data[1,]
metabolite_ids_raw_data <- metabolite_ids_raw_data[-1,]

# The BIO ID column contains the CHEBI ID that can go into pathway analysis

# Saving the IDs for metabolites that were significant at the 0.05 level
significant_metabolites <- p.values[p.values < 0.05]

# Obtaining the significant metabolites' CHEBI ID
significant_metabolites_CHEBI <- metabolite_ids_raw_data %>% 
  filter(`Metabolite Name` %in% names(significant_metabolites)) %>% 
  dplyr::select(`Bio ID`) %>% unlist 

# Check -- 
length(significant_metabolites) == length(significant_metabolites_CHEBI)

Metabolite_ID_file <- file(paste0(currentwd, "Metabolite_IDs_for_Significant_Metabolites.txt"))
writeLines(significant_metabolites_CHEBI, Metabolite_ID_file)
close(Metabolite_ID_file)

# Obtaining the CHEBI ID for all metabolites considered in the analysis
reference_metabolites_CHEBI <- metabolite_ids_raw_data %>% 
  filter(`Metabolite Name` %in% lavage_processed$Metabolite) %>% 
  dplyr::select(`Bio ID`) %>% unlist 

# Check --
length(reference_metabolites_CHEBI) == length(lavage_processed$Metabolite)

# Saving the reference file to upload
Metabolite_ID_reference_file <- file(paste0(currentwd, "Metabolite_IDs_for_Reference.txt"))
writeLines(reference_metabolites_CHEBI, Metabolite_ID_reference_file)
close(Metabolite_ID_reference_file)

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
lavage_log_scaled_pittsburgh <- lavage_log_scaled[,colnames(lavage_log_scaled) %in% pittsburgh_ids]

# Refine the list of Pittsburgh patient IDs since one Pittsburgh individual was not in the Somascan data
pittsburgh_ids <- pittsburgh_ids[pittsburgh_ids %in% colnames(lavage_log_scaled_pittsburgh)]

# Subset the DLCOpp data to contain just the subjects for whom we had Biocrates data
pittsburgh_dlco_avail <- pittsburgh_dlco %>% filter(id %in% pittsburgh_ids)

# Ensure the order of subjects matches with the order of DLCO values
all(pittsburgh_dlco_avail$id == colnames(lavage_log_scaled_pittsburgh)) # TRUE!

# Initialize datafame with columns for protein, test stat, p-value, q-value from 
# Pearson correlation test
dlco_corr <- data.frame(Metabolite = character(nrow(lavage_log_scaled)),
                        TestStat = numeric(nrow(lavage_log_scaled)),
                        PValue = numeric(nrow(lavage_log_scaled)),
                        QValue = numeric(nrow(lavage_log_scaled)))

# Iterate through proteins, test for correlation with DLCO, save in dataframe
for (i in 1:nrow(lavage_log_scaled)) {
  # Calculate Pearson correlation test with DLCO
  res <- cor.test(lavage_log_scaled_pittsburgh[i,],
                  pittsburgh_dlco_avail$DLCO_precent_predicted)
  
  # Save the results
  dlco_corr$Metabolite[i] <- lavage_processed$Metabolite[i]
  dlco_corr$TestStat[i] <- res$statistic
  dlco_corr$PValue[i] <- res$p.value
}

# Calculate multiple comparisons adjustment
dlco_corr$QValue <- p.adjust(dlco_corr$PValue, method = "fdr")

# Order by q-value
dlco_corr <- dlco_corr[order(dlco_corr$QValue, decreasing = FALSE),]

# Investigate specific metabolites we found were significant
metabolites <- c("PC(31:3)", "PC(31:0)", "PC(34:4)", "PC(34:3)", "TG(55:9)",
                 "AC(14:1)", "AC(18:0)", "PC(32:0)", "AC(10:0)", "PC(33:0)",
                 "PC(41:2)", "CE(17:0)", "PC-O(33:6)")

# Select the rows corresponding to above metabolites
dlco_corr_for_metabolites <- dlco_corr %>% filter(Metabolite %in% metabolites)

# Redoing the FDR calculation on just the above subset
dlco_corr_for_metabolites$QValue <- p.adjust(dlco_corr_for_metabolites$PValue, method = "fdr")
