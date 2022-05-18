# -----------------------------------------------------------------------------
# In this script, I process the Lavage Biocrates data and check to make sure
# the processing was done correctly. 
# Sarah Samorodnitsky, with help from Adam Kaplan
# -----------------------------------------------------------------------------

# Setting the working directory
currentwd <- "/Users/sarahsamorodnitsky/Dropbox/HIV-COPD/Lavage_Biocrates/"

# Loading in the data
library(readxl)

# Reading in the lavage data first, this data contains the measurements for
# each of the metabolites. 
lavage_raw_data <- read_excel(paste0(currentwd, 
                                 "Pitt+Van Lavage Biocrates Data.xlsx"))

# Reading in additional subject level data with patient IDs that match the 
# lavage dataset. 
subject_data <- read_excel(paste0(currentwd,
                                   "Pitt-Vancouver FINAL MATCH CORRECTED 7-9-20.xlsx"),
                           n_max = 55)


# Storing metabolite information
metab_info <- lavage_raw_data[-c(1:17), 1:2] # removing blank rows
colnames(metab_info) <- c("Metabolite", "Class") # adding titles

# How many metabolites are there?
nrow(metab_info)

# Extracting the LODs
LODs <- lavage_raw_data[-c(1:16), 17:20]
colnames(LODs) <- LODs[1,] # making the column names the first row
LODs <- LODs[-1,] # removing the first row of column names
LODs <- apply(LODs, 2, as.numeric)

# Adding row names to the LODs so I know which metabolite each corresponds to
rownames(LODs) <- metab_info$Metabolite

# Checking that there is only one non-zero value per row
only_one_non_zero <- apply(LODs, 1, function(row) {
  (sum(row != 0) == 1) | (sum(row != 0) == 0)
})

all(only_one_non_zero)

# There was one row (172) that had all 0s (both the LOD and the metabolite levels). 
# I think this means everything was above the LOD.

# Storing the sample IDs and associated data
type_ID_species <- as.data.frame(t(lavage_raw_data[c(1, 2, 3, 4, 6, 8), ]))
colnames(type_ID_species) <- c("PlateBarCode", "SampleBarCode", "SampleType",
                               "SampleID", "SubmissionName", "Species")

# I didn't remove the first few columns of this data yet because I wanted to keep
# the indices of the subjects correct when I go on to the next step, 
# where I subset the samples out of the quality control samples, etc. 

# Identifying the human samples and selecting them from the dataset
# which samples correspond to humans and lavage samples?
humans_lavage <- which(type_ID_species$SampleType %in% "Sample") 

# Selecting just the human samples corresponding to lavage data (not plasma)
lavage_samples_plus_info <- lavage_raw_data[, humans_lavage] # selecting the actual samples
lavage_pids <- lavage_samples_plus_info[4,] # storing the PIDs
lavage_samples <- lavage_samples_plus_info[-c(1:17),] # removing the first 17 rows

# Naming the columns by the PIDs
colnames(lavage_samples) <- lavage_pids
colnames(lavage_samples_plus_info) <- lavage_pids

# Checking for any duplicate subjects
any(duplicated(colnames(lavage_samples))) # false!
ncol(lavage_samples) # 52 subjects with lavage data

# Converting the lavage_samples data to numeric
# will throw warnings about converting to NAs
lavage_samples <- apply(lavage_samples, 2, as.numeric)

# Have to go down the columns or else it transposes the data. 

# First checking how many people are available in the subject-level data
sum(colnames(lavage_samples) %in% subject_data$id) # 52 subjects are available
sum(subject_data$id %in% colnames(lavage_samples))

# Which aren't shared between the two datasets?
missing_ids <- subject_data$id[!(subject_data$id %in% colnames(lavage_samples))]

# Matching the order of the patients between the two data sets
# Reordering the metabolite data so the case-control pairs are
# together. This will match the order in the subject-level data.
# c(3,1,2)[match(c(1,2,3), c(3,1,2))] what order to put the 2nd argument
# in to match the first. 
id_match_ordering <- na.omit(match(subject_data$id, colnames(lavage_samples))) # need to remove missing obs
lavage_samples_reorder <- lavage_samples[, id_match_ordering]
lavage_samples_plus_info_reorder <- lavage_samples_plus_info[, id_match_ordering]

# Checking that the reordered columns match the subject level ID names
all(colnames(lavage_samples_reorder) == subject_data$id[!(subject_data$id %in% missing_ids)]) # TRUE
all(colnames(lavage_samples_plus_info_reorder) == colnames(lavage_samples_reorder))

# Combining the metabolite data with the subject level data
# Adding columns to the lavage_samples dataset that correspond to the metabolite type
lavage_samples_reorder_metab_info <- cbind(metab_info, lavage_samples_reorder) # combining the metabolite info with the lavage samples

# Subsetting the subject data to just contain the PIDs from the metabolite data
subject_data_avail <- subject_data[!(subject_data$id %in% missing_ids),]
nrow(subject_data_avail) == ncol(lavage_samples_reorder) # making sure same number of obs before adding metab info

# Checked if I removed any subjects within a case-control pair
length(duplicated(subject_data_avail$setnumber...2)[duplicated(subject_data_avail$setnumber...2)]) == nrow(subject_data_avail)/2

# Haven't split up any pairs so don't need to remove any more subjects. 

# -----------------------------------------------------------------------------
# Checking the data has not been changed at all in reformatting. 
# Checking here just the measurement data. 
# -----------------------------------------------------------------------------

# Iterating through each patient and checking that their values still
# match the original data. 

# Just checking lavage_samples_reorder here because the version with the
# metabolite data is the same thing just with a couple extra columns
# appended corresponding to the metabolite data. 

all_match <- c()

for (i in 1:ncol(lavage_samples_reorder)) {
  # storing the current patient i
  patient_id_i <- colnames(lavage_samples_reorder)[i]
  
  # finding that patient id in the raw data
  raw_data_i <- lavage_raw_data[,which(lavage_raw_data[4, ] == patient_id_i)]
  
  # extracting just the lavage data, removing all the extraneous rows
  raw_lavage_data_i <- raw_data_i[-c(1:17), ]
  
  # extracting the processed data corresponding to patient id
  processed_data_i <- lavage_samples_reorder[,colnames(lavage_samples_reorder) == patient_id_i]
  
  # First check all the non-NA values match
  check_matching <- all(as.numeric(raw_lavage_data_i[[1]]) == processed_data_i, na.rm = TRUE)
  
  # Then check the lengths match
  lengths_match <- length(raw_lavage_data_i[[1]]) == length(processed_data_i)
  
  # If all the non-NA values match AND the lengths match, then the remaining values
  # that were checked in line 135 should be NAs and would match. 
  all_match[i] <- check_matching & lengths_match
}

all(all_match)

# -----------------------------------------------------------------------------
# Checking the number of missing and below the LOD (based on the current values 
# of the LOD that I am using)
# -----------------------------------------------------------------------------

# Counting the number of subjects who are missing a metabolite entry
number_metab_with_missing <- apply(lavage_samples_reorder_metab_info, 1, function(metabolite) {
  sum(is.na(metabolite))
})

# Naming the above vector so I know which metabolite it corresponds to 
names(number_metab_with_missing) <- lavage_samples_reorder_metab_info$Metabolite

# What proportion the subjects are missing a metabolite?
prop_metab_with_missing <- number_metab_with_missing/(ncol(lavage_samples_reorder_metab_info)-2) # subtracting off the first two columns which aren't subjects
prop_metab_with_missing[prop_metab_with_missing >= 0.5] # the metabolites for which over 50% of subjects have missing values
metabolites_to_remove <- names(prop_metab_with_missing[prop_metab_with_missing >= 0.5])

# Removing the above metabolites from the lavage samples data:

lavage_samples_reorder_metab_info_no_missing <- 
  lavage_samples_reorder_metab_info[!(lavage_samples_reorder_metab_info$Metabolite %in% metabolites_to_remove),]

# Also removing the above metabolites from the LODs:

LODs_no_missing <- LODs[!(rownames(LODs) %in% metabolites_to_remove), ]

# Checking the lengths match:
nrow(LODs_no_missing) == nrow(lavage_samples_reorder_metab_info_no_missing)

# -----------------------------------------------------------------------------
# Checking how many values were below the LOD. 
# To compute the LOD, I sum over the rows because there is at most 1 non-zero
# entry in each row, corresponding to the operating method used to run the 
# metabolite. 
# -----------------------------------------------------------------------------

# Summing over the rows to get the LOD to use to check which values 
# fell below the LOD. 
LODs_no_missing_total <- apply(LODs_no_missing, 1, sum)

# Just to compare:
LODs_no_missing_with_totals <- cbind.data.frame(LODs_no_missing, LODs_no_missing_total)

# Now, counting the number of subjects whose measurement fell below the LOD:
number_below_LOD <- sapply(1:nrow(lavage_samples_reorder_metab_info_no_missing), function(metab_ind) {
  # removing the first two informational columns
  metab_data <- lavage_samples_reorder_metab_info_no_missing[metab_ind, -c(1:2)] 
  
  # storing the corresponding LOD
  metab_LOD <- LODs_no_missing_total[metab_ind]
  sum(metab_data < metab_LOD, na.rm = TRUE) # how many fell below the LOD?
})

# Naming the above vector with the metabolite type
names(number_below_LOD) <- rownames(LODs_no_missing)

# Computing the proportion of subjects who fell below the LOD on each metabolite
prop_subjects_below_LOD <- number_below_LOD/(ncol(lavage_samples_reorder_metab_info_no_missing) - 2)

# How many metabolites had over 50% below the LOD?
prop_subjects_below_LOD_g0.5 <- prop_subjects_below_LOD[prop_subjects_below_LOD >= 0.5]
length(prop_subjects_below_LOD_g0.5)
metabolites_to_remove_LOD <- names(prop_subjects_below_LOD_g0.5)

# Total number of metabolites after removing the ones with over 50%
# NAs. 
length(metabolites_to_remove_LOD) 

# 151 metabolites had over 50% of their observations < LOD,
# out of 404. 

# Removing the metabolites for which over 50% of subjects were < LOD.
lavage_samples_reorder_metab_info_no_missing_above_LOD <- 
  lavage_samples_reorder_metab_info_no_missing[!(lavage_samples_reorder_metab_info_no_missing$Metabolite %in% metabolites_to_remove_LOD), ]


# -----------------------------------------------------------------------------
# Replacing the remaining NAs with 0s. At this point, there are values that
# are still in the data that are below the LOD but I will keep those as is.
# Then there are the metabolites that have some NAs, which I will set to 0,
# as explained in Adam's documentation. 
# -----------------------------------------------------------------------------

# Creating a new dataframe
lavage_samples_reorder_metab_info_no_missing_above_LOD_no_zeros <- lavage_samples_reorder_metab_info_no_missing_above_LOD

# How many are NAs? Only 34, not that many. 
sum(is.na(lavage_samples_reorder_metab_info_no_missing_above_LOD_no_zeros))

# Replacing the NAs with 0s.
lavage_samples_reorder_metab_info_no_missing_above_LOD_no_zeros[is.na(lavage_samples_reorder_metab_info_no_missing_above_LOD_no_zeros)] <- 0

# -----------------------------------------------------------------------------
# Removing any rows with all 0s. 
# -----------------------------------------------------------------------------

# Removing the metabolite information
lavage_samples_reorder_metab_info_no_missing_above_LOD_no_zeros_no_info <- lavage_samples_reorder_metab_info_no_missing_above_LOD_no_zeros[, -c(1:2)]

# Checking which rows have only 0s
checking_for_all_0s <- apply(lavage_samples_reorder_metab_info_no_missing_above_LOD_no_zeros_no_info, 1, 
      function(row) all(row == 0))

ind_to_remove_all_0s <- which(checking_for_all_0s)

# Removing this row from the cleaned data:
lavage_samples_reorder_metab_info_no_missing_above_LOD_no_zeros_no_all_0_rows <- lavage_samples_reorder_metab_info_no_missing_above_LOD_no_zeros
lavage_samples_reorder_metab_info_no_missing_above_LOD_no_zeros_no_all_0_rows <- lavage_samples_reorder_metab_info_no_missing_above_LOD_no_zeros_no_all_0_rows[-ind_to_remove_all_0s,]

# Checking that no more all 0 rows remain:
any(apply(lavage_samples_reorder_metab_info_no_missing_above_LOD_no_zeros_no_all_0_rows, 1, 
      function(row) all(row == 0)))

# -----------------------------------------------------------------------------
# Saving the final lavage samples with the metabolite info AND the patient
# level data into one .rda file to be accessed in the future. 
# -----------------------------------------------------------------------------

# Renaming to a shorter name for the future
lavage_processed <- lavage_samples_reorder_metab_info_no_missing_above_LOD_no_zeros_no_all_0_rows
subject_processed <- subject_data_avail

save(lavage_processed, subject_processed, file = paste0(currentwd, "BiocratesLavageProcessed.rda"))
