# -----------------------------------------------------------------------------
# In this script, I process the Plasma Biocrates data and check to make sure
# the processing was done correctly. Following similar steps taken as with
# the lavage data. 
# Sarah Samorodnitsky, with help from Adam Kaplan
# -----------------------------------------------------------------------------

# Setting the working directory
currentwd <- "/Users/sarahsamorodnitsky/Dropbox/HIV-COPD/Plasma_Biocrates/"

# Loading in the data
library(readxl)

# Reading in the lavage data first, this data contains the measurements for
# each of the metabolites. 
plasma_raw_data <- read_excel(paste0(currentwd, 
                                     "Pitt+Van Plasma Biocrates Data.xlsx"))

# Reading in additional subject level data with patient IDs that match the 
# lavage dataset. 
subject_data <- read_excel(paste0(currentwd,
                                  "Pitt-Vancouver FINAL MATCH CORRECTED 7-9-20.xlsx"),
                           n_max = 55)

# Storing metabolite information
metab_info <- plasma_raw_data[-c(1:17), 1:2] # removing blank rows
colnames(metab_info) <- c("Metabolite", "Class") # adding titles

# How many metabolites are there?
nrow(metab_info)

# Extracting the LODs
LODs <- plasma_raw_data[-c(1:16), 17:20]
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

# Which row had all zeros?
which(apply(LODs, 1, function(row) {
  all(row == 0)
}))

# LPC(9:0) 
# 172 (the same metabolite as the Lavage data)

# Storing the sample IDs and associated data
type_ID_species <- as.data.frame(t(plasma_raw_data[c(1, 2, 3, 4, 6, 8), ]))
colnames(type_ID_species) <- c("PlateBarCode", "SampleBarCode", "SampleType",
                               "SampleID", "SubmissionName", "Species")

# I didn't remove the first few columns this data yet because I wanted to keep
# the indices of the subjects correct when I go on to the next step, 
# where I subset the samples out of the quality control samples, etc. 

# Identifying the human samples and selecting them from the dataset
# which samples correspond to humans and lavage samples?
humans_plasma <- which(type_ID_species$SampleType %in% "Sample") 

# Selecting just the human samples corresponding to plasma data. 
plasma.1_plus_info <- plasma_raw_data[, humans_plasma] # selecting the actual samples
plasma_pids <- plasma.1_plus_info[4,] # storing the PIDs
plasma.1 <- plasma.1_plus_info[-c(1:17),] # removing the first 17 rows

# Naming the columns by the PIDs
colnames(plasma.1) <- plasma_pids

# Checking for any duplicate subjects
any(duplicated(colnames(plasma.1))) # false!
ncol(plasma.1) # 47 subjects with plasma data

# Converting the plasma.1 data to numeric
# will throw warnings about converting to NAs
plasma.1 <- apply(plasma.1, 2, as.numeric)

# Have to go down the columns or else it transposes the data. 

# First checking how many people are available in the subject-level data
sum(colnames(plasma.1) %in% subject_data$id) # 47 subjects are available
sum(subject_data$id %in% colnames(plasma.1)) # 47 subjects are available

# Which aren't shared between the two datasets?
missing_ids <- subject_data$id[!(subject_data$id %in% colnames(plasma.1))]
colnames(plasma.1)[!(colnames(plasma.1) %in% subject_data$id)] # there are no subjects missing from the metabolite data
# There are some subject IDs available in the subject-level data but whose plasma
# data isn't available. 

# Matching the order of the patients between the two data sets
# Reordering the metabolite data so the case-control pairs are
# together. This will match the order in the subject-level data.
# c(3,1,2)[match(c(1,2,3), c(3,1,2))] what order to put the 2nd argument
# in to match the first. 
id_match_ordering <- na.omit(match(subject_data$id, colnames(plasma.1))) # need to remove missing obs
plasma.2 <- plasma.1[, id_match_ordering]

# Checking that the reordered columns match the subject level ID names
all(colnames(plasma.2) == subject_data$id[!(subject_data$id %in% missing_ids)]) # TRUE

# Combining the metabolite data with the subject level data
# Adding columns to the plasma_samples dataset that correspond to the metabolite type
plasma.2_metab_info <- cbind(metab_info, plasma.2) # combining the metabolite info with the plasma samples

# Subsetting the subject data to just contain the PIDs from the metabolite data
subject_data.1 <- subject_data[!(subject_data$id %in% missing_ids),]
nrow(subject_data.1) == ncol(plasma.2) # making sure same number of obs before adding metab info

# After doing this, I've removed some subjects within a pair so need to remove
# their partner, as well. 

# Checking which subjects no longer have their pair:
pairs_to_remove <- names(which(table(subject_data.1$setnumber...2) == 1))

# IDs corresponding to the above pairs to remove
ids_missing_pair <- subject_data.1$id[subject_data.1$setnumber...2 %in% pairs_to_remove]

# Removing these IDs from the subject-level data
subject_data.2 <- subject_data.1[!(subject_data.1$id %in% ids_missing_pair), ]

# Removing these IDs from the metabolite data
# The data that includes the metabolite info
plasma_samples.3_metab_info <- plasma.2_metab_info
plasma_samples.3_metab_info <- plasma_samples.3_metab_info[,!(colnames(plasma_samples.3_metab_info) %in% ids_missing_pair)]

# Checking the order matches between the metabolite data and the patient-level data
all(colnames(plasma_samples.3_metab_info[,-c(1:2)]) %in% subject_data.2$id)
all(colnames(plasma_samples.3_metab_info[,-c(1:2)]) == subject_data.2$id)

# -----------------------------------------------------------------------------
# Checking the data has not been changed at all in reformatting. 
# Checking here just the measurement data. 
# -----------------------------------------------------------------------------

# Iterating through each patient and checking that their values still
# match the original data. 

# Just checking plasma_samples_reorder here because the version with the
# metabolite data is the same thing just with a couple extra columns
# appended corresponding to the metabolite data. 

all_match <- c()

for (i in 1:ncol(plasma_samples.3_metab_info[,-c(1:2)])) {
  # storing the current patient i
  patient_id_i <- colnames(plasma_samples.3_metab_info[,-c(1:2)])[i]
  
  # finding that patient id in the raw data
  raw_data_i <- plasma_raw_data[,which(plasma_raw_data[4, ] == patient_id_i)]
  
  # extracting just the lavage data, removing all the extraneous rows
  raw_plasma_data_i <- raw_data_i[-c(1:17), ]
  
  # extracting the processed data corresponding to patient id
  processed_data_i <- plasma_samples.3_metab_info[,colnames(plasma_samples.3_metab_info) == patient_id_i]
  
  # First check all the non-NA values match
  check_matching <- all(as.numeric(raw_plasma_data_i[[1]]) == processed_data_i, na.rm = TRUE)
  
  # Then check the lengths match
  lengths_match <- length(raw_plasma_data_i[[1]]) == length(processed_data_i)
  
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
number_metab_with_missing <- apply(plasma_samples.3_metab_info, 1, function(metabolite) {
  sum(is.na(metabolite))
})

# All the missing values correspond to the same row
which(is.na(plasma_samples.3_metab_info), arr.ind = TRUE)[,1]

# Naming the above vector so I know which metabolite it corresponds to 
names(number_metab_with_missing) <- plasma_samples.3_metab_info$Metabolite

# What proportion the subjects are missing a metabolite?
prop_metab_with_missing <- number_metab_with_missing/(ncol(plasma_samples.3_metab_info)-2) # subtracting off the first two columns which aren't subjects
prop_metab_with_missing[prop_metab_with_missing >= 0.5] # the metabolites for which over 50% of subjects have missing values
metabolites_to_remove <- names(prop_metab_with_missing[prop_metab_with_missing >= 0.5])

# Removing the above metabolites from the plasma samples data:
plasma_samples.4_metab_info <- 
  plasma_samples.3_metab_info[!(plasma_samples.3_metab_info$Metabolite %in% metabolites_to_remove),]

# Also removing the above metabolites from the LODs:

LODs_no_missing <- LODs[!(rownames(LODs) %in% metabolites_to_remove), ]

# Checking the lengths match:
nrow(LODs_no_missing) == nrow(plasma_samples.4_metab_info)

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
number_below_LOD <- sapply(1:nrow(plasma_samples.4_metab_info), function(metab_ind) {
  # removing the first two informational columns
  metab_data <- plasma_samples.4_metab_info[metab_ind, -c(1:2)] 
  
  # storing the corresponding LOD
  metab_LOD <- LODs_no_missing_total[metab_ind]
  sum(metab_data < metab_LOD, na.rm = TRUE) # how many fell below the LOD?
})

# Naming the above vector with the metabolite type
names(number_below_LOD) <- rownames(LODs_no_missing)

# Computing the proportion of subjects who fell below the LOD on each metabolite
prop_subjects_below_LOD <- number_below_LOD/(ncol(plasma_samples.4_metab_info) - 2)

# How many metabolites had over 50% below the LOD?
prop_subjects_below_LOD_g0.5 <- prop_subjects_below_LOD[prop_subjects_below_LOD >= 0.5]
length(prop_subjects_below_LOD_g0.5)
metabolites_to_remove_LOD <- names(prop_subjects_below_LOD_g0.5)

# Total number of metabolites after removing the ones with over 50%
# NAs. 
length(prop_subjects_below_LOD) 
length(metabolites_to_remove_LOD)

# 149 metabolites had over 50% of their observations < LOD out of 408. 

# Removing the metabolites for which over 50% of subjects were < LOD.
plasma_samples.5_metab_info <- 
  plasma_samples.4_metab_info[!(plasma_samples.4_metab_info$Metabolite %in% metabolites_to_remove_LOD), ]

# Checking the row numbers
nrow(plasma_samples.5_metab_info) # new data
nrow(plasma_samples.4_metab_info) # old data

# -----------------------------------------------------------------------------
# Replacing the remaining NAs with 0s. At this point, there are values that
# are still in the data that are below the LOD but I will keep those as is.
# Then there are the metabolites that have some NAs, which I will set to 0,
# as explained in Adam's documentation. 
# -----------------------------------------------------------------------------

# Creating a new dataframe
plasma_samples.6_metab_info <- plasma_samples.5_metab_info

# How many are NAs? There are no NAs which is great. 
sum(is.na(plasma_samples.6_metab_info))

# -----------------------------------------------------------------------------
# Removing any rows with all 0s. 
# -----------------------------------------------------------------------------

# Removing the metabolite information
plasma_samples.6_no_info <- plasma_samples.6_metab_info[, -c(1:2)]

# Checking which rows have only 0s
checking_for_all_0s <- apply(plasma_samples.6_no_info, 1, function(row) all(row == 0))

ind_to_remove_all_0s <- which(checking_for_all_0s)

# Removing this row from the cleaned data:
plasma_samples.7_metab_info <- plasma_samples.6_metab_info
plasma_samples.7_metab_info <- plasma_samples.7_metab_info[-ind_to_remove_all_0s,]

# Checking that no more all 0 rows remain:
any(apply(plasma_samples.7_metab_info[,-c(1:2)], 1, function(row) all(row == 0)))

# Checking the number of rows
nrow(plasma_samples.7_metab_info)

# -----------------------------------------------------------------------------
# Saving the final lavage samples with the metabolite info AND the patient
# level data into one .rda file to be accessed in the future. 
# -----------------------------------------------------------------------------

# # Renaming to a shorter name for the future
# plasma_processed <- plasma_samples.7_metab_info
# subject_processed_plasma <- subject_data.2
# 
# save(plasma_processed, subject_processed_plasma, file = paste0(currentwd, "BiocratesPlasmaProcessed.rda"))
#  