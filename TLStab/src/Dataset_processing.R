library(data.table)
library(magrittr)
library(caret)
library(stringr)
library(Peptides)
library(Biostrings)
library(MatchIt)
library(CBPS)
library(designmatch)

data <- fread("./mhc_ligand_full.csv", skip=1, header=T)
names(data)

# Focusing on the Assay Group which si the types of data existing in IEDB: 
data[, `Assay Group`] %>% unique

# We need the values to be half life to predict stability (and also focusing on the class I alleles):
stability_data <- data[`Assay Group` == "half life" & `MHC allele class` == "I"]
rm(data)
gc()
gc()

# Load pseudosequences
# The alleles in the dataset NEED to have pseudosequences!
mhc_pseudoseq <- fread("./updated_pseudosequences.csv")
allele_list <- mhc_pseudoseq[, allele] %>% unique
stability_data <- stability_data[`Allele Name` %in% allele_list]

# All units are in minutes, no further conversion needs to be made
stability_data[, Units] %>% unique

# Looking at the Measurement inequality column, there seem to be a lot of datapoints that do not have the field specified
stability_data[, .N, by = `Measurement Inequality`]

filter_by_equality <- function(data, equalities_vector, study=FALSE, IRI=0) {
  
  # In terms of other specified fields, most are equalities, so it's probably safe to remove all the non-equalities:
  data <- data[`Measurement Inequality` %in% equalities_vector]
  data <- na.omit(data, cols=c("Quantitative measurement"))
  
  # Use the transformation defined in the NetMHCstabpan paper:
  # Note that the 1 in the formula should be different for every allele. For simplification purposes, we are going to leave it as 1.
  # It's something to discuss about in relation to binding affinity, some HLA alleles are more "stable" in general than others.
  if ((IRI != "http://www.iedb.org/reference/1014192") & (IRI != "http://www.iedb.org/reference/1014194")) {
    data <- data[, `Quantitative measurement (h)` := 2^-(1/(`Quantitative measurement`/60))]
  }  
  else {
    data <- data[, `Quantitative measurement (h)` := 2^-(1/`Quantitative measurement`)]
  }
  # Throw away everything else exept the allele name, the peptide and the stability measurement
  train_data <- data[Authors == "Michael Rasmussen M.Sc; Mikkel Nors Harndahl Ph.D.; Anne Bregnballe Kristensen B.Sc.; Ida Kallehauge Nielsen B.Sc.; Kasper W Jorgensen Ph.D.; Anette Stryhn Ph.D.; Morten Nielsen Ph.D.; S&ouml;ren Buus Buus MD, Ph.D."]
  test_data <- data[Authors != "Michael Rasmussen M.Sc; Mikkel Nors Harndahl Ph.D.; Anne Bregnballe Kristensen B.Sc.; Ida Kallehauge Nielsen B.Sc.; Kasper W Jorgensen Ph.D.; Anette Stryhn Ph.D.; Morten Nielsen Ph.D.; S&ouml;ren Buus Buus MD, Ph.D."]
  if (study == TRUE) {
    test_data <- test_data[`Reference IRI` == IRI]
  }
  print(dim(test_data)[1])
  train_data <- train_data[str_length(Description) <= 15 & str_length(Description) > 7, .(`Allele Name`, Description, `Quantitative measurement (h)`)]
  test_data <- test_data[str_length(Description) <= 15 & str_length(Description) > 7, .(`Allele Name`, Description, `Quantitative measurement (h)`)]
  print(dim(test_data)[1])
  # Test data has some duplicates, and duplicates (at least for now) we will remove:
  duplicated_test_data <- test_data[duplicated(test_data, by = c("Allele Name", "Description"))]
  duplicated_test_data <- test_data[duplicated_test_data, on = c("Allele Name", "Description")][, .(`Allele Name`, `Description`, `Quantitative measurement (h)`)]

  # 2a) When we have only equalities, average the results
  duplicated_test_data_to_keep <- duplicated_test_data[, lapply(.SD, mean), 
                                                       by = .(`Allele Name`, `Description`), 
                                                       .SDcols = c("Quantitative measurement (h)")]
  test_data <- fsetdiff(test_data, duplicated_test_data)
  test_data <- rbindlist(list(test_data, duplicated_test_data_to_keep))
  print(dim(test_data)[1])
  existing_in_training <- test_data[train_data, on = .(`Allele Name`, `Description`), nomatch = 0][, .(`Allele Name`, `Description`, `Quantitative measurement (h)`)]
  test_data <- fsetdiff(test_data, existing_in_training)
  print(dim(test_data)[1])
  # Rename the columns as you know: (allele, peptide, measurement_value)
  names(train_data) <- c("allele", "peptide", "measurement_value")
  names(test_data) <- c("allele", "peptide", "measurement_value")
  
  list(train_data, test_data)
}

train_test_list <- filter_by_equality(stability_data, c("="))
train_data <- train_test_list[[1]]
test_data <- train_test_list[[2]]

fwrite(test_data, "./cleaned_data/cleaned_test_data_small.csv")

## Per-study

# Rebecca (199 peptides)

train_test_list <- filter_by_equality(stability_data, c("=",""), study=TRUE, IRI="http://www.iedb.org/reference/1026840")
train_data <- train_test_list[[1]]
test_data <- train_test_list[[2]]

fwrite(test_data, "./cleaned_data/cleaned_test_data_Rebecca.csv")

# Pavlo (91 peptides, maybe leave those out?)

train_test_list <- filter_by_equality(stability_data, c("=",""), study=TRUE, IRI="http://www.iedb.org/reference/1026371")
train_data <- train_test_list[[1]]
test_data <- train_test_list[[2]]

fwrite(test_data, "./cleaned_data/cleaned_test_data_Pavlo.csv")

# Anette (67 peptides, maybe leave those out too?)

train_test_list <- filter_by_equality(stability_data, c("=",""), study=TRUE, IRI="http://www.iedb.org/reference/1037619")
train_data <- train_test_list[[1]]
test_data <- train_test_list[[2]]

fwrite(test_data, "./cleaned_data/cleaned_test_data_Anette.csv")

# Unpublished Ebola

train_test_list <- filter_by_equality(stability_data, c("=",""), study=TRUE, IRI="http://www.iedb.org/reference/1014192")
train_data <- train_test_list[[1]]
test_data <- train_test_list[[2]]

fwrite(test_data, "./cleaned_data/cleaned_test_data_Ebola.csv")

# Unpublished PoxVirus

train_test_list <- filter_by_equality(stability_data, c("=",""), study=TRUE, IRI="http://www.iedb.org/reference/1014194")
train_data <- train_test_list[[1]]
test_data <- train_test_list[[2]]

fwrite(test_data, "./cleaned_data/cleaned_test_data_PoxVirus.csv")

# Whole thing

train_test_list <- filter_by_equality(stability_data, c("=",""))
train_data <- train_test_list[[1]]
test_data <- train_test_list[[2]]

fwrite(train_data, "./cleaned_data/cleaned_train_data.csv")
fwrite(test_data, "./cleaned_data/cleaned_test_data.csv")

all_data <- rbindlist(list(train_data, test_data))

# Are there any duplicates in all of the data?
dim(unique(all_data, by = c("allele", "peptide")))[1]
dim(all_data)[1]

folds <- 10
cvIndex <- createFolds(train_data$`measurement_value`, k = folds, list = FALSE, returnTrain = FALSE)

for (i in seq(folds)) {
  fwrite(train_data[cvIndex == i], paste0("./folds/test_", i, ".csv"))
  fwrite(train_data[cvIndex != i], paste0("./folds/train_", i, ".csv"))
}

for (i in seq(folds)) {
  temp_train_fold <- fread(paste0("./folds/train_", i, ".csv"))
  cvIndex <- createFolds(temp_train_fold$`measurement_value`, k = folds, list = FALSE, returnTrain = FALSE)
  for (i in seq(folds)) {
    fwrite(temp_train_fold[cvIndex != i], paste0("./folds/", i, "/train_", i, ".csv"))
    fwrite(temp_train_fold[cvIndex == i], paste0("./folds/", i, "/val_", i, ".csv"))
  }
}

# Function that encodes the peptide into the MHCFlurry2.0 format:
sequences_to_fixed_length_index_encoded_array <- function(sequence) {
  X_length <- 15 - str_length(sequence)
  if (X_length > 0) {
    return (paste0(sequence, paste0(rep('X', X_length + floor(X_length / 2)), collapse = ""), sequence, 
                   paste0(rep('X', ceiling(X_length / 2) + X_length), collapse = ""), sequence))
  } else {
    return (paste0(sequence, sequence, sequence))
  }
}

# Apply the function to the datasets:
all_data[, peptide := lapply(peptide, sequences_to_fixed_length_index_encoded_array)]

# Featurization function (based on BLOSUM62)
featurization <- function(dataset, pseudosequences, filename) {
  
  # Make the feature string comprising of both the extended peptide sequence and the MHC pseudosequence
  setkey(dataset, allele)
  setkey(pseudosequences, allele)
  dataset <- dataset[mhc_pseudoseq, nomatch=0]
  dataset[, features := paste0(peptide, pseudosequence)]
  
  feature_list <- dataset[, features] %>% unlist
  
  peptide_list <- feature_list %>% strsplit("")
  allele_list <- dataset[, allele]
  
  i <- 0
  features_conc <- sapply(peptide_list, function(aa_seq) {
    i <<- i + 1
    features <- lapply(aa_seq, function(amino_acid) {
      xx <- BLOSUM62[amino_acid, ]
      names(xx) <- NULL
      xx
    }) %>% unlist()
  }) %>% t() %>% as.data.table()
  setnames(features_conc, paste0(rep("V", ncol(features_conc)), seq(ncol(features_conc))))
  features_conc[, peptide := substr(sub("\\X.*", "", feature_list), 1, 15)]
  features_conc[, pseudosequence := substr(feature_list, nchar(feature_list) - 34 + 1, nchar(feature_list))]
  features_conc[, allele := allele_list]
  setcolorder(features_conc, c("peptide", "allele", "pseudosequence", colnames(features_conc)[1:(length(colnames(features_conc)) - 3)]))
  fwrite(features_conc, file = filename, sep = ",")
}

# Load BLOSUM62 matrix
data(BLOSUM62)
BLOSUM62 <- BLOSUM62[-c(21, 22, 23, 25), -c(21, 22, 23, 25)]
BLOSUM62 <- apply(BLOSUM62, MARGIN = 2, FUN = function(X) (X - min(X))/diff(range(X)))
BLOSUM62['X', ] <- rep.int(0, 21)
BLOSUM62[, 'X'] <- rep.int(0, 21)
BLOSUM62['X', 'X'] <- 1

featurization(all_data, mhc_pseudoseq, "./features/Stab_BLOSUM62_features.csv")

# Featurization function (based on BLOSUM62 PLUS Binding Affinity predictions from netMHCpan4.1)
featurization_BA <- function(dataset, pseudosequences, BA_file, filename) {
  
  # Make the feature string comprising of both the extended peptide sequence and the MHC pseudosequence
  setkey(dataset, allele)
  setkey(pseudosequences, allele)
  dataset <- dataset[mhc_pseudoseq, nomatch=0]
  dataset[, features := paste0(peptide, pseudosequence)]
  
  feature_list <- dataset[, features] %>% unlist
  
  peptide_list <- feature_list %>% strsplit("")
  allele_list <- dataset[, allele]
  
  i <- 0
  features_conc <- sapply(peptide_list, function(aa_seq) {
    i <<- i + 1
    features <- lapply(aa_seq, function(amino_acid) {
      xx <- BLOSUM62[amino_acid, ]
      names(xx) <- NULL
      xx
    }) %>% unlist()
  }) %>% t() %>% as.data.table()
  setnames(features_conc, paste0(rep("V", ncol(features_conc)), seq(ncol(features_conc))))
  features_conc[, peptide := substr(sub("\\X.*", "", feature_list), 1, 15)]
  features_conc[, pseudosequence := substr(feature_list, nchar(feature_list) - 34 + 1, nchar(feature_list))]
  features_conc[, allele := allele_list]
  setcolorder(features_conc, c("peptide", "allele", "pseudosequence", colnames(features_conc)[1:(length(colnames(features_conc)) - 3)]))
  
  # Append BA feature
  BA_feature <- fread(BA_file)
  BA_feature <- BA_feature[, .(peptide, allele, BA)]
  names(BA_feature) <- c("peptide", "allele", "V1660")
  features_conc <- features_conc[BA_feature, on = .(peptide, allele)]
  
  fwrite(features_conc, file = filename, sep = ",")
}

featurization_BA(all_data, mhc_pseudoseq, "./cleaned_data/stab_41.csv", "./features/Stab_BLOSUM62_featuresBA.csv")

# Featurization function (based on BLOSUM62 PLUS Eluted Ligand predictions from netMHCpan4.1)
featurization_EL <- function(dataset, pseudosequences, EL_file, filename) {
  
  # Make the feature string comprising of both the extended peptide sequence and the MHC pseudosequence
  setkey(dataset, allele)
  setkey(pseudosequences, allele)
  dataset <- dataset[mhc_pseudoseq, nomatch=0]
  dataset[, features := paste0(peptide, pseudosequence)]
  
  feature_list <- dataset[, features] %>% unlist
  
  peptide_list <- feature_list %>% strsplit("")
  allele_list <- dataset[, allele]
  
  i <- 0
  features_conc <- sapply(peptide_list, function(aa_seq) {
    i <<- i + 1
    features <- lapply(aa_seq, function(amino_acid) {
      xx <- BLOSUM62[amino_acid, ]
      names(xx) <- NULL
      xx
    }) %>% unlist()
  }) %>% t() %>% as.data.table()
  setnames(features_conc, paste0(rep("V", ncol(features_conc)), seq(ncol(features_conc))))
  features_conc[, peptide := substr(sub("\\X.*", "", feature_list), 1, 15)]
  features_conc[, pseudosequence := substr(feature_list, nchar(feature_list) - 34 + 1, nchar(feature_list))]
  features_conc[, allele := allele_list]
  setcolorder(features_conc, c("peptide", "allele", "pseudosequence", colnames(features_conc)[1:(length(colnames(features_conc)) - 3)]))
  
  # Append BA feature
  EL_feature <- fread(EL_file)
  EL_feature <- EL_feature[, .(peptide, allele, EL)]
  names(EL_feature) <- c("peptide", "allele", "V1660")
  features_conc <- features_conc[EL_feature, on = .(peptide, allele)]
  
  fwrite(features_conc, file = filename, sep = ",")
}

featurization_EL(all_data, mhc_pseudoseq, "./cleaned_data/stab_41.csv", "./features/Stab_BLOSUM62_featuresEL.csv")