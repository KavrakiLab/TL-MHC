# Importing libraries
library(data.table)
library(Peptides)
library(stringr)
library(magrittr)
library(Biostrings)
library(caret)

# Load the data
data <- fread("../data/mhcflurry_data.csv")

# Load pseudosequences
# The alleles in the dataset NEED to have pseudosequences!
mhc_pseudoseq <- fread("../data/updated_pseudosequences.csv")
allele_list <- mhc_pseudoseq[, allele] %>% unique
data <- data[allele %in% allele_list]

# For Mass_spec data we will use binary cross-entropy
# For Binding affinity data, we will use MHCflurry loss
# So it's like a combination of the 2 things
data[, measurement_kind := ifelse(measurement_kind == "mass_spec", 1, 0)]
data[, measurement_inequality := ifelse(measurement_kind == 1,
                                        ifelse(measurement_inequality == "<", 1, 0),
                                        ifelse(measurement_inequality == "<", 1, 
                                               ifelse(measurement_inequality == "=", 0, -1)))]

# Throw anything except things we will need
data[, `:=`(original_allele = NULL, measurement_source = NULL, measurement_type = NULL,
            fold_0 = NULL, fold_1 = NULL, fold_2 = NULL, fold_3 = NULL)]

# Load test set for balancing the mass_spec data in terms of hits and misses + extra data to train!
MONOALLELIC_data <- fread("../data/MONOALLELIC_dataset.csv")
MONOALLELIC_data[, `:=`(protein_accession = NULL, sample_id = NULL, sample_group = NULL, 
                        n_flank = NULL, c_flank = NULL, mixmhcpred = NULL, netmhcpan4.ba = NULL, 
                        netmhcpan4.el = NULL, mhcflurry_variant.ba =  NULL)]
MONOALLELIC_data <- unique(MONOALLELIC_data, by = c("peptide", "hla", "hit"))
names(MONOALLELIC_data) <- c("peptide", "measurement_value", "allele")
MONOALLELIC_data[, `:=`(measurement_inequality = measurement_value, 
                        measurement_kind = 1)]
MONOALLELIC_data <- MONOALLELIC_data[allele %in% allele_list]
pos_MONOALLELIC <- MONOALLELIC_data[measurement_value == 1]
neg_MONOALLELIC <- MONOALLELIC_data[measurement_value == 0]

# Duplicates routine:
mass_spec <- data[measurement_kind == 1]
mass_spec[, measurement_value := measurement_inequality]
#mass_spec <- list(mass_spec, pos_MONOALLELIC) %>% rbindlist(use.names = TRUE)

# Subsample based on hits population
mass_spec_positives_by_allele <- mass_spec[measurement_value == 1, .N, by = allele]
allele_population <- lapply(split(mass_spec_positives_by_allele, by="allele"), unlist, use.names=FALSE)
neg_list <- lapply(allele_population, function(allele_pop) {
  print(allele_pop[1])
  print(allele_pop[2])
  allele_size <- dim(neg_MONOALLELIC[allele == allele_pop[1]])[1]
  neg_MONOALLELIC[allele == allele_pop[1]][sample(.N, min(as.numeric(allele_pop[2]), allele_size))]
})
  
neg_mass_spec <- rbindlist(neg_list)
mass_spec <- rbindlist(list(neg_mass_spec, mass_spec), use.names=TRUE)

affinity_equal <- data[measurement_kind == 0 & measurement_inequality == 0]
affinity_non_equal <- data[measurement_kind == 0 & measurement_inequality != 0]

# 1) MASS_SPEC: 
# If there are duplicates with the same measurement_inequality, remove duplicates
mass_spec <- unique(mass_spec, by = c("allele", "peptide", "measurement_inequality"))
# If they are not, remove them all, as they are ambiguous entries!
mass_spec <- mass_spec[!((duplicated(mass_spec, by = c("allele", "peptide"))
                        | duplicated(mass_spec, fromLast = TRUE, by = c("allele", "peptide"))))]

# 2) BA:
# 2a) When we have only equalities, average the results
affinity_equal <- affinity_equal[, lapply(.SD, mean), 
                                   by = .(allele, peptide, measurement_kind, measurement_inequality), 
                                   .SDcols = c("measurement_value")]

# 2b) When we have inequalities, it depends on the inequalities:
#    - If only smaller than , take the smallest value recorded
#    - If bigger than, take the largest value recorder
#    - If both inequalitities exist, probably it's better to discard those data points as ambiguous
affinity_non_equal[, inequalitites := paste(measurement_inequality, collapse = ","), by = .(allele, peptide)]
affinity_non_equal[, inequalitites_unique := paste(unique(strsplit(inequalitites, ",")[[1]]), collapse = ","), by = .(allele, peptide)]
bigger_thans <- affinity_non_equal[inequalitites_unique == '-1', .(max(measurement_value)), by = .(allele, peptide, inequalitites_unique)]
smaller_thans <- affinity_non_equal[inequalitites_unique == '1', .(min(measurement_value)), by = .(allele, peptide, inequalitites_unique)]
mixed <- affinity_non_equal[inequalitites_unique == '-1,1' | inequalitites_unique == '1,-1'] # Not included
affinity_non_equal <- list(bigger_thans, smaller_thans) %>% rbindlist()
affinity_non_equal[, `:=`(measurement_inequality = inequalitites_unique,
                         inequalitites_unique = NULL,
                         measurement_value = V1,
                         V1 = NULL,
                         measurement_kind = 0)]

# Unify equalities with inequalities on affinity data
# If there are duplicates, keep the more specific, equality counterpart
filtered_affinity <- list(affinity_equal, affinity_non_equal) %>% rbindlist(use.names=TRUE)
dups <- filtered_affinity[((duplicated(filtered_affinity, by = c("allele", "peptide"))
                          | duplicated(filtered_affinity, fromLast = TRUE, by = c("allele", "peptide"))))]
dups <- dups[measurement_inequality == 0]
filtered_affinity <- filtered_affinity[!((duplicated(filtered_affinity, by = c("allele", "peptide"))
                                        | duplicated(filtered_affinity, fromLast = TRUE, by = c("allele", "peptide"))))]
filtered_affinity <- list(filtered_affinity, dups) %>% rbindlist(use.names=TRUE)

# Remove datapoints that have bindin affinity equal to zero before the transformation:
filtered_affinity <- filtered_affinity[measurement_value != 0]

filtered_affinity[, measurement_value := ifelse(measurement_value <= 50000, measurement_value, 50000)]
filtered_affinity[, measurement_value := 1 - log(measurement_value, base = 50000)]

# Unify mass_spec and binding affinities again.
# Filter the data points that are duplicates, and they are contradictory (mass_spec == 1 / ba == -1)
# A not-so-good but initial strategy is to also filter out the mass_spec data if I have BA data
data <- list(mass_spec, filtered_affinity) %>% rbindlist(use.names=TRUE)

dups <- data[((duplicated(data, by = c("allele", "peptide"))
             | duplicated(data, fromLast = TRUE, by = c("allele", "peptide"))))]
dups <- dups[measurement_kind == 0]
data <- data[!((duplicated(data, by = c("allele", "peptide"))
              | duplicated(data, fromLast = TRUE, by = c("allele", "peptide"))))]
data <- list(data, dups) %>% rbindlist(use.names=TRUE)
setcolorder(data, c("peptide", "allele", "measurement_kind", "measurement_inequality", "measurement_value"))
fwrite(data, "../cleaned_data/cleaned_train_data.csv")

# Make the ouside validation set
MONOALLELIC_copy <- copy(MONOALLELIC_data) 
MONOALLELIC_copy <- fsetdiff(MONOALLELIC_copy[, .(allele, peptide, measurement_value)], data[, .(allele, peptide, measurement_value)])
pos_MONOALLELIC <- MONOALLELIC_copy[measurement_value == 1]
neg_MONOALLELIC <- MONOALLELIC_copy[measurement_value == 0][sample(.N, 200000)]
MONOALLELIC_copy <- rbindlist(list(pos_MONOALLELIC, neg_MONOALLELIC))
MONOALLELIC_copy[, `:=`(measurement_inequality = measurement_value, 
                        measurement_kind = 1)]

# Making the folds
folds <- 5
cvIndex <- createFolds(data$allele, k = folds, list = FALSE, returnTrain = FALSE)

for (i in seq(folds)) {
  fwrite(data[cvIndex == i], paste0("../folds/val_", i, ".csv"))
  fwrite(data[cvIndex != i], paste0("../folds/train_", i, ".csv"))
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

data[, peptide := lapply(peptide, sequences_to_fixed_length_index_encoded_array)]

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

rm(filtered_affinity)
rm(mass_spec)
rm(neg_mass_spec)
rm(neg_list)
rm(bigger_thans)
rm(smaller_thans)
rm(affinity_non_equal)
rm(affinity_equal)
rm(MONOALLELIC_data)
rm(pos_MONOALLELIC)
rm(neg_MONOALLELIC)
gc()
gc()

# This is done per allele, as the memory is huge
per_allele_list <- split(data, by=c("allele"))
print("Featurization")
res <- lapply(per_allele_list, function(allele_dt) {
  allele_name <- allele_dt[, allele] %>% unique
  print(allele_name)
  featurization(allele_dt, mhc_pseudoseq, paste0("../features/", allele_name, "_Train_BLOSUM62_extended.csv"))
})

# This is done through system calls, I cannot load the whole thing it seems!
system("head -n 1 ../features/'SLA-2*05:02_Train_BLOSUM62_extended.csv' > ../features/all_train_features.csv")
system("tail -q -n +2 ../features/*.csv >> ../features/all_train_features.csv")
system("rm ../features/*_extended.csv")

# Finish by making the validation set
Other <- fread("../data/test_peptides.csv")
Other[, `:=`( peptide = NULL, classification_label = NULL)]
Other[, peptide := Sequence]
Other[, ba_result := as.numeric(ba_result)]

Other[, measurement_value := 1 - log(ba_result, base = 50000)]
Other[, measurement_kind := 0]
Other[, measurement_inequality := 0]
Other[, `:=`(NetMHC4.0 = NULL, `NetMHCpan 4.0` = NULL, MHCFlurry2.0 = NULL, ba_result = NULL, Sequence = NULL)]
Other <- rbindlist(list(Other, MONOALLELIC_copy), use.names=TRUE)
fwrite(Other, "../cleaned_data/test_data.csv")

Other[, peptide := lapply(peptide, sequences_to_fixed_length_index_encoded_array)]
featurization(Other, mhc_pseudoseq, paste0("../features/Test_BLOSUM62_extended.csv"))
