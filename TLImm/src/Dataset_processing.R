library(data.table)
library(caret)
library(stringr)
library(Biostrings)
library(Peptides)
library(magrittr)

reset <- 1 # Continue from here

mhc_pseudoseq <- fread("./data/updated_pseudosequences.csv")
allele_list <- mhc_pseudoseq[, allele] %>% unique

Immuno_dataset <- fread("./data/IEDB_data.csv")
Immuno_dataset <- Immuno_dataset[!str_detect(peptide, "X")]
Immuno_dataset[, `:=`(immunogenicity = ifelse(immunogenicity == "Negative", 0, 1), 
                      test = NULL, respond = NULL,
                      HLA = paste0(substring(HLA, 1, 8), ":", substring(HLA, 9, 10)))]
Immuno_dataset <- Immuno_dataset[HLA %in% allele_list]

Immuno_dataset_Bin <- Immuno_dataset[, .(peptide, HLA, immunogenicity)]
names(Immuno_dataset_Bin) <- c("peptide", "allele", "measurement_value")

Dengue_Bin <- fread("./data/Dengue_data_binary.csv")
Dengue_Bin[, measurement_value := 1] # All are positive actually

Immuno_dataset_Bin <- rbindlist(list(Immuno_dataset_Bin, Dengue_Bin))
fwrite(Immuno_dataset_Bin, "./cleaned_data_v4/Sample_3/Bin_train_data.csv")

Immuno_dataset_Cont <- Immuno_dataset[, .(peptide, HLA, potential)]
names(Immuno_dataset_Cont) <- c("peptide", "allele", "measurement_value")

Dengue_Cont <- fread("./data/Dengue_data.csv")
Immuno_dataset_Cont <- rbindlist(list(Immuno_dataset_Cont, Dengue_Cont))
fwrite(Immuno_dataset_Cont, "./cleaned_data_v4/Sample_3/Cont_train_data.csv")

## BEFORE THAT: The following tests the label distribution percentage that results in the fewest data point loss. 
if (reset == 0) {
  chosen_label_dataset <- copy(Immuno_dataset_Bin)
  chosen_label_dataset[, chosen := 0]
} else {
  chosen_label_dataset <- fread("./data/chosen_data.csv")
}
percentages <- seq(0.01,0.99,0.01)
perc_analysis <- sapply(percentages, function(percentage) {
  allele_population <- split(chosen_label_dataset, by="allele")
  sumsum <- sapply(allele_population, function(allele_pop) {
    #print(dim(allele_pop)[1])
    immunogenic <- allele_pop[measurement_value == 1]
    non_immunogenic <- allele_pop[measurement_value == 0]
    percentage_of_allele <- dim(immunogenic)[1]/dim(allele_pop)[1]
    #print(percentage_of_allele)
    if (percentage_of_allele < percentage) {
      #print(round(dim(immunogenic)[1]*((1 - percentage)/percentage)))
      #non_immunogenic <- non_immunogenic[sample(.N, round(dim(immunogenic)[1]*((1 - percentage)/percentage)))]
      non_immunogenic <- non_immunogenic[order(chosen)][1:(round(dim(immunogenic)[1]*((1 - percentage)/percentage)))]
    }
    else {
      #print(round(dim(non_immunogenic)[1]*(percentage/(1 - percentage))))
      #immunogenic <- immunogenic[sample(.N, round(dim(non_immunogenic)[1]*(percentage/(1 - percentage))))]
      immunogenic <- immunogenic[order(chosen)][1:(round(dim(non_immunogenic)[1]*(percentage/(1 - percentage))))]
    }
    filtered_dataset <- rbindlist(list(immunogenic, non_immunogenic))
    #print(allele_pop[1, HLA])
    #print(dim(filtered_dataset)[1])
    dim(filtered_dataset)[1]
  })
  print(percentage)
  sum(sumsum)
  #print ("Press [enter] to continue")
  #number <- scan(n=1)
})

perc_analysis <- data.table(perc_analysis)
perc_analysis[, perc := percentages*100]
names(perc_analysis) <- c("Data_points", "Percentage")
ggplot(perc_analysis, aes(x=Percentage, y=Data_points)) +
  geom_line() +
  theme_bw() +
  theme(axis.text.x = element_text(size = 16), plot.title = element_text(hjust = 0.5, size = 20, face = "bold"), 
        strip.text.x = element_text(size = 16, face = "italic"), strip.text.y = element_text(size = 16, face = "italic"),
        axis.title.y = element_text(size = rel(1.6)), axis.text.y = element_text(angle = 0, hjust = 1, size = 16), 
        axis.title.x = element_text(size = rel(1.6))) +
  ylab("#Data Points") +
  xlab("Percentage (%)")

# From the above analysis, while 3% of 1's gives the best value, it's the 60% of 1's that makes the most sense. Let's see how both look?

#3%:

allele_population <- split(chosen_label_dataset, by="allele")
percentage <- 0.6
sumsum <- lapply(allele_population, function(allele_pop) {
  immunogenic <- allele_pop[measurement_value == 1]
  non_immunogenic <- allele_pop[measurement_value == 0]
  percentage_of_allele <- dim(immunogenic)[1]/dim(allele_pop)[1]
  if (percentage_of_allele < percentage) {
    #non_immunogenic <- non_immunogenic[sample(.N, round(dim(immunogenic)[1]*((1 - percentage)/percentage)))]
    non_immunogenic <- non_immunogenic[order(chosen)][1:(round(dim(immunogenic)[1]*((1 - percentage)/percentage)))]
  }
  else {
    #immunogenic <- immunogenic[sample(.N, round(dim(non_immunogenic)[1]*(percentage/(1 - percentage))))]
    immunogenic <- immunogenic[order(chosen)][1:(round(dim(non_immunogenic)[1]*(percentage/(1 - percentage))))]
  }
  filtered_dataset <- rbindlist(list(immunogenic, non_immunogenic))
  filtered_dataset
})
Immuno_dataset_expanded <- rbindlist(sumsum)
Immuno_dataset_expanded[, chosen := 1]

chosen_label <- Immuno_dataset_expanded[chosen_label_dataset, on = c("peptide", "allele", "measurement_value")]
chosen_label <- chosen_label[, lapply(.SD, function(x){ ifelse(is.na(x), 0, x) })] # For now deprecated: 
chosen_label[, chosen := ifelse((chosen + i.chosen)/2 > 0, 1, 0)]
chosen_label[, i.chosen := NULL]
fwrite(chosen_label, "./data/chosen_data.csv")

# Subsample based on hits population: If I have more zeroes, remove (randomly?)
#Immunogenic_by_allele <- Immuno_dataset[immunogenicity == 1, .N, by = HLA]
#allele_population <- lapply(split(Immunogenic_by_allele, by="HLA"), unlist, use.names=FALSE)
#neg_list <- lapply(allele_population, function(allele_pop) {
#  allele_size <- dim(Immuno_dataset[immunogenicity == 0 & HLA == allele_pop[1]])[1]
#  Immuno_dataset[immunogenicity == 0 & HLA == allele_pop[1]][sample(.N, min(as.numeric(allele_pop[2]), allele_size))]
#})

#Non_Immuno <- rbindlist(neg_list)
#Immuno <- Immuno_dataset[immunogenicity == 1]
#Immuno_dataset_expanded <- rbindlist(list(Non_Immuno, Immuno))

# ...but if I have more ones, add (maybe experimental non-binders)
# I am commenting this out, as experimental non-binders could probably obfuscate things
# Let's instead remove randomly immunogenic peptides, even though that negates a lot of the dataset.

#ba_data <- fread("../data/mhcflurry_data.csv")
#ba_data <- ba_data[allele %in% allele_list]
#ba_data <- ba_data[measurement_inequality == '>' & measurement_value >= 500]
#ba_data <- ba_data[, `:=`(immunogenicity = 0, potential = 0, HLA = allele, allele = NULL)][, .(peptide, HLA, immunogenicity, potential)]
#ba_data <- ba_data[HLA %in% (Immuno_dataset_expanded[, HLA] %>% unlist)]

#Non_Immunogenic_by_allele <- Immuno_dataset_expanded[immunogenicity == 0, .N, by = HLA]
#allele_population <- lapply(split(Non_Immunogenic_by_allele, by="HLA"), unlist, use.names=FALSE)
#neg_list <- lapply(allele_population, function(allele_pop) {
#  allele_size <- dim(Immuno_dataset_expanded[immunogenicity == 1 & HLA == allele_pop[1]])[1]
#  ba_data[HLA == allele_pop[1]][sample(.N, max(allele_size - as.numeric(allele_pop[2]), 0))]
#})

#Extra_Non_Immuno <- rbindlist(neg_list)

#Immuno_dataset_expanded <- rbindlist(list(Immuno_dataset_expanded, Extra_Non_Immuno))

Immuno_dataset_expanded_Bin <- merge(Immuno_dataset_Bin, Immuno_dataset_expanded, on = c("peptide", "allele", "measurement_value"), nomatch=0)[, .(peptide, allele, measurement_value)]
names(Immuno_dataset_expanded_Bin) <- c("peptide", "allele", "measurement_value")
fwrite(Immuno_dataset_expanded_Bin, "./cleaned_data_v4/Sample_3/Bin_balanced_train_data.csv")
Immuno_dataset_expanded_Cont <- merge(Immuno_dataset_Cont, Immuno_dataset_expanded[, .(peptide, allele)])
names(Immuno_dataset_expanded_Cont) <- c("peptide", "allele", "measurement_value")
fwrite(Immuno_dataset_expanded_Cont, "./cleaned_data_v4/Sample_3/Cont_balanced_train_data.csv")

Immuno_dataset_Bin <- fread("./cleaned_data_v4/Sample_3/Bin_train_data.csv")
Immuno_dataset_Cont <- fread("./cleaned_data_v4/Sample_3/Cont_train_data.csv")
Immuno_dataset_expanded_Bin <- fread("./cleaned_data_v4/Sample_3/Bin_balanced_train_data.csv")
Immuno_dataset_expanded_Cont <- fread("./cleaned_data_v4/Sample_3/Cont_balanced_train_data.csv")
  
folds <- 5

cvIndex <- createFolds(Immuno_dataset_Bin$`measurement_value`, k = folds, list = FALSE, returnTrain = FALSE)
for (i in seq(folds)) {
  fwrite(Immuno_dataset_Bin[cvIndex == i], paste0("./folds_v4/Sample_3/Bin_folds/test_", i, ".csv"))
  fwrite(Immuno_dataset_Bin[cvIndex != i], paste0("./folds_v4/Sample_3/Bin_folds/train_", i, ".csv"))
}

cvIndex <- createFolds(Immuno_dataset_Cont$`measurement_value`, k = folds, list = FALSE, returnTrain = FALSE)
for (i in seq(folds)) {
  fwrite(Immuno_dataset_Cont[cvIndex == i], paste0("./folds_v4/Sample_3/Cont_folds/test_", i, ".csv"))
  fwrite(Immuno_dataset_Cont[cvIndex != i], paste0("./folds_v4/Sample_3/Cont_folds/train_", i, ".csv"))
}

cvIndex <- createFolds(Immuno_dataset_expanded_Bin$`measurement_value`, k = folds, list = FALSE, returnTrain = FALSE)
for (i in seq(folds)) {
  fwrite(Immuno_dataset_expanded_Bin[cvIndex == i], paste0("./folds_v4/Sample_3/Bin_balanced_folds/test_", i, ".csv"))
  fwrite(Immuno_dataset_expanded_Bin[cvIndex != i], paste0("./folds_v4/Sample_3/Bin_balanced_folds/train_", i, ".csv"))
}

cvIndex <- createFolds(Immuno_dataset_expanded_Cont$`measurement_value`, k = folds, list = FALSE, returnTrain = FALSE)
for (i in seq(folds)) {
  fwrite(Immuno_dataset_expanded_Cont[cvIndex == i], paste0("./folds_v4/Sample_3/Cont_balanced_folds/test_", i, ".csv"))
  fwrite(Immuno_dataset_expanded_Cont[cvIndex != i], paste0("./folds_v4/Sample_3/Cont_balanced_folds/train_", i, ".csv"))
}

for (i in seq(folds)) {
  temp_train_fold <- fread(paste0("./folds_v4/Sample_3/Bin_folds/train_", i, ".csv"))
  cvIndex <- createFolds(temp_train_fold$`measurement_value`, k = folds, list = FALSE, returnTrain = FALSE)
  for (i in seq(folds)) {
    fwrite(temp_train_fold[cvIndex != i], paste0("./folds_v4/Sample_3/Bin_folds/", i, "/train_", i, ".csv"))
    fwrite(temp_train_fold[cvIndex == i], paste0("./folds_v4/Sample_3/Bin_folds/", i, "/val_", i, ".csv"))
  }
}

for (i in seq(folds)) {
  temp_train_fold <- fread(paste0("./folds_v4/Sample_3/Cont_folds/train_", i, ".csv"))
  cvIndex <- createFolds(temp_train_fold$`measurement_value`, k = folds, list = FALSE, returnTrain = FALSE)
  for (i in seq(folds)) {
    fwrite(temp_train_fold[cvIndex != i], paste0("./folds_v4/Sample_3/Cont_folds/", i, "/train_", i, ".csv"))
    fwrite(temp_train_fold[cvIndex == i], paste0("./folds_v4/Sample_3/Cont_folds/", i, "/val_", i, ".csv"))
  }
}

for (i in seq(folds)) {
  temp_train_fold <- fread(paste0("./folds_v4/Sample_3/Bin_balanced_folds/train_", i, ".csv"))
  cvIndex <- createFolds(temp_train_fold$`measurement_value`, k = folds, list = FALSE, returnTrain = FALSE)
  for (i in seq(folds)) {
    fwrite(temp_train_fold[cvIndex != i], paste0("./folds_v4/Sample_3/Bin_balanced_folds/", i, "/train_", i, ".csv"))
    fwrite(temp_train_fold[cvIndex == i], paste0("./folds_v4/Sample_3/Bin_balanced_folds/", i, "/val_", i, ".csv"))
  }
}

for (i in seq(folds)) {
  temp_train_fold <- fread(paste0("./folds_v4/Sample_3/Cont_balanced_folds/train_", i, ".csv"))
  cvIndex <- createFolds(temp_train_fold$`measurement_value`, k = folds, list = FALSE, returnTrain = FALSE)
  for (i in seq(folds)) {
    fwrite(temp_train_fold[cvIndex != i], paste0("./folds_v4/Sample_3/Cont_balanced_folds/", i, "/train_", i, ".csv"))
    fwrite(temp_train_fold[cvIndex == i], paste0("./folds_v4/Sample_3/Cont_balanced_folds/", i, "/val_", i, ".csv"))
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

Immuno_dataset_Bin <- fread("./cleaned_data_v4/Sample_3/Bin_train_data.csv")
Immuno_dataset_Cont <- fread("./cleaned_data_v4/Sample_3/Cont_train_data.csv")
Immuno_dataset_expanded_Bin <- fread("./cleaned_data_v4/Sample_3/Bin_balanced_train_data.csv")
Immuno_dataset_expanded_Cont <- fread("./cleaned_data_v4/Sample_3/Cont_balanced_train_data.csv")

Immuno_dataset_Bin[, peptide := lapply(peptide, sequences_to_fixed_length_index_encoded_array)]
Immuno_dataset_Cont[, peptide := lapply(peptide, sequences_to_fixed_length_index_encoded_array)]
Immuno_dataset_expanded_Bin[, peptide := lapply(peptide, sequences_to_fixed_length_index_encoded_array)]
Immuno_dataset_expanded_Cont[, peptide := lapply(peptide, sequences_to_fixed_length_index_encoded_array)]

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

featurization(Immuno_dataset_Bin, mhc_pseudoseq, "./features_v4/Sample_3/Bin_BLOSUM62_features.csv")
featurization(Immuno_dataset_Cont, mhc_pseudoseq, "./features_v4/Sample_3/Cont_BLOSUM62_features.csv")
featurization(Immuno_dataset_expanded_Bin, mhc_pseudoseq, "./features_v4/Sample_3/Bin_balanced_BLOSUM62_features.csv")
featurization(Immuno_dataset_expanded_Cont, mhc_pseudoseq, "./features_v4/Sample_3/Cont_balanced_BLOSUM62_features.csv")

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
  features_conc <- features_conc[BA_feature, on = .(peptide, allele), nomatch = 0]
  
  fwrite(features_conc, file = filename, sep = ",")
}

featurization_BA(Immuno_dataset_Bin, mhc_pseudoseq, "./data/Bin_new_41.csv", "./features_v4/Sample_3/Bin_BLOSUM62_featuresBA.csv")
featurization_BA(Immuno_dataset_Cont, mhc_pseudoseq, "./data/Bin_new_41.csv", "./features_v4/Sample_3/Cont_BLOSUM62_featuresBA.csv")
featurization_BA(Immuno_dataset_expanded_Bin, mhc_pseudoseq, "./data/Bin_new_41.csv", "./features_v4/Sample_3/Bin_balanced_BLOSUM62_featuresBA.csv")
featurization_BA(Immuno_dataset_expanded_Cont, mhc_pseudoseq, "./data/Bin_new_41.csv", "./features_v4/Sample_3/Cont_balanced_BLOSUM62_featuresBA.csv")

# Featurization function (based on BLOSUM62 PLUS Binding Affinity predictions from netMHCpan4.1)
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
  features_conc <- features_conc[EL_feature, on = .(peptide, allele), nomatch = 0]
  
  fwrite(features_conc, file = filename, sep = ",")
}

featurization_EL(Immuno_dataset_Bin, mhc_pseudoseq, "./data/Bin_new_41.csv", "./features_v4/Sample_3/Bin_BLOSUM62_featuresEL.csv")
featurization_EL(Immuno_dataset_Cont, mhc_pseudoseq, "./data/Bin_new_41.csv", "./features_v4/Sample_3/Cont_BLOSUM62_featuresEL.csv")
featurization_EL(Immuno_dataset_expanded_Bin, mhc_pseudoseq, "./data/Bin_new_41.csv", "./features_v4/Sample_3/Bin_balanced_BLOSUM62_featuresEL.csv")
featurization_EL(Immuno_dataset_expanded_Cont, mhc_pseudoseq, "./data/Bin_new_41.csv", "./features_v4/Sample_3/Cont_balanced_BLOSUM62_featuresEL.csv")