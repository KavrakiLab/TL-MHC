library(data.table)
library(Peptides)
library(stringr)
library(magrittr)
library(Biostrings)

`%notin%` <- Negate(`%in%`)

# Load the data
data <- fread("./data/MHCflurry_data/mhcflurry_data.csv")

# Load pseudosequences
# The alleles in the dataset NEED to have pseudosequences!
mhc_pseudoseq <- fread("./data/updated_pseudosequences.csv")
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
#MONOALLELIC_data <- fread("./data/MHCflurry_data/MONOALLELIC_dataset.csv")
#MONOALLELIC_data[, `:=`(protein_accession = NULL, sample_id = NULL, sample_group = NULL, 
#                        n_flank = NULL, c_flank = NULL, mixmhcpred = NULL, netmhcpan4.ba = NULL, 
#                        netmhcpan4.el = NULL, mhcflurry_variant.ba =  NULL)]
#MONOALLELIC_data <- unique(MONOALLELIC_data, by = c("peptide", "hla", "hit"))
#names(MONOALLELIC_data) <- c("peptide", "measurement_value", "allele")
#MONOALLELIC_data[, `:=`(measurement_inequality = measurement_value, 
#                        measurement_kind = 1)]
#MONOALLELIC_data <- MONOALLELIC_data[allele %in% allele_list]
#pos_MONOALLELIC <- MONOALLELIC_data[measurement_value == 1]
#neg_MONOALLELIC <- MONOALLELIC_data[measurement_value == 0]

# Duplicates routine:
mass_spec <- data[measurement_kind == 1]
mass_spec[, measurement_value := measurement_inequality]
mass_spec <- mass_spec[measurement_value == 1]

# Subsample based on hits population
#mass_spec_positives_by_allele <- mass_spec[measurement_value == 1, .N, by = allele]
#allele_population <- lapply(split(mass_spec_positives_by_allele, by="allele"), unlist, use.names=FALSE)
#neg_list <- lapply(allele_population, function(allele_pop) {
#  print(allele_pop[1])
#  print(allele_pop[2])
#  allele_size <- dim(neg_MONOALLELIC[allele == allele_pop[1]])[1]
#  neg_MONOALLELIC[allele == allele_pop[1]][sample(.N, min(as.numeric(allele_pop[2]), allele_size))]
#})
#neg_mass_spec <- rbindlist(neg_list)
#mass_spec <- rbindlist(list(neg_mass_spec, mass_spec), use.names=TRUE)

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
mixed <- affinity_non_equal[inequalitites_unique == '-1,1' | inequalitites_unique == '1,-1'] # Total: 337, not included
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

# MONOALLELIC MHCFlurry
MONOALLELIC_data <- fread("./data/MHCflurry_data/MONOALLELIC_dataset.csv")
MONOALLELIC_data[, `:=`(protein_accession = NULL, sample_id = NULL, sample_group = NULL, 
                        n_flank = NULL, c_flank = NULL, mixmhcpred = NULL, netmhcpan4.ba = NULL, 
                        netmhcpan4.el = NULL, mhcflurry_variant.ba =  NULL)]
MONOALLELIC_data <- unique(MONOALLELIC_data, by = c("peptide", "hla", "hit"))
names(MONOALLELIC_data) <- c("peptide", "measurement_value", "allele")
MONOALLELIC_data[, `:=`(measurement_inequality = measurement_value, 
                        measurement_kind = 1)]
MONOALLELIC_data <- MONOALLELIC_data[allele %in% allele_list]
pos_MONOALLELIC <- MONOALLELIC_data[measurement_value == 1]
mass_spec <- rbindlist(list(mass_spec, pos_MONOALLELIC), use.names=TRUE)
mass_spec <- unique(mass_spec)
rm(MONOALLELIC_data)
rm(pos_MONOALLELIC)
gc()

# MULTIALLELIC OLD
MULTIALLELIC_OLD <- fread("./data/MHCflurry_data/Deconvoluted_old.csv")
MULTIALLELIC_OLD <- MULTIALLELIC_OLD[, .(Peptide, MHC, hit)]
names(MULTIALLELIC_OLD) <- c("peptide", "allele", "measurement_value")
MULTIALLELIC_OLD[, `:=`(measurement_inequality = measurement_value, 
                        measurement_kind = 1)]
MULTIALLELIC_OLD <- MULTIALLELIC_OLD[allele %in% allele_list]
MULTIALLELIC_OLD <- MULTIALLELIC_OLD[measurement_value == 1]
mass_spec <- rbindlist(list(mass_spec, MULTIALLELIC_OLD), use.names=TRUE)
mass_spec <- unique(mass_spec)
rm(MULTIALLELIC_OLD)
gc()

# MULTIALLELIC RECENT
MULTIALLELIC_RECENT <- fread("./data/MHCflurry_data/Deconvoluted_recent.csv")
MULTIALLELIC_RECENT <- MULTIALLELIC_RECENT[, .(Peptide, MHC, hit)]
names(MULTIALLELIC_RECENT) <- c("peptide", "allele", "measurement_value")
MULTIALLELIC_RECENT[, `:=`(measurement_inequality = measurement_value, 
                        measurement_kind = 1)]
MULTIALLELIC_RECENT <- MULTIALLELIC_RECENT[allele %in% allele_list]
MULTIALLELIC_RECENT <- MULTIALLELIC_RECENT[measurement_value == 1]
mass_spec <- rbindlist(list(mass_spec, MULTIALLELIC_RECENT), use.names=TRUE)
mass_spec <- unique(mass_spec)
rm(MULTIALLELIC_RECENT)
gc()

mass_spec_positives_by_allele <- mass_spec[measurement_value == 1, .N, by = allele]

rm(affinity_equal)
rm(affinity_non_equal)
rm(bigger_thans)
rm(data)
rm(dups)
rm(mixed)
rm(smaller_thans)
gc()

# SAMPLE NEGATIVES!
MULTIALLELIC_OLD <- fread("./data/MHCflurry_data/Deconvoluted_old.csv")
MULTIALLELIC_OLD <- MULTIALLELIC_OLD[, .(Peptide, MHC, hit)]
names(MULTIALLELIC_OLD) <- c("peptide", "allele", "measurement_value")
MULTIALLELIC_OLD[, `:=`(measurement_inequality = measurement_value, 
                        measurement_kind = 1)]
MULTIALLELIC_OLD <- MULTIALLELIC_OLD[allele %in% allele_list]
MULTIALLELIC_OLD <- MULTIALLELIC_OLD[measurement_value == 0]

mass_spec_positives_by_allele <- mass_spec[measurement_value == 1, .N, by = allele]
allele_population <- lapply(split(mass_spec_positives_by_allele, by="allele"), unlist, use.names=FALSE)
neg_list <- lapply(allele_population, function(allele_pop) {
  allele_size <- dim(MULTIALLELIC_OLD[allele == allele_pop[1]])[1]
  if (allele_size >= as.numeric(allele_pop[2])*3) {
    print(allele_pop[1])
    print(allele_pop[2])
    MULTIALLELIC_OLD[allele == allele_pop[1]][sample(.N, as.numeric(allele_pop[2])*3)]
  }
})
neg_mass_spec <- rbindlist(neg_list)
mass_spec <- rbindlist(list(neg_mass_spec, mass_spec), use.names=TRUE)

rm(MULTIALLELIC_OLD)
rm(neg_list)
rm(neg_mass_spec)
gc()


MULTIALLELIC_RECENT <- fread("./data/MHCflurry_data/Deconvoluted_recent.csv")
MULTIALLELIC_RECENT <- MULTIALLELIC_RECENT[, .(Peptide, MHC, hit)]
names(MULTIALLELIC_RECENT) <- c("peptide", "allele", "measurement_value")
MULTIALLELIC_RECENT[, `:=`(measurement_inequality = measurement_value, 
                        measurement_kind = 1)]
MULTIALLELIC_RECENT <- MULTIALLELIC_RECENT[allele %in% allele_list]
MULTIALLELIC_RECENT <- MULTIALLELIC_RECENT[measurement_value == 0]

filter_those_out <- mass_spec[measurement_value == 0, allele] %>% unique()
mass_spec_positives_by_allele <- mass_spec[measurement_value == 1 & allele %notin% filter_those_out, .N, by = allele]
allele_population <- lapply(split(mass_spec_positives_by_allele, by="allele"), unlist, use.names=FALSE)
neg_list <- lapply(allele_population, function(allele_pop) {
  allele_size <- dim(MULTIALLELIC_RECENT[allele == allele_pop[1]])[1]
  if (allele_size >= as.numeric(allele_pop[2])*3) {
    print(allele_pop[1])
    print(allele_pop[2])
    MULTIALLELIC_RECENT[allele == allele_pop[1]][sample(.N, as.numeric(allele_pop[2])*3)]
  }
})
neg_mass_spec <- rbindlist(neg_list)
mass_spec <- rbindlist(list(neg_mass_spec, mass_spec), use.names=TRUE)

rm(MULTIALLELIC_RECENT)
rm(neg_list)
rm(neg_mass_spec)
gc()


MONOALLELIC_data <- fread("./data/MHCflurry_data/MONOALLELIC_dataset.csv")
MONOALLELIC_data[, `:=`(protein_accession = NULL, sample_id = NULL, sample_group = NULL, 
                        n_flank = NULL, c_flank = NULL, mixmhcpred = NULL, netmhcpan4.ba = NULL, 
                        netmhcpan4.el = NULL, mhcflurry_variant.ba =  NULL)]
MONOALLELIC_data <- unique(MONOALLELIC_data, by = c("peptide", "hla", "hit"))
names(MONOALLELIC_data) <- c("peptide", "measurement_value", "allele")
MONOALLELIC_data[, `:=`(measurement_inequality = measurement_value, 
                        measurement_kind = 1)]
MONOALLELIC_data <- MONOALLELIC_data[allele %in% allele_list]
MONOALLELIC_data <- MONOALLELIC_data[measurement_value == 0]

filter_those_out <- mass_spec[measurement_value == 0, allele] %>% unique()
mass_spec_positives_by_allele <- mass_spec[measurement_value == 1 & allele %notin% filter_those_out, .N, by = allele]
allele_population <- lapply(split(mass_spec_positives_by_allele, by="allele"), unlist, use.names=FALSE)
neg_list <- lapply(allele_population, function(allele_pop) {
  allele_size <- dim(MONOALLELIC_data[allele == allele_pop[1]])[1]
  if (allele_size >= as.numeric(allele_pop[2])*3) {
    print(allele_pop[1])
    print(allele_pop[2])
    MONOALLELIC_data[allele == allele_pop[1]][sample(.N, as.numeric(allele_pop[2])*3)]
  }
})
neg_mass_spec <- rbindlist(neg_list)
mass_spec <- rbindlist(list(neg_mass_spec, mass_spec), use.names=TRUE)

rm(MONOALLELIC_data)
rm(neg_list)
rm(neg_mass_spec)
gc()

human_proteome <- readAAStringSet("./data/MHCflurry_data/uniprot-compressed_true_download_true_format_fasta_query__28_28prote-2022.11.20-21.11.23.09.fasta")
protein_sequences <- as.data.frame(human_proteome)$x
sequence_length <- c(8, 9, 10, 11)
i <- 0 
negatives_list <- lapply(protein_sequences, function(sequence) {
  len <- sample(sequence_length, 1)
  peptide_list <- strsplit(sequence, paste0("(?<=.{", len , "})"), perl = TRUE)[[1]]
  peptide_list <- peptide_list[!grepl(c("U|X"), peptide_list)]
  peptide_list <- peptide_list[nchar(peptide_list)==max(nchar(peptide_list))]
  i <<- i + 1
  message(paste(i, "/", length(protein_sequences)),"\r",appendLF=FALSE)
  peptide_list
})
negatives_list <- unique(unlist(negatives_list))

filter_those_out <- mass_spec[measurement_value == 0, allele] %>% unique()
mass_spec_positives_by_allele <- mass_spec[measurement_value == 1 & allele %notin% filter_those_out, .N, by = allele]
allele_population <- lapply(split(mass_spec_positives_by_allele, by="allele"), unlist, use.names=FALSE)
neg_list <- lapply(allele_population, function(allele_pop) {
  print(allele_pop[1])
  print(allele_pop[2])
  data.table(allele = rep(allele_pop[1], as.numeric(allele_pop[2])*3), 
             peptide = sample(negatives_list, as.numeric(allele_pop[2])*3),
             measurement_value = rep(0, as.numeric(allele_pop[2])*3),
             measurement_inequality = rep(0, as.numeric(allele_pop[2])*3),
             measurement_kind = rep(1, as.numeric(allele_pop[2])*3))
})
neg_mass_spec <- rbindlist(neg_list)
mass_spec <- rbindlist(list(neg_mass_spec, mass_spec), use.names=TRUE)

sample <- sample.int(n = nrow(filtered_affinity), size = floor(.85*nrow(filtered_affinity)), replace = F)
filtered_affinity_train <- filtered_affinity[sample(.N, dim(filtered_affinity)[1]*0.85)]
filtered_affinity_test  <- fsetdiff(filtered_affinity, filtered_affinity_train)

rm(human_proteome)
rm(protein_sequences)
rm(sequence_length)
rm(negatives_list)
rm(filter_those_out)
rm(i)
rm(mass_spec_positives_by_allele)
rm(allele_population)
rm(neg_list)
rm(neg_mass_spec)
rm(filtered_affinity)
gc()

data <- list(mass_spec, filtered_affinity_train) %>% rbindlist(use.names=TRUE)
dups <- data[((duplicated(data, by = c("allele", "peptide"))
               | duplicated(data, fromLast = TRUE, by = c("allele", "peptide"))))]
dups <- dups[measurement_kind == 0]
data <- data[!((duplicated(data, by = c("allele", "peptide"))
                | duplicated(data, fromLast = TRUE, by = c("allele", "peptide"))))]
data <- list(data, dups) %>% rbindlist(use.names=TRUE)
setcolorder(data, c("peptide", "allele", "measurement_kind", "measurement_inequality", "measurement_value"))

rm(mass_spec)
rm(filtered_affinity_train)
rm(dups)
gc()

fwrite(data, "./data/MHCflurry_data/train_data_v3.csv")

# Make the validation set
MONOALLELIC_Pyke <- fread("./data/MHCflurry_data/MONOALLELIC_Pyke.csv")
MONOALLELIC_Pyke <- MONOALLELIC_Pyke[, .(peptide, allele, label)]
MONOALLELIC_Pyke[, allele := paste0(substr(allele, 1, 5), '*', substr(allele, 6, 10))]
names(MONOALLELIC_Pyke) <- c("peptide", "allele", "measurement_value")
MONOALLELIC_Pyke[, `:=`(measurement_inequality = measurement_value, 
                        measurement_kind = 1)]
MONOALLELIC_Pyke <- MONOALLELIC_Pyke[allele %in% allele_list]
MONOALLELIC_Pyke <- fsetdiff(MONOALLELIC_Pyke[, .(allele, peptide, measurement_value)], data[, .(allele, peptide, measurement_value)])
mass_spec_positives_by_allele <- MONOALLELIC_Pyke[, .N, by = allele]

human_proteome <- readAAStringSet("./data/MHCflurry_data/uniprot-compressed_true_download_true_format_fasta_query__28_28prote-2022.11.20-21.11.23.09.fasta")
protein_sequences <- as.data.frame(human_proteome)$x
sequence_length <- c(8, 9, 10, 11)
i <- 0 
negatives_list <- lapply(protein_sequences, function(sequence) {
  len <- sample(sequence_length, 1)
  peptide_list <- strsplit(sequence, paste0("(?<=.{", len , "})"), perl = TRUE)[[1]]
  peptide_list <- peptide_list[!grepl(c("U|X"), peptide_list)]
  peptide_list <- peptide_list[nchar(peptide_list)==max(nchar(peptide_list))]
  i <<- i + 1
  message(paste(i, "/", length(protein_sequences)),"\r",appendLF=FALSE)
  peptide_list
})
negatives_list <- unique(unlist(negatives_list))

allele_population <- lapply(split(mass_spec_positives_by_allele, by="allele"), unlist, use.names=FALSE)
neg_list <- lapply(allele_population, function(allele_pop) {
  print(allele_pop[1])
  print(allele_pop[2])
  data.table(allele = rep(allele_pop[1], as.numeric(allele_pop[2])*10), 
             peptide = sample(negatives_list, as.numeric(allele_pop[2])*10),
             measurement_value = rep(0, as.numeric(allele_pop[2])*10),
             measurement_inequality = rep(0, as.numeric(allele_pop[2])*10),
             measurement_kind = rep(1, as.numeric(allele_pop[2])*10))
})
neg_mass_spec <- rbindlist(neg_list)
MONOALLELIC_Pyke[, `:=`(measurement_inequality = measurement_value, 
                        measurement_kind = 1)]
MONOALLELIC_Pyke <- rbindlist(list(neg_mass_spec, MONOALLELIC_Pyke), use.names=TRUE)

Other <- fread("./data/MHCflurry_data/test_peptides.csv")
Other[, `:=`( peptide = NULL, classification_label = NULL)]
Other[, peptide := Sequence]
Other[, ba_result := as.numeric(ba_result)]

Other[, measurement_value := 1 - log(ba_result, base = 50000)]
Other[, measurement_kind := 0]
Other[, measurement_inequality := 0]
#Other[, NetMHC4.0 := 1 - log(NetMHC4.0, base = 50000)]
#Other[, `NetMHCpan 4.0` := 1 - log(`NetMHCpan 4.0`, base = 50000)]
#Other[, MHCFlurry2.0 := 1 - log(MHCFlurry2.0, base = 50000)]
Other[, `:=`(NetMHC4.0 = NULL, `NetMHCpan 4.0` = NULL, MHCFlurry2.0 = NULL, ba_result = NULL, Sequence = NULL)]
Other <- rbindlist(list(Other, filtered_affinity_test, MONOALLELIC_Pyke), use.names=TRUE)

fwrite(Other, "./data/MHCflurry_data/validation_data_v3.csv")

rm(MONOALLELIC_Pyke)
rm(data)
rm(Other)
rm(mass_spec_positives_by_allele)
rm(allele_population)
rm(neg_list)
rm(human_proteome)
rm(negatives_list)
rm(neg_mass_spec)
rm(protein_sequences)
rm(filtered_affinity_test)
gc()

# Make test set too!
MULTIALLELIC_Pyke <- fread("./data/MHCflurry_data/Pyke_multiallelic.csv")
MULTIALLELIC_Pyke <- MULTIALLELIC_Pyke[modification == "na", .(peptide, source)]
MULTIALLELIC_Pyke2 <- fread("./data/MHCflurry_data/Pyke_multiallelic_2.csv")
MULTIALLELIC_Pyke2[, allele := paste0(substr(allele, 1, 5), '*', substr(allele, 6, 10))]
MULTIALLELIC_Pyke2 <- MULTIALLELIC_Pyke2[MULTIALLELIC_Pyke, on = c("peptide", "source"), nomatch=0]
MULTIALLELIC_Pyke2 <- MULTIALLELIC_Pyke2[order(peptide)]
reps <- MULTIALLELIC_Pyke2[, .N, by = .(peptide, source)][, N]
MULTIALLELIC_Pyke2[, Indexes := rep(seq(1:length(reps)), reps)]
MULTIALLELIC_Pyke2[, measurement_value := 1]
MULTIALLELIC_Pyke2[, measurement_inequality := 1]
MULTIALLELIC_Pyke2[, measurement_kind := 1]
MULTIALLELIC_Pyke2[, source := NULL]

human_proteome <- readAAStringSet("./data/MHCflurry_data/uniprot-compressed_true_download_true_format_fasta_query__28_28prote-2022.11.20-21.11.23.09.fasta")
protein_sequences <- as.data.frame(human_proteome)$x
sequence_length <- c(8, 9, 10, 11)
i <- 0 
negatives_list <- lapply(protein_sequences, function(sequence) {
  len <- sample(sequence_length, 1)
  peptide_list <- strsplit(sequence, paste0("(?<=.{", len , "})"), perl = TRUE)[[1]]
  peptide_list <- peptide_list[!grepl(c("U|X"), peptide_list)]
  peptide_list <- peptide_list[nchar(peptide_list)==max(nchar(peptide_list))]
  i <<- i + 1
  message(paste(i, "/", length(protein_sequences)),"\r",appendLF=FALSE)
  peptide_list
})
negatives_list <- unique(unlist(negatives_list))

i <- 1
neg_list <- lapply(seq(1:length(reps)), function(Index) {
  message(paste(i, "/", length(reps)),"\r",appendLF=FALSE)
  dt <- MULTIALLELIC_Pyke2[Indexes == Index]
  i <<- i + 1
  data.table(allele = rep(dt[, allele], 3), 
             peptide = rep(sample(negatives_list, 3), times = 1, each = dim(dt)[1]),
             Indexes = rep(Index, dim(dt)[1]*3),
             measurement_value = rep(0, dim(dt)[1]*3),
             measurement_inequality = rep(0, dim(dt)[1]*3),
             measurement_kind = rep(1, dim(dt)[1]*3))
})
neg_mass_spec <- rbindlist(neg_list)
MULTIALLELIC_Pyke2 <- rbindlist(list(neg_mass_spec, MULTIALLELIC_Pyke2), use.names=TRUE)

fwrite(MULTIALLELIC_Pyke2, "./data/MHCflurry_data/test_data_v3.csv")

rm(MULTIALLELIC_Pyke)
rm(MULTIALLELIC_Pyke2)
rm(neg_list)
rm(human_proteome)
rm(negatives_list)
rm(neg_mass_spec)
rm(protein_sequences)
rm(reps)
rm(i)
gc()

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
  c(0)
}

# Load BLOSUM62 matrix
data(BLOSUM62)
BLOSUM62 <- BLOSUM62[-c(21, 22, 23, 25), -c(21, 22, 23, 25)]
BLOSUM62 <- apply(BLOSUM62, MARGIN = 2, FUN = function(X) (X - min(X))/diff(range(X)))
BLOSUM62['X', ] <- rep.int(0, 21)
BLOSUM62[, 'X'] <- rep.int(0, 21)
BLOSUM62['X', 'X'] <- 1

data <- fread("./data/MHCflurry_data/train_data_v3.csv")

data[, peptide := lapply(peptide, sequences_to_fixed_length_index_encoded_array)]
alleles_already_done <- list.files("./data/MHCFlurry_features/")
alleles_already_done <- sapply(strsplit(alleles_already_done, '_'), `[`, 1)
alleles_already_done <- append(alleles_already_done, "HLA-B*27:05")
data <- data[allele %notin% alleles_already_done]

# This is done per allele, as the memory is huge
per_allele_list <- split(data, by=c("allele"))
rm(data)
gc()
print("Featurization")
res <- lapply(per_allele_list, function(allele_dt) {
  allele_name <- allele_dt[, allele] %>% unique
  print(allele_name)
  featurization(allele_dt, mhc_pseudoseq, paste0("./data/MHCFlurry_features/", allele_name, "_Train_BLOSUM62_extended_v3.csv"))
})

data <- fread("./data/MHCflurry_data/train_data_v3.csv")
data[, peptide := lapply(peptide, sequences_to_fixed_length_index_encoded_array)]
data <- data[allele == "HLA-B*27:05"]
data_0 <- data[measurement_inequality == 0]
data_00 <- data_0[1:75000]
data_01 <- data_0[75001:dim(data_0)[1]]
data_1 <- data[measurement_inequality != 0]
rm(data)
gc()
featurization(data_1, mhc_pseudoseq, paste0("./data/MHCFlurry_features/", "HLA-B*27:05_1_Train_BLOSUM62_extended_v3.csv"))
featurization(data_00, mhc_pseudoseq, paste0("./data/MHCFlurry_features/", "HLA-B*27:05_00_Train_BLOSUM62_extended_v3.csv"))
featurization(data_01, mhc_pseudoseq, paste0("./data/MHCFlurry_features/", "HLA-B*27:05_01_Train_BLOSUM62_extended_v3.csv"))
gc()

# This is done through system calls, I cannot load the whole thing it seems!
system("head -n 1 ./data/MHCFlurry_features/'SLA-2*05:02_Train_BLOSUM62_extended_v3.csv' > ./data/MHCFlurry_features/all_train_features_v3.csv")
system("tail -q -n +2 ./data/MHCFlurry_features/*.csv >> ./data/MHCFlurry_features/all_train_features_v3.csv")
system("rm ./MHCFlurry_features/*_extended_v3.csv")

rm(data_0)
rm(data_00)
rm(data_01)
rm(data_1)
rm(per_allele_list)
rm(res)
gc()

Other <- fread("./data/MHCflurry_data/validation_data_v3.csv")
Other[, peptide := lapply(peptide, sequences_to_fixed_length_index_encoded_array)]

# This is done per allele, as the memory is huge
per_allele_list <- split(Other, by=c("allele"))
print("Featurization")
res <- lapply(per_allele_list, function(allele_dt) {
  allele_name <- allele_dt[, allele] %>% unique
  print(allele_name)
  featurization(allele_dt, mhc_pseudoseq, paste0("./data/MHCFlurry_features/", allele_name, "_Train_BLOSUM62_extended.csv"))
})

Other <- fread("./data/MHCflurry_data/test_data_v3.csv")
Other[, peptide := lapply(peptide, sequences_to_fixed_length_index_encoded_array)]

# This is done per allele, as the memory is huge
per_allele_list <- split(Other, by=c("allele"))
print("Featurization")
res <- lapply(per_allele_list, function(allele_dt) {
  allele_name <- allele_dt[, allele] %>% unique
  print(allele_name)
  featurization(allele_dt, mhc_pseudoseq, paste0("./data/MHCFlurry_features/", allele_name, "_Train_BLOSUM62_extended.csv"))
})
