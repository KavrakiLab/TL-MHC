library(data.table)
library(magrittr)
library(ranger)
library(caret)
library(e1071)
library(stringr)
library(viridis)
library(tuneRanger)

rf_and_plot <- function(train_data, classifier, plot_limits) {
  set.seed(111)
  if (classifier)
    x_labels <- c(paste0("peptide", seq(9)), paste0("allele", seq(34)))
  else
    x_labels <- c(paste0("peptide", seq(9)), paste0("allele", seq(34)), "BA")

  rf_model <- ranger(measurement_value ~ ., data = train_data, importance='impurity_corrected', 
                     probability = classifier, classification = classifier, num.trees=500)
  
  importances <- ranger::importance(rf_model)
  average_pooling <- sapply(split(importances, (seq_along(importances) - 1) %/% 21), mean) 
  max_pooling <- sapply(split(importances, (seq_along(importances) - 1) %/% 21), max) 
  names(average_pooling) <- x_labels
  names(max_pooling) <- x_labels
  rf_importance <- data.table(
    labels = x_labels,
    values = average_pooling
  )
  rf_importance$labels <- factor(rf_importance$labels, levels = x_labels)
  
  ggplot(rf_importance, aes(x=labels, y=values)) +
    geom_bar(stat="identity") +
    theme_bw() +
    geom_vline(xintercept = 9.5, linetype="dashed", color = "blue") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 12, vjust = 0.5), plot.title = element_text(hjust = 0.5, size = 20, face = "bold"), 
          strip.text.x = element_text(size = 16, face = "italic"), strip.text.y = element_text(size = 16, face = "italic"),
          axis.title.y = element_text(size = rel(1.3)), axis.title.x = element_text(size = rel(1.3)),
          axis.text.y = element_text(angle = 0, hjust = 1, size = 12)) +
    ylab("Mean Decrease in impurity") +
    xlab("Features") +
    scale_y_continuous(expand = c(0, 0), limits = plot_limits) 
}

prepare_dataset <- function(data_file, feature_file) {
  features <- fread(feature_file)
  data <- fread(data_file)
  data[, peptide := substr(sub("\\X.*", "", peptide), 1, 15)]
  
  setkey(data, peptide, allele)
  setkey(features, peptide, allele)
  
  all_data <- features[data, nomatch = 0] %>% unique(by = c("peptide", "allele", "measurement_value"))
  setcolorder(all_data, c("peptide", "allele", "pseudosequence", "measurement_value", colnames(all_data)[4:(length(colnames(all_data)) - 1)]))
  all_data
}

columns_of_interest <- c("measurement_value", 
                         sapply(seq(1, 189), function(x) {paste0("V", x)}), 
                         sapply(seq(946, 1659), function(x) {paste0("V", x)}))
train_data <- prepare_dataset("./cleaned_data_v2/PRIME_train_data.csv", "./features_v2/PRIME_BLOSUM62_features.csv")
train_data <- train_data[str_length(peptide) == 9]
train_data <- as.data.frame(train_data[, .SD, .SDcols = columns_of_interest])
rf_and_plot(train_data, TRUE, c(-0.2, 16))

columns_of_interest <- c("measurement_value", 
                         sapply(seq(1, 189), function(x) {paste0("V", x)}), 
                         sapply(seq(946, 1659), function(x) {paste0("V", x)}))
train_data <- prepare_dataset("./cleaned_data_v2/PRIME_balanced_train_data.csv", "./features_v2/PRIME_balanced_BLOSUM62_features.csv")
train_data <- train_data[str_length(peptide) == 9]
train_data <- as.data.frame(train_data[, .SD, .SDcols = columns_of_interest])
rf_and_plot(train_data, TRUE, c(-0.07, 2))

train_data <- fread("./cleaned_data_v2/PRIME_train_data.csv")
allele_list <- train_data[, .N, by = .(allele)][N > 100, allele] 
stats <- train_data[, .N, by = .(allele,measurement_value)]
stats[, measurement_value := factor(ifelse(measurement_value == 1, "Positive", "Negative"), levels = c("Positive", "Negative"))]
names(stats) <- c("allele", "Immunogenicity", "value")
ggplot(stats[allele %in% allele_list], aes(x=allele, y=value, fill=Immunogenicity)) +
  geom_bar(position="dodge", stat="identity") +
  theme_bw() +
  geom_vline(xintercept = 6.5, linetype="dashed", color = "blue") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 12, vjust = 0.5), plot.title = element_text(hjust = 0.5, size = 20, face = "bold"), 
        strip.text.x = element_text(size = 16, face = "italic"), strip.text.y = element_text(size = 16, face = "italic"),
        axis.title.y = element_text(size = rel(1.3)), axis.title.x = element_text(size = rel(1.3)),
        axis.text.y = element_text(angle = 0, hjust = 1, size = 12)) +
  ylab("#Data points") +
  scale_fill_viridis(discrete = TRUE) +
  xlab("Alleles") +
  scale_y_continuous(expand = c(0, 0))

train_data <- fread("./cleaned_data_v2/DeepImmuno_train_data.csv")
names(train_data) <- c("peptide", "allele", "Immunogenicity")
ggplot(train_data[allele %in% allele_list & substr(allele, 5, 5) == "A"], aes(x=Immunogenicity, y = ..scaled.., fill=allele)) +
  stat_density(alpha=0.4, colour = "black") + 
  theme_bw() +
  facet_wrap(vars(allele)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 12, vjust = 0.5), plot.title = element_text(hjust = 0.5, size = 20, face = "bold"), 
        strip.text.x = element_text(size = 16, face = "italic"), strip.text.y = element_text(size = 16, face = "italic"),
        axis.title.y = element_text(size = rel(1.3)), axis.title.x = element_text(size = rel(1.3)),
        axis.text.y = element_text(angle = 0, hjust = 1, size = 12), legend.position = "none",) +
  ylab("Density") +
  xlab("Immunogenicity Potential") 