library(data.table)
library(Peptides)
library(stringr)
library(magrittr)
library(Biostrings)
library(plyr)

# Load the data
data <- fread("./data/MHCflurry_data/MULTIALLELIC_OLD_dataset.csv")
data <- data[, Indexx := seq(1, dim(data)[1])]
data <- data[, .(V1, V2, V7, V10, V12, V14, V16, Indexx)]
names(data) <- c("Cell_Line", "Peptide", "Label", "Predictor1", "Predictor2", "Predictor3", "Predictor4", "Indexx")
data <- melt.data.table(data, id.vars = c("Cell_Line", "Peptide", "Label", "Indexx"))
data <- data[, .(Cell_Line, Peptide, Label, Indexx, value)]
data[, value := as.factor(value)]
levelss <- data[, value] %>% unique()
data <- ddply(data, .(Cell_Line, Peptide, Label, Indexx), function(x) which.max(tabulate(x$value)))
data <- as.data.table(data)
data[, V1 := levels(levelss)[V1]]
data <- data[, .(Peptide, V1, Label)]
names(data) <- c("Peptide", "MHC", "hit")
fwrite(data, "./data/MHCflurry_data/Deconvoluted_old.csv")
