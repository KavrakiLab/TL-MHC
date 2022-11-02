Some explanation of the files (very broadly, more details will follow):

Dataset_processing.R : Receives the mhcflurry BA/EL data + the MONOALLELIC cell line data from HLAThensa +Abelin and makes the training/validation folds + the left-out test set. It also featurizes each of those datasets using BLOSUM62.