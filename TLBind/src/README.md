### Datasets

The datasets that were used to pre-train TLBind were the following:

- All datasets that were used for training MHCFlurry2.0PS were used for the training process.
- MONOALLELIC datasets from the publication by Pyke et al. were used for validation and early stopping.
- MULTIALLELIC datasets from the publication by Pyke et al. were used for testing purposes. 

All of the datasets were used as is, and can be retireved from the MHCFlurry2.0 and Pyke et al. publications. Then, the 2 scripts in the folder can be used to process these datasets. If there is any issue with using/understanding the code provided, please reach out to the authors of the TLBind + TLStab + TLImm. 

### Some explanation of the files (very broadly, more details will follow suit):

- Dataset_processing.R : Receives the mhcflurry BA/EL data + the MONOALLELIC cell line data from HLAThensa +Abelin and makes the training/validation folds + the left-out test set. It also featurizes each of those datasets using BLOSUM62.

- Deconvolute_multiallelic_data.R : This file was used to take the MHCFlurry MULTIALLELIC-OLD and MULTIALLELIC-RECENT datasets, and deconvolutes these so that we can easily use them for training. It does this by taking the ensemble of predictions from 4 different state-of-the-art binding prediction models, and performs majority voting to assign an allele to a peptide from the ones in the cell line. 