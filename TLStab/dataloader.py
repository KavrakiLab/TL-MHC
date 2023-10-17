import numpy as np
import torch
import pandas as pd
from torch.utils.data import DataLoader
from torch.utils.data import Dataset
from Bio.Align import substitution_matrices
from sklearn.preprocessing import MinMaxScaler

import sys

class EarlyStopping():
    """
    Early stopping to stop the training when the loss does not improve after
    certain epochs.
    """
    def __init__(self, patience=500, min_delta=0):
        """
        :param patience: how many epochs to wait before stopping when loss is
               not improving
        :param min_delta: minimum difference between new loss and old loss for
               new loss to be considered as an improvement
        """
        self.patience = patience
        self.min_delta = min_delta
        self.counter = 0
        self.best_loss = None
        self.early_stop = False
        
    def __call__(self, val_loss):
        if self.best_loss == None:
            self.best_loss = val_loss
        elif self.best_loss - val_loss > self.min_delta:
            self.best_loss = val_loss
            # reset counter if validation loss improves
            self.counter = 0
        elif self.best_loss - val_loss < self.min_delta:
            self.counter += 1
            print(f"INFO: Early stopping counter {self.counter} of {self.patience}")
            if self.counter >= self.patience:
                print('INFO: Early stopping')
                self.early_stop = True

class FeaturePairDataset_train(Dataset):
    # Class that defines how DataLoader should be working during training
    
    def __init__(self, idx_dict, feats_dict):
        self.idx_dict = idx_dict 
        self.feats_dict = feats_dict
        
    def __len__(self):
        return len(self.idx_dict)

    def __getitem__(self, idx):
        (hlapeptide, measurement_kind, measurement_inequality, measurement_value) = self.idx_dict[idx]
        data = torch.from_numpy(np.asarray(self.feats_dict[hlapeptide]).astype('float32'))
        measurement_kind = torch.from_numpy(np.asarray(measurement_kind).astype('bool_'))
        measurement_inequality = torch.from_numpy(np.asarray(measurement_inequality).astype('int32'))
        measurement_value = torch.from_numpy(np.asarray(measurement_value).astype('float32'))
        return data, measurement_kind, measurement_inequality, measurement_value

class FeaturePairDataset_eval(Dataset):
    # Class that defines how DataLoader should be working during evaluation
    
    def __init__(self, idx_dict, feats_dict):
        self.idx_dict = idx_dict 
        self.feats_dict = feats_dict
        
    def __len__(self):
        return len(self.idx_dict)

    def __getitem__(self, idx):
        (hlapeptide, peptide, allele) = self.idx_dict[idx]
        data = torch.from_numpy(np.asarray(self.feats_dict[hlapeptide]).astype('float32'))
        return data, peptide, allele
    
def get_dataset(path):
    
    # Returns dictionary of dataset given path
    
    DT = pd.read_csv(path)
    num_of_MS = DT[DT['measurement_kind'] == 1].shape[0]
    num_of_BA = DT[(DT['measurement_kind'] == 0) & (DT['measurement_inequality'] == 0)].shape[0]
    num_of_MR = DT[(DT['measurement_kind'] == 0) & (DT['measurement_inequality'] != 0)].shape[0]
    
    DT['idx'] = range(0, DT.shape[0])
    DT = DT.set_index('idx')
    DT['key'] = DT['peptide'] + '-' + DT['allele']
    idx_dict = DT[['key', 'measurement_kind', 
                   'measurement_inequality', 'measurement_value']].to_dict('index')
    for key, value in idx_dict.items():
        idx_dict[key] = list(value.values())
    
    return idx_dict, num_of_BA, num_of_MS, num_of_MR

def get_features(train_features_path, peptide_list):
    
    # Returns the features of the dataset
    
    dtype_dict = {}
    for i in range(1, 1659):
        dtype_dict["V" + str(i)] = np.float16
    DT_feat = pd.read_csv(train_features_path, dtype = dtype_dict).drop_duplicates(subset=['peptide', 'allele'])
    
    #Filtering out stuff that is not going to be used:
    DT_feat['key'] = DT_feat['peptide'] + '-' + DT_feat['allele']
    DT_feat = DT_feat[DT_feat['key'].isin(peptide_list)]
    
    peptide_allele_list = list(DT_feat['key'])
    feats = np.asarray(DT_feat.drop(columns=['key', 'peptide', 'allele', 'pseudosequence']))
    num_features = feats.shape[1]
    feats_dict = {}
    for i, peptide in enumerate(peptide_allele_list):
        feats_dict[peptide] = feats[i, :]
    
    return feats_dict, num_features

def get_dataloader_train(idx_dict, feats_dict, batch_size, shuffle):
    
    # Returns dataloader given dataset
    
    data = FeaturePairDataset_train(idx_dict=idx_dict, feats_dict=feats_dict)
    loader = DataLoader(data, batch_size=batch_size, shuffle=shuffle, collate_fn=None, num_workers=0)
    
    return loader

def get_dataloader_eval(path, BA):

    # Given a .csv file, return a data point dictionary and a feature dictionary
    DT = pd.read_csv(path)
    if 'peptide' not in DT or 'allele' not in DT:
        print('Error: Please include a file with a peptide and an allele column')
        sys.exit(0)

    # First, the data instance dictionary
    DT = pd.read_csv(path)
    DT['idx'] = range(0, DT.shape[0])
    DT = DT.set_index('idx')
    DT['key'] = DT['peptide'] + '-' + DT['allele']
    idx_dict = DT[['key', 'peptide', 'allele']].to_dict('index')
    for key, value in idx_dict.items():
        idx_dict[key] = list(value.values())

    # Control:
    if BA:
        if 'BA' not in DT:
            print('Error: Please include a file with a BA column')
            sys.exit(0)

    # Second, generate the feature dictionary
    blosum62 = np.array(substitution_matrices.load("BLOSUM62"))[:20, :20]
    aa_codes = {
        'A' : 0,
        'R' : 1,
        'N' : 2,
        'D' : 3,
        'C' : 4,
        'Q' : 5,
        'E' : 6,
        'G' : 7,
        'H' : 8,
        'I' : 9,
        'L' : 10,
        'K' : 11,
        'M' : 12,
        'F' : 13,
        'P' : 14,
        'S' : 15,
        'T' : 16,
        'W' : 17,
        'Y' : 18,
        'V' : 19,
        'X' : 20
    }
    scaler = MinMaxScaler().fit(blosum62)
    scaled_blosum62 = scaler.transform(blosum62)
    scaled_blosum62 = np.vstack([scaled_blosum62, np.array([0]*20)])
    scaled_blosum62 = np.hstack([scaled_blosum62, np.array([0]*20 + [1]).reshape(-1, 1)])
    pseudoDT = pd.read_csv("./misc/updated_pseudosequences.csv")

    peptide_allele_list = DT['key'].tolist()
    if BA: BA_list = DT['BA'].tolist()
    feats_dict = {}
    for i, key in enumerate(peptide_allele_list):
        
        peptide = key.split('-')[0]
        X_length = 15 - len(peptide)
        if X_length > 0:
            aligned_peptide = peptide + ''.join(['X']*(X_length + int(np.floor(X_length / 2)))) + peptide + ''.join(['X']*(X_length + int(np.ceil(X_length / 2)))) + peptide
        else:
            aligned_peptide = ''.join([peptide]*3)
        allele = key.split('-', 1)[1]
        feature_list = []
        for amino_acid in list(aligned_peptide):
            feature_list.append(scaled_blosum62[aa_codes[amino_acid]])
        allele_pseudosequence = pseudoDT[pseudoDT['allele'] == allele]['pseudosequence'].values[0]
        for amino_acid in list(allele_pseudosequence):
            feature_list.append(scaled_blosum62[aa_codes[amino_acid]])
        if BA: feature_list.append(BA_list[i]) 
        feats_dict[key] = np.hstack(feature_list)

    # Finally, return dataloader given dataset

    data = FeaturePairDataset_eval(idx_dict=idx_dict, feats_dict=feats_dict)
    loader = DataLoader(data, batch_size=1, shuffle=False, collate_fn=None, num_workers=0)
    num_features = feats_dict[key].shape[0]

    return loader, num_features
