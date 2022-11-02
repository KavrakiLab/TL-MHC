import pandas as pd
import numpy as np
from math import sqrt
from scipy import stats

import torch
from torch import optim, nn
import torch.nn.functional as F
from torch.utils.data import DataLoader
from torch.utils.data import Dataset

# Local imports
from dataloader import get_dataset, get_features, get_dataloader_train, EarlyStopping
from model import Pretrained_BAEL

import gc
import sys

# Basic routine for forwarding
def forwarding(model, loader, optimizer, num_of_BA, num_of_MS, num_of_MR, action='train', effort=''):
    loss = 0
    bce_loss = 0
    mse_loss = 0
    mr_loss = 0
    if action == 'train':
        model.train()
    else:
        model.eval()
    for features, measurement_kind, measurement_inequality, measurement_value in loader:
        if action=='train': optimizer.zero_grad() 
        out = model(features.cuda())
        bce_el_loss = F.binary_cross_entropy(torch.sigmoid(out[measurement_kind, 1].cpu()), 
                                             measurement_value[measurement_kind],
                                             reduction='sum')
        bce_el_loss = torch.nan_to_num(bce_el_loss)
        mean_bce_el_loss = bce_el_loss / measurement_value[measurement_kind].shape[0]
        mse_ba_loss = F.mse_loss(out[~measurement_kind & (measurement_inequality == 0), 0].cpu(),
                                measurement_value[~measurement_kind & (measurement_inequality == 0)],
                                reduction='sum')
        mse_ba_loss = torch.nan_to_num(mse_ba_loss)
        mean_mse_ba_loss = mse_ba_loss / measurement_value[~measurement_kind & (measurement_inequality == 0)].shape[0]
        mr = F.margin_ranking_loss(out[~measurement_kind & (measurement_inequality != 0), 0].cpu(),
                        measurement_value[~measurement_kind & (measurement_inequality != 0)],
                        (-measurement_inequality[~measurement_kind & (measurement_inequality != 0)]),
                        reduction='sum')
        mr = torch.nan_to_num(mr)
        mean_mr = mr / (-measurement_inequality[~measurement_kind & (measurement_inequality != 0)]).shape[0]
        loss = mean_bce_el_loss + mean_mse_ba_loss + mean_mr
        if action == 'train':
            loss.backward()
            optimizer.step()
        bce_loss += bce_el_loss.item()
        mse_loss += mse_ba_loss.item()
        mr_loss += mr.item()
    bce_loss = bce_loss / num_of_MS
    mse_loss = mse_loss / num_of_BA
    try: 
        mr_loss = mr_loss / num_of_MR
    except ZeroDivisionError:
        mr_loss = 0
    regression_loss = bce_loss + mse_loss + mr_loss
    return regression_loss, bce_loss, mse_loss, mr_loss

def main(args):
    
    batch_size = int(args[0])
    lr = float(args[1])
    latent_dim = int(args[2])
    weight_decay = float(args[3])
    epochs = 3000
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    #device = 'cpu'
    effort = 'BA_EL_' + str(batch_size) + '_' + str(lr) + '_' + str(latent_dim) + '_' + str(weight_decay) + '_'
    
    # Load data and make the loaders
    idx_dict_train, train_BA, train_MS, train_MR = get_dataset("../cleaned_data/cleaned_train_data.csv")
    idx_dict_valid, valid_BA, valid_MS, valid_MR = get_dataset("../cleaned_data/test_data.csv")

    train_peptide_list = [x[0] for x in list(idx_dict_train.values())]
    valid_peptide_list = [x[0] for x in list(idx_dict_valid.values())]
    peptide_list = train_peptide_list + valid_peptide_list

    feats_dict, num_features = get_features(train_features_path="../features/all_train_features.csv",
                                            peptide_list=train_peptide_list)
    
    test_feats_dict, num_features = get_features(train_features_path="../features/Test_BLOSUM62_extended.csv", peptide_list=valid_peptide_list)

    train_loader = get_dataloader_train(idx_dict_train, 
                                        {key: feats_dict[key] for key in train_peptide_list},
                                        batch_size=batch_size, shuffle=True)

    valid_loader = get_dataloader_train(idx_dict_valid, 
                                        {key: test_feats_dict[key] for key in valid_peptide_list},
                                        batch_size=batch_size, shuffle=False)
    
    del train_peptide_list
    del valid_peptide_list
    del peptide_list
    del feats_dict
    del test_feats_dict
    del idx_dict_train
    del idx_dict_valid
    gc.collect()
    gc.collect()
    
    #Define Model 
    model = Pretrained_BAEL(num_features=num_features, latent_dim=latent_dim)
    model.cuda()
    optimizer = optim.SGD(model.parameters(), lr=lr, weight_decay=weight_decay)
    early_stopping = EarlyStopping()
    best_val_loss = None

    for epoch in range(epochs):
    
        # Training section
        regression_loss, bce_loss, mse_loss, mr_loss = forwarding(model, train_loader, optimizer,                                                                       train_BA, train_MS, train_MR,
                                                                  action='train', effort=effort)
        print('Epoch: {}, Train Average loss: {:.4f}'.format(epoch, regression_loss))
        print('BCE Loss: {:.4f}, MSE loss: {:.4f}, Margin Loss {:.4f}'.format(bce_loss, 
                                                                              mse_loss,
                                                                              mr_loss))
        
        with open('../models/logs/' + str(effort) + '.txt', 'a') as f:
            f.write('Train epoch: {}, Train Average loss: {:.4f}\n'.format(epoch, regression_loss))
            f.write('BCE Loss: {:.4f}, MSE loss: {:.4f}, Margin Loss {:.4f}\n'.format(bce_loss, 
                                                                                      mse_loss,
                                                                                      mr_loss))
            
        # Validation section
        val_regression_loss, bce_loss, mse_loss, mr_loss = forwarding(model, valid_loader, optimizer, 
                                                                      valid_BA, valid_MS, valid_MR,
                                                                      action='valid', effort=effort)
        print('Validation Average Regression loss: {:.4f}'.format(val_regression_loss))
        print('BCE Loss: {:.4f}, MSE loss: {:.4f}, Margin Loss {:.4f}'.format(bce_loss,
                                                                              mse_loss,
                                                                              mr_loss))
        with open('../models/logs/' + str(effort) + '.txt', 'a') as f:
            f.write('Validation Average Regression loss: {:.4f}\n'.format(val_regression_loss))
            f.write('BCE Loss: {:.4f}, MSE loss: {:.4f}, Margin Loss {:.4f}\n'.format(bce_loss, 
                                                                                      mse_loss,
                                                                                      mr_loss))
        # Save the best-so-far model
        if best_val_loss == None or val_regression_loss < best_val_loss:
            best_val_loss = val_regression_loss
            torch.save(model.state_dict(), '../models/torch_models/' + str(effort) + '.pt')
        
        # Early Stopping
        early_stopping(val_regression_loss)
        if early_stopping.early_stop:
            break

if __name__ == "__main__":
    main(sys.argv[1:])