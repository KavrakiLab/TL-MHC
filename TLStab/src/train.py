import sys
import gc
import numpy as np
import random

import torch
import torch.nn.functional as F
from torch import optim

from immuno_dataloader import get_immuno_dataset, get_features, get_immuno_dataloader, EarlyStopping
from immuno_models import MLP_multitasker

from sklearn.metrics import mean_squared_error, average_precision_score, log_loss
from math import sqrt
from scipy import stats

def forwarding(model, loader, optimizer, size, action='train', effort='', frozen = ''):

    # Continue with the forwarding operation:
    forwarding_loss = 0
    preds_list = []
    labels_list = []
    if action == 'train':
        model.train()
        torch.set_grad_enabled(True)
        if frozen in ["freeze", "new_layer"]:
            model.fc1.weight.requires_grad = False
            model.fc1.bias.requires_grad = False
            model.fc2.weight.requires_grad = False
            model.fc2.bias.requires_grad = False
    else:
        model.eval()
        torch.set_grad_enabled(False)
    for features, measurement_value in loader:
        if action=='train': optimizer.zero_grad()
        out = model(features.cuda())
        loss = F.mse_loss(out[0][0].cpu(), measurement_value[0], reduction = 'mean')
        adding_loss = F.mse_loss(out[0][0].cpu(), measurement_value[0], reduction = 'sum')
        if action == 'train':
            loss.backward()
            optimizer.step()
        else:
            preds_list.append(out[0][0].cpu())
            labels_list.append(measurement_value[0])
        forwarding_loss += adding_loss.item()
    forwarding_loss = forwarding_loss / size
    return forwarding_loss, np.array(labels_list), np.array(preds_list)

def main(args):

    # Initialization
    np.random.seed(0)
    random.seed(0)
    torch.manual_seed(0)
    torch.cuda.manual_seed(0)
    torch.cuda.manual_seed_all(0)
    torch.backends.cudnn.enabled = False
    torch.backends.cudnn.benchmark = False
    torch.backends.cudnn.deterministic = True

    # Configuration (configuration file in the future!!)
    print(args[0])
    print(args[1])
    print(args[2])
    print(args[3])
    print(args[4])
    print(args[5])
    print(args[6])
    batch_size = int(args[0])
    lr = float(args[1])
    latent_dim = int(args[2])
    weight_decay = float(args[3])
    model_type = args[4]
    fold = int(args[5])
    weight_init = args[6]
    epochs = 10000
    device = torch.device('cuda')

    # Make selection based on the dataset

    if model_type == "BA":
        feature_type = "BA"
    else:
        feature_type = ""

    # Identifier of process
    effort = 'Stab_' + str(batch_size) + '_' + str(lr) + '_' + str(latent_dim) + '_' + str(weight_decay) + '_' + str(model_type) + '_fold' + str(fold)

    idx_dict_train = get_immuno_dataset("./Nested_CVs_stab/folds/train_" + str(fold) + ".csv")
    idx_dict_val = get_immuno_dataset("./Nested_CVs_stab/folds/test_" + str(fold) + ".csv")

    feats_dict, num_features = get_features("./Nested_CVs_stab/features/Stab_BLOSUM62_features" + feature_type + ".csv")

    train_peptide_list = [x[0] for x in list(idx_dict_train.values())]
    valid_peptide_list = [x[0] for x in list(idx_dict_val.values())]

    train_size = len(train_peptide_list)
    valid_size = len(valid_peptide_list)

    train_loader = get_immuno_dataloader(idx_dict_train,
                                         {key: feats_dict[key] for key in train_peptide_list},
                                         batch_size=batch_size, shuffle=True)

    valid_loader = get_immuno_dataloader(idx_dict_val,
                                         {key: feats_dict[key] for key in valid_peptide_list},
                                         batch_size=1, shuffle=False)

    del train_peptide_list
    del valid_peptide_list
    del feats_dict
    del idx_dict_train
    del idx_dict_val
    gc.collect()
    gc.collect()

    model = MLP_multitasker(num_features=num_features, latent_dim=latent_dim, output_dim = 2)
    model_path = './Stab_models_paper_revision_v3/' + str(fold) + '/torch_models/' + weight_init + '_fold' + str(fold) + '.pt'
    model.load_state_dict(torch.load(model_path))
    if model_type in ["freeze", "new_layer"]:
        model.fc1.weight.requires_grad = False
        model.fc1.bias.requires_grad = False
        model.fc2.weight.requires_grad = False
        model.fc2.bias.requires_grad = False
    model.cuda()
    optimizer = optim.SGD(filter(lambda p: p.requires_grad, model.parameters()), lr=lr, weight_decay=weight_decay)
    early_stopping = EarlyStopping()
    best_val_loss = None
    
    for epoch in range(epochs):

        # Training section
        train_loss, labels_list, preds_list = forwarding(model, train_loader, optimizer, train_size, action='train', effort=effort, frozen = model_type)
        print('Epoch: {}, Train Average loss: {:.4f}'.format(epoch, train_loss))
        with open('/home/rf27/Stab_models_whole_paper_revision_v3/' + str(fold) + '/logs/' + str(effort) + '.txt', 'a') as f:
            f.write('Train epoch: {}, Train Average loss: {:.4f}\n'.format(epoch, train_loss))

        # Validation section
        val_loss, labels_list, preds_list = forwarding(model, valid_loader, optimizer, valid_size, action='valid', effort=effort, frozen = model_type)
        print('Validation Average loss: {:.4f}'.format(val_loss))
        with open('/home/rf27/Stab_models_whole_paper_revision_v3/' + str(fold) + '/logs/' + str(effort) + '.txt', 'a') as f:
            f.write('Validation Average loss: {:.4f}\n'.format(val_loss))

        print("MSE: ", mean_squared_error(labels_list, preds_list))
        print("Spearman CC: ", stats.spearmanr(labels_list, preds_list)[0])

        with open('/home/rf27/Stab_models_whole_paper_revision_v3/' + str(fold) + '/logs/' + str(effort) + '.txt', 'a') as f:
            f.write("MSE: " + str(mean_squared_error(labels_list, preds_list)) + ", Spearman CC: " + str(stats.spearmanr(labels_list, preds_list)[0]) + "\n")

        # Save the best-so-far model
        if best_val_loss == None or val_loss < best_val_loss:
            best_val_loss = val_loss
            torch.save(model.state_dict(), '/home/rf27/Stab_models_whole_paper_revision_v3/' + str(fold) + '/torch_models/' + str(effort) + '.pt')

        # Early Stopping
        early_stopping(val_loss)
        if early_stopping.early_stop:
            break

if __name__ == "__main__":
    main(sys.argv[1:])
