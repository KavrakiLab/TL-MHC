import pandas as pd
import numpy as np
import math

import torch

# Local imports
from dataloader import get_dataloader_eval
from model import Pretrained_BAEL

import argparse
import sys

# Basic routine for forwarding
def forwarding(model, loader, fold):

	model.eval()
	pred_list = []
	for features, peptide, allele in loader: 
		out = model(features)
		predicted_value = torch.sigmoid(out[:, 0]).item()
		pred_list.append((peptide[0], allele[0], predicted_value))

	return pd.DataFrame(pred_list, columns =['peptide', 'allele', 'Imm' + str(fold)])

def main(args):

	# Input arguments
	parser = argparse.ArgumentParser(description="TLImm arguments", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('input_file', type=str, nargs=1, help='Input_file (required!)')
	parser.add_argument('--kt_method', default="finetune", choices=["baseline", "consensus", "BA", "finetune"], help='Choice of knowledge transfer method (default is finetune, corresponding to TLImm)', type=str)
	parser.add_argument('--ablation', default="Balanced_continuous", choices=["Unbalanced_binary", "Balanced_binary", "Unbalanced_continuous", "Balanced_continuous"], help='Choice of model in regards to the dataset that the model was was trained (default is finetune, corresponding to TLImm + Bal. Cont.)', type=str)
	parser.add_argument('--l', default=0.14, help='When using Consensus scroing method, define the mixing lambda', type=float)
	parser.add_argument('--out', default="TLImm_out.csv", help='File with predictions from TLImm', type=str)
	args = parser.parse_args()

	# Initialize
	BA = False
	if args.kt_method == "BA": BA = True
	folder = args.kt_method
	if folder == "consensus": folder = "baseline"
	latent_dim = 512
	if args.kt_method == "BA": latent_dim = 128
	loader, num_features = get_dataloader_eval(args.input_file[0], BA)

	if args.ablation == "Unbalanced_binary": args.l = 0.0
	if args.ablation == "Balanced_binary": args.l = 0.33
	if args.ablation == "Unbalanced_continuous": args.l = 0.0
	if args.ablation == "Balanced_continuous": args.l = 0.14

	# Load model and predict for each fold
	prediction_list = []
	for fold in range(1, 11):
		model = Pretrained_BAEL(num_features = num_features, latent_dim=latent_dim) # Latent dimension is 512 for the best model
		model.load_state_dict(torch.load('./models/' + folder + '/' + args.ablation + '/IMM_MLP_' + str(fold) + '.pt', map_location=torch.device('cpu')))

		prediction_list.append(forwarding(model, loader, fold))

	predictions = prediction_list[0]
	for fold in range(1, 10):
		predictions = predictions.merge(prediction_list[fold], on = ['peptide', 'allele'])
	
	predictions['Immunogenicity'] = predictions.iloc[:, 2:].mean(axis=1)	
	predictions.loc[predictions['Immunogenicity'] < 0, 'Immunogenicity'] = 0.0
	
	if args.kt_method == "consensus":
		data = pd.read_csv(args.input_file[0])

		if 'peptide' not in data or 'allele' not in data:
			print('Error: Please include a file with a peptide,allele and a BA column')
			sys.exit(0)

		predictions = predictions.merge(data, on=['peptide', 'allele'])
		predictions['Immunogenicity'] = (1 - args.l)*predictions['Immunogenicity'] + args.l*predictions['BA'] 

	output_file = predictions[['peptide', 'allele', 'Immunogenicity']].to_csv(args.out, index=False)
if __name__ == "__main__":
    main(sys.argv[1:])
