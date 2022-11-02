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
		predicted_value = out[:, 0].item()
		if predicted_value > 1:
			pred_list.append((peptide[0], allele[0], (-1)/math.log2(0.99)))
		elif predicted_value > 0:
			pred_list.append((peptide[0], allele[0], (-1)/math.log2(predicted_value)))
		else:
			pred_list.append((peptide[0], allele[0], 0))

	return pd.DataFrame(pred_list, columns =['peptide', 'allele', 'Stab' + str(fold)])

def main(args):

	# Input arguments
	parser = argparse.ArgumentParser(description="TLStab arguments", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('input_file', type=str, nargs=1, help='Input_file (required!)')
	parser.add_argument('--kt_method', default="finetune", choices=["baseline", "consensus", "BA", "finetune"], help='Choice of knowledge transfer method (default is finetune, corresponding to TLStab)', type=str)
	parser.add_argument('--l', default=0.26, help='When using Consensus scroing method, define the mixing lambda', type=float)
	parser.add_argument('--out', default="TLStab_out.csv", help='File with predictions from TLStab', type=str)
	args = parser.parse_args()

	# Initialize
	BA = False
	if args.kt_method == "BA": BA = True
	folder = args.kt_method
	if folder == "consensus": folder = "baseline"
	latent_dim = 512
	if args.kt_method == "BA": latent_dim = 128
	loader, num_features = get_dataloader_eval(args.input_file[0], BA)

	# Load model and predict for each fold
	prediction_list = []
	for fold in range(1, 11):
		model = Pretrained_BAEL(num_features = num_features, latent_dim=latent_dim) # Latent dimension is 512 for the best model
		model.load_state_dict(torch.load('./models/' + folder + '/STAB_MLP_' + str(fold) + '.pt', map_location=torch.device('cpu')))

		prediction_list.append(forwarding(model, loader, fold))

	predictions = prediction_list[0]
	for fold in range(1, 10):
		predictions = predictions.merge(prediction_list[fold], on = ['peptide', 'allele'])
	
	predictions['Half-life (h)'] = predictions.iloc[:, 2:].mean(axis=1)	
	predictions.loc[predictions['Half-life (h)'] < 0, 'Half-life (h)'] = 0.0
	
	if args.kt_method == "consensus":
		data = pd.read_csv(args.input_file[0])

		if 'peptide' not in data or 'allele' not in data:
			print('Error: Please include a file with a peptide,allele and a BA column')
			sys.exit(0)

		predictions = predictions.merge(data, on=['peptide', 'allele'])
		predictions['Half-life (h)'] = (1 - args.l)*predictions['Half-life (h)'] + args.l*predictions['BA'] 

	output_file = predictions[['peptide', 'allele', 'Half-life (h)']].to_csv(args.out, index=False)
if __name__ == "__main__":
    main(sys.argv[1:])
