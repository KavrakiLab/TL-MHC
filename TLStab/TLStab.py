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
	parser.add_argument('--out', default="TLStab_out.csv", help='File with predictions from TLStab', type=str)
	args = parser.parse_args()

	# Initialize
	loader, num_features = get_dataloader_eval(args.input_file[0], BA)

	# Load model and predict for each fold
	prediction_list = []
	for fold in range(1, 11):
		model = Pretrained_BAEL(num_features = num_features, latent_dim=512) # Latent dimension is 512 for the best model
		model.load_state_dict(torch.load('./models/weights/TLStab_' + str(fold) + '.pt', map_location=torch.device('cpu')))

		prediction_list.append(forwarding(model, loader, fold))

	predictions = prediction_list[0]
	for fold in range(1, 10):
		predictions = predictions.merge(prediction_list[fold], on = ['peptide', 'allele'])
	
	predictions['Half-life (h)'] = predictions.iloc[:, 2:].mean(axis=1)	
	predictions.loc[predictions['Half-life (h)'] < 0, 'Half-life (h)'] = 0.0

	output_file = predictions[['peptide', 'allele', 'Half-life (h)']].to_csv(args.out, index=False)
if __name__ == "__main__":
    main(sys.argv[1:])
