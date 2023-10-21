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
		pred_list.append((peptide[0], allele[0], predicted_value))

	return pd.DataFrame(pred_list, columns =['peptide', 'allele', 'Imm' + str(fold)])

def main(args):
	
	print("Initializing TLImm!")
	# Input arguments
	parser = argparse.ArgumentParser(description="TLImm arguments", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('input_file', type=str, nargs=1, help='Input_file (required!)')
	parser.add_argument('--out', default="TLImm_out.csv", help='File with predictions from TLImm', type=str)
	args = parser.parse_args()

	# Initialize
	print("Reading input file")
	loader, num_features = get_dataloader_eval(args.input_file[0], False)

	# Load model and predict for each fold
	print("Generating predictions...")
	prediction_list = []
	for fold in range(1, 16):
		model = Pretrained_BAEL(num_features = num_features, latent_dim=512) # Latent dimension is 512 for the best model
		model.load_state_dict(torch.load('./models/weights/TLImm_' + str(fold) + '.pt', map_location=torch.device('cpu')))

		prediction_list.append(forwarding(model, loader, fold))

	predictions = prediction_list[0]
	for fold in range(1, 15):
		predictions = predictions.merge(prediction_list[fold], on = ['peptide', 'allele'])
	
	predictions['Immunogenicity'] = predictions.iloc[:, 2:].mean(axis=1)	
	predictions.loc[predictions['Immunogenicity'] < 0, 'Immunogenicity'] = 0.0
	output_file = predictions[['peptide', 'allele', 'Immunogenicity']].to_csv(args.out, index=False)
	print("Output generated in " + args.out)
	print("End of TLImm!")

if __name__ == "__main__":
    main(sys.argv[1:])
