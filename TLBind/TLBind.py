import pandas as pd
import numpy as np

import torch

# Local imports
from dataloader import get_dataloader_eval
from model import Pretrained_BAEL

import argparse
import sys

# Basic routine for forwarding
def forwarding(model, loader):

	model.eval()
	pred_list = []
	for features, peptide, allele in loader: 
		out = model(features)
		pred_list.append((peptide[0], allele[0], 50000**(1 - out[:, 0].item()), torch.sigmoid(out[:, 1]).item()))

	return pd.DataFrame(pred_list, columns =['peptide', 'allele', 'BA_score', 'EL_score'])

def main(args):

	# Input arguments
	parser = argparse.ArgumentParser(description="TLBind arguments", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('input_file', type=str, nargs=1, help='Input_file (required!)')
	parser.add_argument('--out', default="TLBind_out.csv", help='File with predictions from TLBind', type=str)
	args = parser.parse_args()

	# Initialize loader
	loader, num_features = get_dataloader_eval(args.input_file[0])

	# Load model
	model = Pretrained_BAEL(num_features = num_features, latent_dim=512) # Latent dimension is 512 for the best model
	model.load_state_dict(torch.load('./models/Pretrained_BA-EL.pt', map_location=torch.device('cpu')))

	predictions = forwarding(model, loader)

	output_file = predictions.to_csv(args.out, index=False)

if __name__ == "__main__":
    main(sys.argv[1:])