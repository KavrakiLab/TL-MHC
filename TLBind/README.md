# TLBind

TLBind is our pretrained binding affinity/eluted ligard predictor. As, during development, we focused on the downstream tasks of pMHC stability and peptide immunogenicity predictions, TLBind could be underperforming in comparison to state-of-the-art binding affinity predictors. However, it is fully usable if you desire to use it.

## Usage

To use TLBind, you will first need to be under the `TLenv` environment (see instructions on the top README file on how to install and activate the `TLenv` conda environment. 

An input .csv file that has at least an `allele` and a `peptide` column is needed (an example of SARS-CoV-2 peptides can be found in the `examples` folder). 

Given the .csv, TLBind can be used as follows:

```
  python TLBind.py examples/SARS-CoV-2_peptides_example.csv
```
This will result in an output file named `TLBind_out.csv` with binding affinity predictions (3rd column) as well as peptide presentation probabilities (4th column).

To specify a specific output file name, one can use something like the following:
```
  python TLBind.py examples/SARS-CoV-2_peptides_example.csv --out Desired_output_filename.csv
```
