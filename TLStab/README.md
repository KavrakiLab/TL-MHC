# TLStab

TLStab is our transfer learning-based pMHC stability predictor.

## Usage

To use TLStab, you will first need to be under the `TLenv` environment (see instructions on the top README file on how to install and activate the `TLenv` conda environment. 

An input .csv file that has at least an `allele` and a `peptide` column is needed (an example of SARS-CoV-2 peptides can be found in the `examples` folder). 

Given the .csv, TLStab can be used as follows:

```
  python TLStab.py examples/SARS-CoV-2_peptides_example.csv
```
This will result in an output file named `TLStab_out.csv` with pMHC stability predictors for each peptide-allele pair (in half life minutes) in the 3rd column.

Finally, to specify a specific output file name, one can use something like the following:
```
  python TLStab.py examples/SARS-CoV-2_peptides_example.csv --out Desired_output_filename.csv
```