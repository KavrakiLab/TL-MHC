# TLImm

TLImm is our transfer learning-based peptide immunogenicity predictor.

## Usage

To use TLImm, you will first need to be under the `TLenv` environment (see instructions on the top README file on how to install and activate the `TLenv` conda environment. 

An input .csv file that has at least an `allele` and a `peptide` column is needed (an example of SARS-CoV-2 peptides can be found in the `examples` folder). 

Given the .csv, TLImm can be used as follows:

```
  python TLImm.py examples/SARS-CoV-2_peptides_example.csv
```
This will result in an output file named `TLImm_out.csv` with peptide immunogenicity probabilities for each peptide-allele pair in the 3rd column.

Finally, to specify a specific output file name, one can use something like the following:
```
  python TLImm.py examples/SARS-CoV-2_peptides_example.csv --out Desired_output_filename.csv
```
