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

The default usage of TLStab is with the transfer learning protocol that we discuss in our paper. However, we also provide all the knowledge transfer models discussed in the paper:
- A) Baseline (`baseline`)
- B) Consensus scoring (`consensus`)
- C) Binding Affinity as a feature (`BA`)
- D) Transfer learning (`finetune`) (default usage)

To use these models, one needs to change the `--kt_method` flag like so:
```
  python TLStab.py examples/SARS-CoV-2_peptides_example.csv --kt_method baseline
```
In the case of either Consensus Scoring or Binding Affinity as a feature, an additional `BA` column has to be provided in the .csv file (see `examples` folder). As the models were trained on NetMHCpan4.1BA predictions, one should be using those for the best performance, however, any Binding Affinity prediction will be sufficient for usage. 

In the case of Consensus Scoring, an additional lambda parameter is needed to weight Binding Affinity and pMHC Stability predictions. The default is equal to 0.26, as this is the best one found in the cross validation process described in the paper), but other ones can be used also by setting the parameter accordingly:
```
  python TLStab.py examples/SARS-CoV-2_peptides_BA_example.csv --kt_method consensus --l 0.5
```

Finally, to specify a specific output file name, one can use something like the following:
```
  python TLStab.py examples/SARS-CoV-2_peptides_example.csv --out Desired_output_filename.csv
```
