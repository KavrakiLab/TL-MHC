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

The default usage of TLImm is with the transfer learning protocol that we discuss in our paper. However, we also provide all the knowledge transfer models discussed in the paper:
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
  python TLImm.py examples/SARS-CoV-2_peptides_BA_example.csv --kt_method consensus --l 0.5
```

In addition to that, we also provide the different instances of TLImm (+ other knowledge transfer methods) when trained on different datasets (Balanced-unbalanced or binary labeled-continuously labeled). This analysis can be seen in the ablation section in our paper. You can use those instances by changing the `--ablation` flag:
- A) **Unbalanced\_binary**
- B) **Balanced\_binary**
- C) **Unbalanced\_continuous**
- D) **Balanced\_continuous**

```
  python TLImm.py examples/SARS-CoV-2_peptides_example.csv --ablation Unbalanced_binary
```

Finally, to specify a specific output file name, one can use something like the following:
```
  python TLImm.py examples/SARS-CoV-2_peptides_example.csv --out Desired_output_filename.csv
```
