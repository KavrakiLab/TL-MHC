# TL-MHC

The following repo consists of 2 main tools + an auxiliary one, namely:
- TLStab
- TLImm
- *TLBind

## Usage

To use these tools, you will need to have `conda` installed. When conda is sucessfully install, you can create an enviroment by simply running:
```
  conda env create -f TLenv.yml
```
in the folder, followed by
```
  conda activate TLenv
```
in order to run the tools. 

## TLStab

TLStab is our transfer learning-based pMHC stability predictor. You will find the corresponding files and instructions on how to use in the `TLStab` folder. 


## TLImm

TLImm is our transfer learning-based peptide immunogenicity predictor. You will find the corresponding files and instructions on how to use in the `TLImm` folder.

## *TLBind

TLBind is our pretrained binding affinity/eluted ligard predictor. As, during development, we focused on the downstream tasks of pMHC stability and peptide immunogenicity predictions, TLBind could be underperforming in comparison to state-of-the-art binding affinity predictors. However, it is fully usable if you desire to use it. 
