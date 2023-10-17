# TL-MHC

The following repo consists of 2 main tools + an auxiliary one, namely:
- TLStab
- TLImm
- TLBind

## Usage

To use these tools, you will need to have either `conda` or `mamba` installed. We highly recommend mamba, as it handles dependencies much faster and more reliably. To install `mamba` follow the instructions [here](https://mamba.readthedocs.io/en/latest/mamba-installation.html). When either `conda` or `mamba` are sucessfully installed, you can create an enviroment by simply running:

```
  mamba env create -f TLenv.yml
```
or
```
  conda env create -f TLenv.yml
```

in the folder, followed by
```
  mamba activate TLenv
```
or
```
  conda activate TLenv
```
in order to run the tools. 

## TLStab

TLStab is our transfer learning-based pMHC stability predictor. You will find the corresponding files and instructions on how to use in the `TLStab` folder. 


## TLImm

TLImm is our transfer learning-based peptide immunogenicity predictor. You will find the corresponding files and instructions on how to use in the `TLImm` folder.

## TLBind

TLBind is our pretrained binding affinity/eluted ligard predictor. You will find the corresponding files and instructions on how to use in the `TLBind` folder.
