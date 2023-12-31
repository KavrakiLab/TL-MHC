# TL-MHC

The following repo consists of 2 main tools + an auxiliary one, namely:
- TLStab
- TLImm
- TLBind

## Usage

To use these tools, you will need to have either `conda` or `mamba` installed. We highly recommend mamba, as it handles dependencies much faster and more reliably. To install `mamba` follow the instructions [here](https://mamba.readthedocs.io/en/latest/mamba-installation.html). If you prefer to use conda instead, you can find instructions [here](https://conda.io/projects/conda/en/latest/user-guide/install/index.html), if it's not already installed in your machine. 

First, clone the repo, and get into the main folder of the repo (where this README file is found). When either `conda` or `mamba` are sucessfully installed, you can create an enviroment by simply running:

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

#### Important note for Windows Installations

It seems that, when using `conda` on Windows, there is a conflict between the `pytorch` installation and with some of the other packages. This happened to us on Windows only. A workaround that worked for us was to remove the `pytorch` entry from the YAML environment file, and install `pytorch` later when into the environment, using a `pytorch` installation command taken from [here](https://pytorch.org/get-started/locally/). Please reach out to us if you have any issues installing the required packages. 

## TLStab

TLStab is our transfer learning-based pMHC stability predictor. You will find the corresponding files and instructions on how to use in the `TLStab` folder. 


## TLImm

TLImm is our transfer learning-based peptide immunogenicity predictor. You will find the corresponding files and instructions on how to use in the `TLImm` folder.

## TLBind

TLBind is our pretrained binding affinity/eluted ligand predictor. You will find the corresponding files and instructions on how to use in the `TLBind` folder.
