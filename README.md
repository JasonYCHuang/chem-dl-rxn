# chem-dl-rxn


This repository containing [1] training deep learning models [2] a flask inference server.

We use the similar method in this [paper](https://pubs.acs.org/doi/10.1021/acscentsci.9b00576).


### 1. Training deep learning models

we use pytorch and fairseq for the transformer model.

```bash
$ cd training_forward_retro
```

Then, please follow steps in:
1. `training_forward_retro/INSTALL_TRAINING.md` for installation
2. `training_forward_retro/README.md` for training reaction prediction
3. `training_forward_retro/README-RETRO.md` for training retrosynthesis


### 2. Run a flask inference server

Please follow steps in:
1. `INSTALL_SERVER.md` for installation


### 3. Some detail

Dataset of reaction prediction
- USPTO_STEREO dataset, mixed starting materials and reactants SMILES as the model input; products SMILES as the model output.
- Train/Valid/Test = 902K / 50K / 50K

Dataset of retrosynthesis
- USPTO_STEREO dataset, products SMILES as the model input; starting materials SMILES as the model output.
- Train/Valid/Test = 902K / 50K / 50K

Accuracy of reaction prediction = 72.2%

Accuracy of retrosynthesis = 41.2% 
