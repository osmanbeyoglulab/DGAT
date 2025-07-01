## Introduction

This repository is for DGAT Demo, including Model training, and Prediction. 

DGAT is a deep learning model for predicting protein expression from Spatial Transcriptomics (ST) data. It uses a graph attention network to learn the relationships between genes and proteins, and then predicts protein expression from ST data.

## Data 
**The training datasets for DGAT model, pre-trained model and Spatial Transcriptomics for predictions can be downloaded [here](https://drive.google.com/drive/folders/1M9gIOFjK4wBk7DEj9MrwAlV8OHEgS2BA?usp=sharing).**

Training datasets include 2 Tonsil samples, 1 GBM sample, 1 Breast sample, and 2 Mesothelioma samples.

ST datasets include 5 samples, which are 2 Breast samples, a Melanoma sample, a Ovarian sample, and a Prostate sample.


## Prerequisites
The code runs on Python 3.10 and above.

Please find the requirements.txt file in this repository for the list of required packages. You can install them using pip:

```bash
pip install -r requirements.txt
```
****

## Run

Training notebook will lead you through the training process of DGAT model. Before running the notebook, please make sure you have downloaded the training datasets and put them in the `DGAT_datasets` folder. Detailed instructions are provided in the notebook.

## Prediction
To run the prediction, you can follow the Predict_ST notebook. This notebook will guide you through the process of using the pre-trained DGAT model to predict gene expression from Spatial Transcriptomics data. Please make sure you have the pre-trained model or the model you trained in the previous step. The ST datasets should also be placed in the `ST_data` folder.

## Key functions

### Train

**preprocess_train_list(adata_list, pdata_list)**
: This function performs quality control and normalization on the training datasets, finally it will return **common_genes** and **common_proteins** for you to check. **adata_list** should be the list of ST datasets, and **pdata_list** should be the list of protein datasets (all the datasets should be in the form of `anndata.AnnData`).

**train(train_adata_list, train_pdata_list, processed_data_dir)**
: This function trains the DGAT model using the preprocessed training datasets. It takes in the list of processed ST datasets and protein datasets, and the directory to save the processed graph data (for more efficient reproduction).

### Predict

**preprocess_ST(adata)**
: This function performs quality control and normalization on a Spatial Transcriptomics dataset in the form of `anndata.AnnData`. It returns a preprocessed one.

**fill_genes(test_adata, common_gene)**
: This function fills in the missing genes in the test dataset with zeros, **ensuring that the test dataset has the same genes as the training dataset**. It takes in the test dataset and the list of common genes.

**protein_predict(adata, common_gene, common_protein, model_save_dir, pyg_data_dir)**
: This function predicts the protein expression from the Spatial Transcriptomics dataset using the pre-trained DGAT model located in `model_save_dir`. It takes in the preprocessed test dataset, the list of common genes, the list of common proteins, and the directory for saving the graph data. It returns an `anndata.AnnData` object, in which the adata.X is the predicted protein expression, and the adata.var_names are the common proteins.
