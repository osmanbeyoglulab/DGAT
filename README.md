## Introduction
Spatial transcriptomics (ST) technologies provide genome-wide mRNA profiles in tissue context 
but lack direct protein-level measurements, which are critical for interpreting cellular function 
and microenvironmental organization. We present DGAT (Dual-Graph Attention Network), a 
deep learning framework that imputes spatial protein expression from transcriptomics-only ST 
data by learning RNA–protein relationships from spatial CITE-seq datasets. DGAT constructs 
heterogeneous graphs integrating transcriptomic, proteomic, and spatial information, encoded 
using graph attention networks. Task-specific decoders reconstruct mRNA and predict protein 
abundance from a shared latent representation. Benchmarking across public and in-house 
datasets—including tonsil, breast cancer, glioblastoma, and malignant mesothelioma—
demonstrates that DGAT outperforms existing methods in protein imputation accuracy. Applied 
to ST datasets lacking protein measurements, DGAT reveals spatially distinct cell states, 
immune phenotypes, and tissue architectures not evident from transcriptomics alone. DGAT 
enables proteome-level insights from transcriptomics-only data, bridging a critical gap in spatial 
omics and enhancing functional interpretation in cancer, immunology, and precision medicine.

![img.png](imgs/img.png)
## Data 
**The training datasets for DGAT model, pre-trained model and Spatial Transcriptomics for predictions can be downloaded [here](https://drive.google.com/drive/folders/1M9gIOFjK4wBk7DEj9MrwAlV8OHEgS2BA?usp=sharing).**

## Installation
The code runs on Python 3.11.

You can install the required packages using pip:

  - For CUDA 
    ```
    pip install -r requirements_CUDA.txt
    ```
  - For already installed torch 
    ```
    pip install -r requirements_torch_ready.txt
    ```
  - For CPU-only
    ```
    pip install -r requirements_CPU.txt
    ```



## Usage

### Train

[Training notebook](Demo_Train.ipynb) will lead you through the training process of DGAT model. Before running the notebook, please make sure you have downloaded the training datasets and put them in the `DGAT_datasets` folder. Detailed instructions are provided in the notebook.

### Prediction
To run the prediction, you can follow the [Predict_ST](Predict_ST.ipynb) notebook. This notebook will guide you through the process of using the pre-trained DGAT model to predict gene expression from Spatial Transcriptomics data. Please make sure you have the pre-trained model or the model you trained in the previous step. The ST datasets should also be placed in the `ST_data` folder.
