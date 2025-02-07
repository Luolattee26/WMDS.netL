# WMDS.netL: Advanced Cancer-Driving LncRNA Identification

## Table of Contents

- [Introduction](#Introduction)
- [About](#About)
- [Usage](#Usage)
  - [Installation](#Installation)
  - [Run](#Run)
- [Analysis](#Analysis)
- [System Information](#System-Information)
- [Contact](#Contact)

## Introduction

`WMDS.net` is an algorithm based on network control theory for identifying cancer driver genes. Compared with other methods and traditional differential gene statistical tests, `WMDS.net` offers higher accuracy, thereby reducing false positives (https://github.com/chaofen123/WMDS.net, https://doi.org/10.1093/bioinformatics/btad071). `WMDS.netL` is an improved and optimized version of `WMDS.net`, focusing specifically on the identification of cancer-driving lncRNAs during tumorigenesis and progression.
![Workflow of WMDS.netL](workflow.png)

## About

This repository includes the deployment code for `WMDS.netL` and related integration analysis codes. For the data used, if the file size meets GitHub's upload restrictions, it will also be included here (for files exceeding the size limit, acquisition methods will be provided). You can reproduce the results presented in our paper (to be published) using these codes. All codes are organized according to the sequence of figures in the paper, with brief comments at the beginning of each code file explaining its purpose and the final output.

## Usage

The WMDS.netL algorithm is built via **MATLAB** and you can find it in `./code/WMDS.netL_algorithm/` where you can find the relevant source code as well as a priori data.

### Installation

```
git clone git@github.com:Luolattee26/WMDS.netL.git
cd WMDS.netL
conda env create -f environment.yml
conda activate lncRNA_W
```

### Run  

- Before running the algorithm, ensure that **MATLAB** is installed on your system. The `WMDS.netL` algorithm requires two input files:  
  - A normal expression matrix (`TYPEnormal.txt`)  
  - A tumor expression matrix (`TYPEtumor.txt`)  
  These files should be placed in the `./data/` directory, where each row represents a gene and each column represents a sample.  

- Once the input files are prepared and the runtime environment is set up, update the `cancer_type` variable in the `WMDS.netL` code located at `/code/WMDS.netL_algorithm/WMDS.netL.m`.  
  - Set `cancer_type` to match the `TYPE` specified in your input file names.  
  - Additionally, in the `for` loop (`for o=1:14`), adjust the values as needed based on the specific cancer type you are analyzing.  

- For a better understanding of the required input data format, you may refer to TCGA data available on the Xena platform.  
  ![Input format](input_format.jpg)  

To run `WMDS.netL`, use the following command in **BASH**:  
```bash
matlab -nodisplay -nosplash -r code/WMDS.netL_algorithm/WMDS.netL
```  
or directly run in **MATLAB**:
```MATLAB
code/WMDS.netL_algorithm/WMDS.netL
```  




## Analysis  

If you are interested in our analytical approach or would like to reproduce our results, please refer to the `./code/` directory.  

- Each code file is organized in the order corresponding to the figures in the manuscript, with filenames indicating their respective figures.  
- Some original data files have not been uploaded due to size limitations. However, all data are publicly available, and you can refer to our manuscript for access details.  
- Certain figures or panels require multiple scripts for analysis and visualization. The execution order of these scripts is indicated in their filenames.  

### Example: Figure 3D  
For Figure 3D, we used the phastCon score to assess the conservation of different lncRNA transcripts. This analysis involves the following steps:  

1. **Prepare input files** using `./code/fig.3/fig3.d_0_gtf2bed.sh`  
2. **Perform conservation analysis** with `./code/fig.3/fig3.d_1_get_phastCon_res.py`  
3. **Generate plots** using `./code/fig.3/fig3.d&e_plot_Con_TS.R`  

- The numerical prefix (e.g., `0`, `1`) in the filenames denotes the order in which the scripts should be executed.  
- Similar multi-script workflows apply to the following analyses: *Figure 3F, 4A, 4B, 6A, 6B*.  




## System Information  

- The `WMDS.netL` algorithm was implemented and tested on a PC with **Intel 4-core CPUs (3.40 GHz × 4) and 24 GB of RAM**. Systems with higher hardware specifications can efficiently run `WMDS.netL` and the other analyses provided in this repository.  

- We strongly recommend conducting the subsequent analyses on **Linux, macOS, or WSL (Windows Subsystem for Linux)**. This is because the `pyBigWig` package, which is required for conservation analysis, is only supported in Unix-like environments.  

- For detailed information about the computational environment, refer to the configuration files included in this repository:  
  - `./sessionInfo.txt` contains the session details for the R analysis scripts.  
  - `./environment.yml` provides the Conda virtual environment configuration.  




## Contact

We welcome anyone to use `WMDS.netL` for academic exploration in cancer biology, please cite our latest publication (to be published).

If you have any questions or would like to discuss ideas, feel free to contact us at: luo_itm@zju.edu.cn.
