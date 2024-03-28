<div align="center">

# Map of epigenetic age acceleration: a worldwide meta-analysis

</div>

<br>

## Description
Repository with scripts used in the generation of the manuscript [Map of epigenetic age acceleration: a worldwide meta-analysis](https://www.biorxiv.org/content/10.1101/2024.03.17.585398).
> **Note**: _No new or unpublished methods were used in the study._

## Abstract
We present a systematic meta-analysis of epigenetic age acceleration based on by far the largest collection publicly
available DNA methylation data for healthy samples (93 datasets, 23K samples), focusing on the geographic (25 countries) 
and ethnic (31 ethnicities) aspects around the world. We employed the most popular epigenetic tools for assessing age 
acceleration and examined their quality metrics and ability to extrapolate to epigenetic data from different tissue 
types and age ranges different from the training data of these models. In most cases, the models proved to be inconsistent 
with each other and showed different signs of age acceleration, with the PhenoAge model tending to systematically underestimate 
and different versions of the GrimAge model tending to systematically overestimate the age prediction of healthy subjects. 
Referring to data availability and consistency, most countries and populations are still not represented in GEO, moreover, 
different datasets use different criteria for determining healthy controls. Because of this, it is difficult to fully isolate 
the contribution of "geography/environment", "ethnicity" and "healthiness" to epigenetic age acceleration.  
Among the explored metrics, only the DunedinPACE, which measures aging rate, appears to adequately reflect the standard 
of living and socioeconomic indicators in countries, although it has a limited application to blood methylation data only. 
Invariably, by epigenetic age acceleration, males age faster than females in most of the studied countries and populations.

## Project Structure

```
├── notebooks              <- Jupyter notebooks and R scripts have been used for DNAm data processing and figure generation.
│   ├── GPL13534                 <- Scripts for GPL13534 Illumina standard.
│   │   ├── GSE...                     <- Script for certain Dataset.
│   │   └── ...  
│   ├── GSE16304                 <- Scripts for GSE16304 Illumina standard.
│   │   ├── GSE...                     <- Script for certain Dataset.
│   │   └── ...  
│   ├── GPL21145                 <- Scripts for GPL21145 Illumina standard.
│   │   ├── GSE...                     <- Script for certain Dataset.
│   │   └── ...  
│   ├── GPL23976                 <- Scripts for GPL23976 Illumina standard.
│   │   ├── GSE...                     <- Script for certain Dataset.
│   │   └── ...  
│   │
│   ├── plots_epi_age.ipynb      <- Script for plotting almost all figures related to epigenetic age acceleration.
│   ├── plots_geomap.ipynb       <- Script for plotting geographical maps.
│   ├── plots_metrics.ipynb      <- Script with Machine Learning metrics analysis.
│   └── plots_workflow.py        <- Script with workflow plots for Figure 1.
│
├── src                    <- Source code for auxiliary functions
│
├── .gitignore                <- List of files ignored by git
├── .project-root             <- File for inferring the position of project root directory
├── requirements.txt          <- File for installing python dependencies
└── LICENSE                   <- MIT license
└── README.md
```

## Project Environment

Install python dependencies:

```bash
# clone project
git clone https://github.com/GillianGrayson/MetaEpiAge
cd MetaEpiAge

# [OPTIONAL] create conda environment
conda create -n myenv python=3.9
conda activate myenv

# install requirements
pip install -r requirements.txt
```

Install required R packages:

```r
if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install("ChAMP")
BiocManager::install("methylGSA")
BiocManager::install("preprocessCore")
install.packages("devtools")
install.packages("splitstackshape")
install.packages("reticulate")
devtools::install_github("danbelsky/DunedinPACE")
devtools::install_github("https://github.com/regRCPqn/regRCPqn")
```

For using [PC clocks](https://www.nature.com/articles/s43587-022-00248-2) you need to setup path to 3 files from [Original Repository](https://github.com/MorganLevineLab/PC-Clocks):
1) `run_calcPCClocks.R`
2) `run_calcPCClocks_Accel.R`
3) `CalcAllPCClocks.RData` (very big file, but it is nesessary)
4) Apply changes from [This Issue](https://github.com/MorganLevineLab/PC-Clocks/issues/10)

## License

MetaEpiAge is licensed under the MIT License.

```
The MIT License (MIT)

Copyright (c) 2024 Alena Kalyakulina, Igor Yusipov

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```