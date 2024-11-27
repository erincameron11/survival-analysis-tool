# TCGA SIGvival
A TCGA survival analysis tool, using `ssGSEA` to calculate gene signature scores, and `kaplanmeier` python package to plot kaplan meier plots.

---

## Datasets
Datasets were obtained from `UCSC Xena` GDC Pan-Cancer (PANCAN) cohort:
* RNA matrix & ID/Gene Mapping matrix: [here](https://xenabrowser.net/datapages/?dataset=GDC-PANCAN.htseq_fpkm-uq.tsv&host=https%3A%2F%2Fgdc.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443)
* Phenotype matrix: [here](https://xenabrowser.net/datapages/?dataset=GDC-PANCAN.basic_phenotype.tsv&host=https%3A%2F%2Fgdc.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443)
* Survival matrix: [here](https://xenabrowser.net/datapages/?dataset=GDC-PANCAN.survival.tsv&host=https%3A%2F%2Fgdc.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443)

## Introduction
The `TCGA SIGvival` tool is a Python-based web application designed to make cancer data exploration accessible and interactive. Built with Streamlit to enhance the user experience, and deployed on Streamlit Community Cloud, the app allows users to perform survival analysis using datasets from The Cancer Genome Atlas (TCGA). By entering specific gene names and cancer types, users can calculate ssGSEA scores and visualize the impact of different genetic factors on patient survival through Kaplan-Meier plots. This streamlined approach enables researchers and data enthusiasts alike to gain insights into cancer outcomes without the need for extensive coding or data preprocessing.

Visit the `TCGA SIGvival` application to explore it yourself at [https://tcga-sigvival.streamlit.app](https://tcga-sigvival.streamlit.app).   
Visit the `TCGA SIGvival: Survival Analysis - Python` blog post at [https://erincameron11.github.io/posts/2024/11/tcga-sigvival/](https://erincameron11.github.io/posts/2024/11/tcga-sigvival/).

## Objective
Develop a Python application, using Streamlit Community Cloud to perform interactive survival analysis on TCGA cancer datasets.

Behind the scenes, the application performs the following tasks:
1. Loads all relevant datasets to extract gene names, and cancer type data. Streamlit's `@st.cache_data` decorator was used in order to cache the data after the first data read-in.
2. Uses the user-entered data to calculate ssGSEA scores.
3. Plots the results on a Kaplan Meier plot and outputs the plot.

## Implementation
### Files:
The TCGA Survival Analysis Tool is composed of several key Python files, each serving a distinct function in the process. Below is an overview of the contents of each file and its role in the project:
* **DataPreprocessing.ipynb**: This notebook handles the preprocessing of TCGA data, including formatting, cleanup, and gene name mapping. It ensures that the data is ready for downstream analysis by converting it into a usable format and resolving any inconsistencies in gene annotations.
* **DataAnalysis.ipynb**: This notebook focuses on testing and validating the core calculations for the app. It contains code for performing GSVA (Gene Set Variation Analysis) or ssGSEA (single-sample Gene Set Enrichment Analysis), as well as generating Kaplan-Meier plots to visualize survival data based on gene expression signatures.
* **SurvivalAnalysisTool.py**: This is the main application file, where the Streamlit UI/UX is structured. It defines the layout, handles user inputs, and integrates the functionality provided by the other modules. This is the heart of the app, allowing users to interact with the data and visualize results.
* **data.py**: This file is responsible for loading all the necessary datasets into the application. It provides the data infrastructure for the app, ensuring that the correct files are accessible and properly formatted for analysis.
* **helpers.py**: This utility file handles various functions to support the app’s operation. It includes form validation, submission handling, garbage collection for unused variables, and the creation of the RNA expression dataframe from the selected cancer types. Additionally, it calculates ssGSEA scores, plots Kaplan-Meier curves, and provides functionality for users to download the results.
* **styling.py**: This file implements custom JavaScript and CSS to enhance the user interface of the app. It includes styling for various elements, as well as functionality for auto-scrolling the page in response to user actions, improving the overall user experience.
* **./data folder**: This folder includes all TCGA expression matrices, separated out by cancer types.

### Entering Data into the App
To use the TCGA Survival Analysis Tool, follow these steps on the app's interface:
1. **Enter Signature Name**: Provide the custom name of the gene signature you wish to analyze.
2. **Select Gene Names**: Choose the specific genes you want to include in the analysis. These can be selected from the list of available genes based on the signature.
3. **Select Cancer Types**: Pick one or more cancer types from the available options.
4. **Select Cut-Point**: Choose the cut-point for the survival analysis. This value will be used to categorize patients into high or low expression groups for the Kaplan-Meier plot.
5. **Click "Create KM Plot"**: Once all the form fields are filled out, press the button to generate the Kaplan-Meier plot, which visualizes the survival curves based on the selected gene signature and cancer types. 

## Technology
* Python
* Jupyter Notebook
* Streamlit

## Packages
* `streamlit` -- For the UI of the application
* `pandas` -- For data analysis
* `numpy` -- For scientific computing
* `gseapy` -- For the ssGSEA calculations
* `kaplanmeier` -- For the creation of a Kaplan Meier plot
* `matplotlib` -- For the output of the Kaplan Meier plot
* `datetime` -- For file naming convention for exports
* `threading` -- For accelerating the ssGSEA calculation
* `statsmodels.api` -- For hazard ratio calculations 
* `os` -- For KM plot downloading
* `gc` -- For garbage collection of unused objects
* `io` -- For data export
* `zipfile` -- To zip files for data exports

## Run - Developer
1. Open Terminal and navigate to project folder
2. Type `streamlit run SurvivalAnalysisTool.py` to open the streamlit app in the browser
3. Use preferred editing software to edit the file and see changes reflected in the streamlit app browser