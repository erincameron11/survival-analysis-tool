# Survival Analysis Tool

## Datasets
Datasets were obtained from `UCSC Xena` GDC Pan-Cancer (PANCAN) cohort:
* RNA matrix & ID/Gene Mapping matrix: [here](https://xenabrowser.net/datapages/?dataset=GDC-PANCAN.htseq_fpkm-uq.tsv&host=https%3A%2F%2Fgdc.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443)
* Phenotype matrix: [here](https://xenabrowser.net/datapages/?dataset=GDC-PANCAN.basic_phenotype.tsv&host=https%3A%2F%2Fgdc.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443)
* Survival matrix: [here](https://xenabrowser.net/datapages/?dataset=GDC-PANCAN.survival.tsv&host=https%3A%2F%2Fgdc.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443)

## Packages
* `pip install streamlit` -- For the UI of the application
* `pip install pandas` -- For data analysis
* `pip install numpy` -- For scientific computing
* `pip install gseapy` -- For the GSVA calculation
* `pip install kaplanmeier` -- For the creation of a Kaplan Meier plot
* `pip install matplotlib` -- For the output of the Kaplan Meier plot

## Development
Data was originally transformed from `csv` files to `parquet` files for faster overall application runtime. In addition to this, Streamlit's `@st.cache_data` decorator was used in order to cache the data after the first data read-in. 

## Run - Developer
1. Open Terminal and navigate to project folder
2. Type `streamlit run SurvivalAnalysisTool.py` to open the streamlit app in the browser
3. Use preferred editing software to edit the file and see changes reflected in the streamlit app browser