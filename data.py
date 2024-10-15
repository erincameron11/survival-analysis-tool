# ------------------------------------ IMPORTS ------------------------------------
import streamlit as st
import pandas as pd # for data set analysis and manipulation
from helpers import *



# ------------------------------------ DATA ------------------------------------
@st.cache_data
def load_data():
    """
    Loads gene names and cancer types, as well as survival, and phenotype data files. Uses st.cache_data decorator to cache the DataFrames.

    Parameters
    ----------
    None

    Returns
    -------
    gene_names : list (str)
        A list of all gene names from the RNA dataset.
    cancer_types : list (str)
        A list of all cancer types from the phenotype dataset.
    phenotype_df : pandas DataFrame
        Phenotype DataFrame filtered for common samples, and reordered to RNA ordering.
    survival_df : pandas DataFrame
        Survival DataFrame filtered for common samples, and reordered to RNA ordering.
    """
    # Load the smallest cancer type dataset to gather the gene names
    gene_names_df = pd.read_parquet('./data/GDC-PANCAN.htseq_fpkm-uq_TCGA-CHOL.parquet')
    gene_names = gene_names_df.index.tolist()
    
    # Load the phenotype dataset to gather the cancer types
    phenotype_df = pd.read_parquet('./data/GDC-PANCAN.basic_phenotype_processed.parquet')
    cancer_types = phenotype_df['project_id'].unique()

    # Load the survival dataset
    survival_df = pd.read_parquet('./data/GDC-PANCAN.survival_processed.parquet')

    # Garbage collection of unused objects 
    garbage_collection(gene_names_df)

    return gene_names, cancer_types, phenotype_df, survival_df
