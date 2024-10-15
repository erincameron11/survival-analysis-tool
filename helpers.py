# ------------------------------------ IMPORTS ------------------------------------
import streamlit as st # for UI
import pandas as pd # for data set analysis and manipulation
from datetime import datetime # for file naming convention for exports
import numpy as np # for scientific calculations
import gseapy as gp # for ssGSEA calculation
import threading # for accelerating the ssGSEA calculation
import kaplanmeier as km # for kaplan meier plotting
import statsmodels.api as sm # for hazard ratio calculations 
import os # for KM plot downloading
from pathlib import Path # for KM plot downloading
import psutil # TESTING -- for memory logging
import gc # TESTING -- for garbage collection of unused objects



# ------------------------------------ HELPER FUNCTIONS ------------------------------------
def handle_submit():
    """
    Sets the session_state of the form to True once the form is validated.

    Parameters
    ----------
    None

    Returns
    -------
    None
    """
    # If the entire form is filled out, submit
    if validate_form():
        # Mark the form as submitted
        st.session_state.form_submitted = True


def validate_form():
    """
    Validates that all form fields are filled out with data.

    Parameters
    ----------
    None

    Returns
    -------
    bool
        True if the form is valid, False if the form is invalid.
    """
    # Use state sessions to get form values
    signature_name, genes_entered, cancer_types_entered, cut_point_entered = get_form_values()
    
    # If all form fields filled out, return True, else False
    if signature_name and genes_entered and cancer_types_entered and cut_point_entered:
        # Garbage collection of unused objects 
        garbage_collect_form_values(signature_name, genes_entered, cancer_types_entered, cut_point_entered)
        st.session_state.form_submitted = True
        return True
    else:
        # Garbage collection of unused objects 
        garbage_collect_form_values(signature_name, genes_entered, cancer_types_entered, cut_point_entered)
        return False


def get_form_values():
    """
    Obtains form values entered by user using Streamlit session_state.

    Parameters
    ----------
    None

    Returns
    -------
    signature_name : str
        Custom name of the signature entered by user.
    genes_entered : list (str)
        1 or more genes selected by the user.
    cancer_types_entered : list (str)
        1 or more cancer types selected by user.
    cut_point_entered : str
        Cut-point selected by user.
    """
    # Use state sessions to get form values
    signature_name = st.session_state.get('signature_name', '')
    genes_entered = st.session_state.get('genes_entered', '')
    cancer_types_entered = st.session_state.get('cancer_types_entered', '')
    cut_point_entered = st.session_state.get('cut_point_entered', '')
    return signature_name, genes_entered, cancer_types_entered, cut_point_entered


# TEST - to collect and remove unused objects and items
def garbage_collection(garbage):
    del garbage
    gc.collect()

def garbage_collect_form_values(signature_name, genes_entered, cancer_types_entered, cut_point_entered):
    # Garbage collection of unused objects 
    garbage_collection(signature_name)
    garbage_collection(genes_entered)
    garbage_collection(cancer_types_entered)
    garbage_collection(cut_point_entered)
    

def create_rna_dataframe(cancer_types_entered):
    """
    Generates an RNA DataFrame by reading in and concatenating datasets for the user-selected cancer types.

    Parameters
    ----------
    cancer_types_entered : list (str)
        A list of cancer types selected by the user.

    Returns
    -------
    df : pandas DataFrame
        The RNA DataFrame constructed using selected cancer types.
    """
    # Identify folder where the files are stored
    data_folder = './data/'
    # Define an empty list to hold all the loaded DataFrames
    df_list = []
    
    # Loop through each cancer type
    for cancer_type in cancer_types_entered:
        # Construct the file name pattern to look for
        file_name = f'GDC-PANCAN.htseq_fpkm-uq_{cancer_type}.parquet'
        file_path = os.path.join(data_folder, file_name)

        # TCGA-BRCA was separated into 2 separate files for file size considerations
        if cancer_type == 'TCGA-BRCA':
            # Read the two parquet files
            df = pd.read_parquet(f'./data/GDC-PANCAN.htseq_fpkm-uq_{cancer_type}_1.parquet')
            df_list.append(df)
            df = pd.read_parquet(f'./data/GDC-PANCAN.htseq_fpkm-uq_{cancer_type}_2.parquet')
            df_list.append(df)
        else:
            df = pd.read_parquet(f'./data/GDC-PANCAN.htseq_fpkm-uq_{cancer_type}.parquet')
            df_list.append(df)
            
    df = pd.concat(df_list, axis=1)

    # Garbage collection of unused objects 
    garbage_collection(df_list)

    return df


def calculate_ssgsea(df, phenotype_df):
    """
    Calculates ssGSEA scores using GSEAPY ssGSEA.

    Parameters
    ----------
    df : pandas DataFrame
        RNA DataFrame.
    phenotype_df : pandas DataFrame
        Phenotype DataFrame.

    Returns
    -------
    scores.res2d : pandas DataFrame
        ssGSEA DataFrame output with score values.
    """
    # Use state sessions to get form values
    signature_name, genes_entered, cancer_types_entered, cut_point_entered = get_form_values()

    # Create a dictionary of signature and gene names
    signature = {signature_name: genes_entered}
    
    # Determine the number of threads to run the calculations on
    n_threads=threading.active_count()-1
    
    # Calculate the ssGSEA scores
    scores = gp.ssgsea(data=df, gene_sets=signature, outdir=None, 
               sample_norm_method='rank', threads=n_threads, min_size=1,
               verbose=True)
    
    # Garbage collection of unused objects 
    garbage_collect_form_values(signature_name, genes_entered, cancer_types_entered, cut_point_entered)
    
    return scores.res2d


def create_km_plot(ssgsea_scores, survival_df):
    """
    Creates a Kaplan Meier plot output.

    Parameters
    ----------
    ssgsea_scores : pandas DataFrame
        ssGSEA DataFrame output with score values.
    survival_df : pandas DataFrame
        Survival DataFrame.

    Returns
    -------
    km_plot_figure : matplotlib.figure.Figure
        The Kaplan Meier plot figure object.
    """
    # Use state sessions to get form values
    signature_name, genes_entered, cancer_types_entered, cut_point_entered = get_form_values()
    cut_point_entered = cut_point_entered.lower()
    
    # Define the number of quantile cuts to make
    if 'median' in cut_point_entered:
        n = 2
        labels = ['Low: below median', 'High: above median']
    if 'tertile' in cut_point_entered:
        n = 3
        labels = ['Low: bottom tertile', 'Medium: middle tertile', 'High: top tertile']
    if 'quartile' in cut_point_entered:
        n = 4
        labels = ['Low: bottom quartile', 'Medium1: second quartile', 'Medium2: third quartile', 'High: top quartile']
    
    # Make the quantile cuts & label samples by the scoring grouping
    nes_scores = ssgsea_scores['NES']
    km_groups = pd.qcut(nes_scores, n, labels=labels)
    
    # Bind KM groups to survival dataframe by aligning the indices of both the Series and the DataFrame
    # Create a new dataframe by copying the original survival_df
    km_df = survival_df.copy()
    # Add the 'NES_group' column from km_groups to the new dataframe
    km_df['NES_group'] = km_groups
    
    # Drop any 'NES_group' null values
    km_df = km_df.dropna(subset=['NES_group'])
    
    # BUT the user might not want all groups (quantiles) on the plot (eg; top & bottom only)
    if 'top' in cut_point_entered:
        # Subset km_df to include only rows where 'NES_group' contains 'top' or 'bottom'
        km_df = km_df[km_df['NES_group'].str.contains('Top|Bottom', case=False)]

    # Plot the KM plot
    time_event = km_df['OS.time']
    censoring = km_df['OS'] # Alive / Dead
    y = km_df['NES_group']
    
    # Compute Survival
    results = km.fit(time_event, censoring, y)
    
    # Locate P value
    p_value = results['logrank_P']
    
    # Compute hazard ratio
    hazard_df = km_df.copy()
    hazard_df['NES_group'] = hazard_df['NES_group'].cat.codes
    cox_model = sm.PHReg(hazard_df['OS.time'], hazard_df[['NES_group']], status=hazard_df['OS'])
    hazard_results = cox_model.fit()
    # Locate the log hazard ratio (log HR)
    log_hazard_ratio = hazard_results.params[0]
    # Calculate the Hazard Ratio (HR)
    hazard_ratio = np.exp(log_hazard_ratio)
    
    # Plot with P value, hazard ratio, and signature name
    title = f'{signature_name}\nP={round(p_value, 4)}, HR={round(hazard_ratio, 4)}'
    km_plot = km.plot(results, title=title, dpi=300, visible=True)
    # Extract the figure and legend from the plot output
    km_plot_figure = km_plot[0]
    ax = km_plot[1]
    ax.legend(title='NES')
    # Adjust the margins and legend
    km_plot_figure.subplots_adjust(top=0.9, bottom=0.1, left=0.1, right=0.9)
    
    return km_plot_figure


def download_output(ssgsea_scores, km_plot_figure):
    """
    Downloads ssGSEA data to CSV file and KM plot to PNG on the users' local machine Downloads folder.

    Parameters
    ----------
    ssgsea_scores : pandas DataFrame
        ssGSEA DataFrame output with score values.
    km_plot_figure : matplotlib.figure.Figure
        The Kaplan Meier plot figure object.

    Returns
    -------
    None
    """
    # Get the current date and time for file naming
    today = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    # Get the Downloads folder path
    downloads_folder = str(Path.home() / "Downloads")
    
    # ssGSEA export into CSV format
    # Create the full path for the ssGSEA scores and output to CSV
    ssgsea_file_path = os.path.join(downloads_folder, f'ssgsea_scores_{today}.csv')
    ssgsea_scores.to_csv(ssgsea_file_path, index=False)

    # KM plot export
    # Create the full path for the KM plot and output file
    km_file_path = os.path.join(downloads_folder, f'km_plot_{today}.png')
    km_plot_figure.savefig(km_file_path, bbox_inches='tight')