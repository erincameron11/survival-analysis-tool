# Import statements
import streamlit as st # for UI
import pandas as pd # for data set analysis and manipulation
import kaplanmeier as km
import matplotlib.pyplot as plt # TESTING for image export of KM plots
from datetime import datetime # for file naming convention for exports
import numpy as np # TESTING for st.pyplot
import streamlit.components.v1 as components # for KM plot page anchor
import time # for page anchor scrolling
import gseapy as gp # for GSVA calculation
import threading # for accelerating the GSVA calculation
import os # for KM plot downloading
from pathlib import Path # for KM plot downloading



# ------------------------------------ DATA ------------------------------------
@st.cache_data
def load_data(rna_filename, gene_mapping_filename, survival_filename, phenotype_filename):
    """
    Loads all RNA, Gene Mapping, Survival, and Phenotype data files. Uses st.cache_data decorator to cache the dataframes.

    Parameters
    ----------
    rna_filename : str
        RNA filename string.
    gene_mapping_filename : str
        Gene mapping filename string.
    survival_filename : str
        Survival filename string.
    phenotype_filename : str
        Phenotype filename string.

    Returns
    -------
    merged_trimmed_filtered_df : pandas DataFrame
        RNA dataframe with mapped gene names, and filtered for common samples
    survival_filtered_ordered_df : pandas DataFrame
        Survival dataframe filtered for common samples, and reordered to RNA ordering
    phenotype_filtered_ordered_df : pandas DataFrame
        Phenotype dataframe filtered for common samples, and reordered to RNA ordering
    """
    # Load in the RNA matrix 
    df = pd.read_parquet(rna_filename)
    
    # Load in the ID/Gene Mapping file and create a merged dataframe to map the RNA IDs to the gene names
    mapping = pd.read_csv(gene_mapping_filename, sep='\t')
    mapping.rename(columns={'id': 'xena_sample'}, inplace=True) # rename column to merge on
    merged_df = pd.merge(mapping, df, on='xena_sample', how='outer', indicator=True) # merge the dataframes to map gene names

    # Ensure the expression dataframe is in the format: indexed on gene names column labels as sample ids
    merged_trimmed_df = merged_df
    merged_trimmed_df.drop(columns=['xena_sample', 'chrom', 'chromStart', 'chromEnd', 'strand', '_merge'], axis=1, inplace=True)
    merged_trimmed_df.set_index('gene', inplace=True)
    
    # Load in survival and phenotype dataframes
    survival_df = pd.read_parquet(survival_filename)
    phenotype_df = pd.read_parquet(phenotype_filename)

    # Create a list of all samples in each dataframe
    samples_rna = list(merged_trimmed_df.columns)
    samples_pheno = list(phenotype_df['sample'].values)
    samples_survival = list(survival_df['sample'].values)
    
    # Find all common samples in all three lists
    common_samples = list(set(samples_rna) & set(samples_pheno) & set(samples_survival))
    
    # Subset and reorder all three datasets by common_samples
    merged_trimmed_filtered_df = merged_trimmed_df[common_samples] # Filter merged_df by columns in common_samples
    phenotype_filtered_df = phenotype_df[phenotype_df['sample'].isin(common_samples)]
    survival_filtered_df = survival_df[survival_df['sample'].isin(common_samples)]
    
    # Reorder phenotype and survival dataframes to match the columns of the rna matrix
    column_order = list(merged_trimmed_filtered_df.columns)
    phenotype_filtered_ordered_df = phenotype_filtered_df.set_index('sample').loc[column_order].reset_index()
    survival_filtered_ordered_df = survival_filtered_df.set_index('sample').loc[column_order].reset_index()
    
    return merged_trimmed_filtered_df, survival_filtered_ordered_df, phenotype_filtered_ordered_df



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
    signature_name = st.session_state.get('signature_name', '')
    genes_entered = st.session_state.get('genes_entered', '')
    cancer_types_entered = st.session_state.get('cancer_types_entered', '')
    cut_point_entered = st.session_state.get('cut_point_entered', '')
    
    # If all form fields filled out, return True, else False
    if signature_name and genes_entered and cancer_types_entered and cut_point_entered:
        return True
    else:
        return False


def calculate_gsva(df, phenotype_df):
    """
    Calculates GSVA scores using GSEAPY ssGSEA.

    Parameters
    ----------
    df : pandas DataFrame
        RNA dataframe.

    Returns
    -------
    scores : gseapy.ssgsea.SingleSampleGSEA
        ssGSEA dataframe output with score values.
    """
    # Use state sessions to get form values
    signature_name = st.session_state.get('signature_name', '')
    genes_entered = st.session_state.get('genes_entered', '')
    cancer_types_entered = st.session_state.get('cancer_types_entered', '')
    cut_point_entered = st.session_state.get('cut_point_entered', '')

    # Create a dictionary of signature and gene names
    signature = {signature_name: genes_entered}
    
    # Determine the number of threads to run the calculations on
    n_threads=threading.active_count()-1

    # Get a list of samples for cancer types
    cancer_type_samples = list(phenotype_df.query('project_id in @cancer_types_entered')['sample'])
    # Subset the original RNA count matrix to only have samples in selected cancer types
    subset_counts = df.loc[: , cancer_type_samples]
    
    # Calculate the GSVA scores
    scores = gp.ssgsea(data=subset_counts, gene_sets=signature, outdir=None, 
               sample_norm_method='rank', threads=n_threads, min_size=0,
               verbose=True)
    return scores


def create_km_plot(path):
    """
    Creates a Kaplan Meier plot output and saves to file.

    Parameters
    ----------
    path : str
        RNA filename string.

    Returns
    -------
    None
    """
    # Create and show KM plot
    # TESTING -- Example dataframe from kaplanmeier package
    example_df = km.example_data()
    # Kaplan-Meier fit
    time_event = example_df['time']
    censoring = example_df['Died']
    group = example_df['group']
    # Compute survival Kaplan-Meier estimates
    results = km.fit(time_event, censoring, group)
    # Save the plot as an image, don't display
    km.plot(results, savepath=path, visible=False, dpi=200)


def download_output():
    """
    Downloads GSVA data to CSV file and KM plot to PNG on the users' local machine.

    Parameters
    ----------
    None

    Returns
    -------
    None
    """
    # Get the current date and time for file naming
    today = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    # Get the Downloads folder path
    downloads_folder = str(Path.home() / "Downloads")
    
    # GSVA export into CSV format
    # TESTING -- gsva calculations for testing purposes
    gsva_testing = pd.read_parquet('./data/gsva_scores.parquet')
    # Create the full path for the GSVA scores
    gsva_file_path = os.path.join(downloads_folder, f'gsva_scores_{today}.csv')
    # Output to CSV to the specified filepath
    gsva_testing.to_csv(gsva_file_path, index=False)

    # KM plot export
    # Create the full path for the KM plot
    km_file_path = os.path.join(downloads_folder, f'km_plot_{today}.png')
    # Create KM plot and save to the file path
    create_km_plot(km_file_path)


def block_form_submit():
    """
    Blocks the Streamlit form from submitting on Enter press with text_input components, using JavaScript injection.

    Parameters
    ----------
    None

    Returns
    -------
    None
    """
    components.html("""
    <script>
        const inputs = window.parent.document.querySelectorAll('input');
        inputs.forEach(input => {
            input.addEventListener('keydown', function(event) {
                if (event.key === 'Enter') {
                    event.preventDefault();
                }
            });
        });
        </script>
    """, height=0)


def auto_scroll():
    """
    Auto-scrolls user screen down 500 pixels on each method call.

    Parameters
    ----------
    None

    Returns
    -------
    None
    """
    # Define JavaScript code for auto-scroll to the results section
    st.components.v1.html("""
        <script>
            console.log(window.parent.document.querySelector(".main"));
            window.parent.document.querySelector(".main").scrollTo({top: 500, behavior: 'smooth'});
        </script>""", height=0)



# ------------------------------------ STYLING FUNCTIONS ------------------------------------
# Function to alter CSS styling for multiselect, text input, and buttons
def custom_css():
    """
    Applies all custom CSS to the Streamlit application.

    Parameters
    ----------
    None

    Returns
    -------
    None
    """
    st.markdown("""
       <style>
        /* Multiselect initial border colour */
        div[data-baseweb="select"] > div {
            border: 0em solid white !important;
            border-radius: 4px !important;
            box-shadow: 0 0 0 1px white !important; /* Maintain custom box-shadow */
        }
        /* Multiselect outline colour when focused */
        div[data-baseweb="select"] > div:focus-within {
            border: 0em solid #c4c1c1;
            box-shadow: 0 0 0 1px #c4c1c1 !important;  
            border-color: #c4c1c1 !important;
        }
        /* Multiselect tag styling */
        .stMultiSelect div[data-baseweb="select"] span[data-baseweb="tag"] {
            background-color: #4A9661 !important;
            color: white !important;
        }
        /* Text Input outline styling */
        .stTextInput > div[class]:focus-within {
            border-bottom-color: #c4c1c1 !important;
            border-top-color: #c4c1c1 !important;
            border-left-color: #c4c1c1 !important;
            border-right-color: #c4c1c1 !important;
        }
        [data-testid="InputInstructions"] { 
            display: None;
        }
        /* Text Input placeholder text opacity */
        .stTextInput > div > div > input::placeholder {
            opacity: 1.0;
        }
        /* Text Input border radius */
        .stTextInput > div {
            border-radius: 4px;
        }
        /* Button initial style */
        div.stButton > button, div.stFormSubmitButton > button {
            background-color: white !important;
            border-color: #D5D6D8 !important;
            color: #31333F !important;
        }
        /* Button hover effect */
        div.stButton > button:hover, div.stFormSubmitButton > button:hover {
            background-color: #1f77b4 !important;
            border-color: #1f77b4 !important;
            color: white !important;
        }
        /* Button active effect */
        div.stButton > button:active, div.stFormSubmitButton > button:active {
            background-color: #5a9bd4 !important;
            border-color: #5a9bd4 !important;
            color: white !important;
        }
        </style>
    """, unsafe_allow_html=True)


# ------------------------------------ APP ------------------------------------
def main():
    """
    Main method to run the application.

    Parameters
    ----------
    None

    Returns
    -------
    None
    """
    # App title
    st.title("SURVIVAL ANALYSIS")
    st.divider()
    
    # Create a field for informational text
    # st.write("Describe KM plots, GSVA and what to do on this webpage...")

    # Call the load data method
    df, survival_df, phenotype_df = load_data('./data/GDC-PANCAN.htseq_fpkm-uq.parquet', 
                                              './data/gencode.v22.annotation.gene.probeMap',
                                              './data/GDC-PANCAN.survival.parquet', 
                                              './data/GDC-PANCAN.basic_phenotype.parquet',
                                             )
    # Locate all gene names in a list
    gene_names = df.index.unique()
    
    # Create a form for data input
    with st.form("km_plot_form", clear_on_submit=False):
        # Create a placeholder for form validation error correction
        form_validation_placeholder = st.empty()

        # Text input field for custom signature name
        signature_name = st.text_input("Signature Name:", value="", placeholder="Enter signature name", key='signature_name')
        
        # Multiselect element for gene selection
        genes_entered = st.multiselect(
            "Gene Names:",
            gene_names,
            placeholder="Enter gene names",
            key='genes_entered',
        )
        
        # Dropdown for cancer type
        cancer_types = phenotype_df['project_id'].unique()
        # TODO: add a PAN-Cancer option for all cancers?
        cancer_types_entered = st.multiselect(
            "Cancer Type:",
            (cancer_types),
            placeholder="Select cancer type",
            key='cancer_types_entered',
            # on_change=calculate_gsva,
        )
    
        # Dropdown for cut-point
        cut_point_entered = st.selectbox(
            "Cut-Point:",
            ("Median", "Tertile - ALL", "Tertile - Top and Bottom only", "Quartile - ALL", "Quartile - Top and Bottom only"),
            index=None,
            placeholder="Select cut-point",
            key='cut_point_entered',
        )
        
        # Submit button
        submit_button = st.form_submit_button(":chart_with_downwards_trend: Create KM Plot", on_click=handle_submit)
        
        # If the submit button was clicked, check that all fields are filled in
        if submit_button:
            # If not all fields are filled in
            if not validate_form():
                # Display red text form validation
                with form_validation_placeholder:
                    st.markdown("""
                        <p style='color:#cc0000; text-align:center;'>Please fill out all fields</p>
                    """, unsafe_allow_html=True)
            else:
                # Auto scroll to GSVA calculation info message
                auto_scroll()

    # Block the form from submitting on Enter press with text_input (built-in streamlit functionality)
    block_form_submit()

    # Apply CSS for custom styling
    custom_css()

    # If the submit button was pressed and submitted successfully
    if st.session_state.get('form_submitted', False):
        # Calculate GSVA
        gsva_info = st.info('Calculating GSVA scores...', icon="🔄")
        gsva = calculate_gsva(df, phenotype_df)
        gsva_info.empty()
        
        # Create the kaplan meier results
        create_km_plot('km_plot.png')

        # Scroll down once calculations complete
        auto_scroll()
        
        # Display the results inside a container
        with st.container(border=True):
            # Display results subheader
            st.subheader("Results")
            
            # Use state sessions to get form values
            signature_name = st.session_state.get('signature_name', '')
            genes_entered = st.session_state.get('genes_entered', '')
            cancer_types_entered = st.session_state.get('cancer_types_entered', '')
            cut_point_entered = st.session_state.get('cut_point_entered', '')
            
            # Display the entered form values for user
            genes_entered_str = ", ".join(genes_entered)
            cancer_types_entered_str = ", ".join(cancer_types_entered)
            st.write("**Signature Name**: ", signature_name, "  \n**Gene Names**: ", genes_entered_str, "  \n**Cancer Types**: ", cancer_types_entered_str, "  \n**Cut-point**: ", cut_point_entered)
            st.divider()

            # Create placeholders to hold the GSVA, KM plot and download button content
            # gsva_placeholder = st.empty()
            km_plot_placeholder = st.empty()
            download_results_placeholder = st.empty()
            # with gsva_placeholder:
                # Display GSVA output
                # st.write(gsva.res2d.head())
            with km_plot_placeholder:
                # Display the KM plot image created
                st.image("km_plot.png")
            with download_results_placeholder:
                st.button(":arrow_down: Download Results", on_click=download_output)
    

# ------------------------------------ RUN THE APP ------------------------------------
if __name__ == "__main__":
    main()