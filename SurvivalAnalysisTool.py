# Import statements
import streamlit as st # for UI
import pandas as pd # for data set analysis and manipulation
import matplotlib.pyplot as plt # for KM plots
from datetime import datetime # for file naming convention for exports
import numpy as np # for scientific calculations
import streamlit.components.v1 as components # for KM plot page anchor
import gseapy as gp # for ssGSEA calculation
import threading # for accelerating the ssGSEA calculation
import kaplanmeier as km # for kaplan meier plotting
import statsmodels.api as sm # for hazard ratio calculations 
import os # for KM plot downloading



# ------------------------------------ DATA ------------------------------------
@st.cache_data
def load_data():
    """
    Loads all RNA, Gene Mapping, Survival, and Phenotype data files. Uses st.cache_data decorator to cache the dataframes.

    Parameters
    ----------
    None.

    Returns
    -------
    merged_trimmed_filtered_df : pandas DataFrame
        RNA dataframe with mapped gene names, and filtered for common samples
    survival_filtered_ordered_df : pandas DataFrame
        Survival dataframe filtered for common samples, and reordered to RNA ordering
    phenotype_filtered_ordered_df : pandas DataFrame
        Phenotype dataframe filtered for common samples, and reordered to RNA ordering
    """
    # Load the smallest cancer type dataset to gather the gene names
    gene_names_df = pd.read_parquet('./data/GDC-PANCAN.htseq_fpkm-uq_TCGA-CHOL.parquet')
    gene_names = gene_names_df.index.tolist()
    
    # Load the phenotype dataset to gather the cancer types
    phenotype_df = pd.read_parquet('./data/GDC-PANCAN.basic_phenotype_processed.parquet')
    cancer_types = phenotype_df['project_id'].unique()

    # Load the survival dataset
    survival_df = pd.read_parquet('./data/GDC-PANCAN.survival_processed.parquet')




    
    
    # # Read in both RNA matrices (split for data storage purposes)
    # df_1 = pd.read_parquet('./data/GDC-PANCAN.htseq_fpkm-uq_1.parquet')
    # df_2 = pd.read_parquet('./data/GDC-PANCAN.htseq_fpkm-uq_2.parquet')
    
    # # Concatenate the two RNA files into one
    # df = pd.concat([df_1, df_2], axis=1)
    
    # # Load in the ID/Gene Mapping file and create a merged dataframe to map the RNA IDs to the gene names
    # mapping = pd.read_parquet('./data/gencode.v22.annotation.gene.parquet')
    # mapping.rename(columns={'id': 'xena_sample'}, inplace=True) # rename column to merge on
    # merged_df = pd.merge(mapping, df, on='xena_sample', how='outer', indicator=True) # merge the dataframes to map gene names

    # # Ensure the expression dataframe is in the format: indexed on gene names column labels as sample ids
    # merged_trimmed_df = merged_df.copy()
    # merged_trimmed_df.drop(columns=['xena_sample', 'chrom', 'chromStart', 'chromEnd', 'strand', '_merge'], axis=1, inplace=True)
    # merged_trimmed_df.set_index('gene', inplace=True)
    
    # # Load in survival and phenotype dataframes
    # survival_df = pd.read_parquet('./data/GDC-PANCAN.survival.parquet')
    # phenotype_df = pd.read_parquet('./data/GDC-PANCAN.basic_phenotype.parquet')

    # # Create a list of all samples in each dataframe
    # samples_rna = list(merged_trimmed_df.columns)
    # samples_pheno = list(phenotype_df['sample'].values)
    # samples_survival = list(survival_df['sample'].values)
    
    # # Find all common samples in all three lists
    # common_samples = list(set(samples_rna) & set(samples_pheno) & set(samples_survival))
    
    # # Subset and reorder all three datasets by common_samples
    # merged_trimmed_filtered_df = merged_trimmed_df[common_samples] # Filter merged_df by columns in common_samples
    # phenotype_filtered_df = phenotype_df[phenotype_df['sample'].isin(common_samples)]
    # survival_filtered_df = survival_df[survival_df['sample'].isin(common_samples)]
    
    # # Reorder phenotype and survival dataframes to match the columns of the rna matrix
    # column_order = list(merged_trimmed_filtered_df.columns)
    # phenotype_filtered_ordered_df = phenotype_filtered_df.set_index('sample').loc[column_order].reset_index()
    # survival_filtered_ordered_df = survival_filtered_df.set_index('sample').loc[column_order].reset_index()
    
    # return merged_trimmed_filtered_df, survival_filtered_ordered_df, phenotype_filtered_ordered_df

    return gene_names, cancer_types, phenotype_df, survival_df



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
        return True
    else:
        return False


def create_rna_matrix(cancer_types_entered):
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

    return df


def calculate_ssgsea(df, phenotype_df):
    """
    Calculates ssGSEA scores using GSEAPY ssGSEA.

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
    signature_name, genes_entered, cancer_types_entered, cut_point_entered = get_form_values()

    # Create a dictionary of signature and gene names
    signature = {signature_name: genes_entered}
    
    # Determine the number of threads to run the calculations on
    n_threads=threading.active_count()-1

    # Get a list of samples for cancer types
    # cancer_type_samples = list(phenotype_df.query('project_id in @cancer_types_entered')['sample'])
    # Subset the original RNA count matrix to only have samples in selected cancer types
    # subset_counts = df.loc[: , cancer_type_samples]
    
    # Calculate the ssGSEA scores
    scores = gp.ssgsea(data=df, gene_sets=signature, outdir=None, 
               sample_norm_method='rank', threads=n_threads, min_size=3,
               verbose=True)
    return scores.res2d


def create_km_plot(ssgsea_scores, survival_df):
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
    # legend.set_title('NES')
    # Save the plot - using bbox_inches='tight' to avoid cutting off any content
    return km_plot_figure


def download_output(ssgsea_scores, km_plot_figure):
    """
    Downloads ssGSEA data to CSV file and KM plot to PNG on the users' local machine.

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
    
    # ssGSEA export into CSV format
    # Create the full path for the ssGSEA scores
    ssgsea_file_path = os.path.join(downloads_folder, f'ssgsea_scores_{today}.csv')
    # Output to CSV to the specified filepath
    ssgsea_scores.to_csv(ssgsea_file_path, index=False)

    # KM plot export
    # Create the full path for the KM plot
    km_file_path = os.path.join(downloads_folder, f'km_plot_{today}.png')
    # Save KM plot to the file path
    km_plot_figure.savefig(km_file_path, bbox_inches='tight')


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
    genes_entered : list
        1 or more genes selected by the user.
    cancer_types_entered : list
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
    st.title(":dna: TCGA SIGvival")
    st.write(f"*Erin Cameron &nbsp;&nbsp;&nbsp; | &nbsp;&nbsp;&nbsp; [GitHub](https://github.com/erincameron11/survival-analysis-tool)*")
    st.divider()
    
    # Create a field for informational text
    st.write("Enter signature name, gene names, cancer types, and cut-point to generate ssGSEA scores and visualize survival outcomes with a Kaplan-Meier plot based on TCGA RNA and phenotype survival data.")

    # Call the load data method
    # df, survival_df, phenotype_df = load_data()
    gene_names, cancer_types, phenotype_df, survival_df = load_data()
    
    # Locate all gene names in a list
    # gene_names = df.index.unique()
    
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
        # cancer_types = phenotype_df['project_id'].unique()
        # TODO: add a PAN-Cancer option for all cancers?
        cancer_types_entered = st.multiselect(
            "Cancer Type:",
            (cancer_types),
            placeholder="Select cancer type",
            key='cancer_types_entered',
        )
    
        # Dropdown for cut-point
        cut_point_entered = st.selectbox(
            "Cut-Point:",
            ("Median", "Tertile", "Tertile - Top & Bottom only", "Quartile", "Quartile - Top & Bottom only"),
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
                # Auto scroll to ssGSEA calculation info message
                auto_scroll()

    # Block the form from submitting on Enter press with text_input (built-in streamlit functionality)
    block_form_submit()

    # Apply CSS for custom styling
    custom_css()

    # If the submit button was pressed and submitted successfully
    if st.session_state.get('form_submitted', False):
        df = create_rna_matrix(cancer_types_entered)
        
        # Calculate ssGSEA
        ssgsea_info = st.info('Calculating ssGSEA scores...', icon="🔄")
        ssgsea_scores = calculate_ssgsea(df, phenotype_df)
        ssgsea_info.empty()
        
        # Create the kaplan meier results
        km_plot_figure = create_km_plot(ssgsea_scores, survival_df)

        # Scroll down once calculations complete
        auto_scroll()
        
        # Display the results inside a container
        with st.container(border=True):
            # Display results subheader
            st.subheader("Results")
            
            # Use state sessions to get form values
            signature_name, genes_entered, cancer_types_entered, cut_point_entered = get_form_values()
            
            # Display the entered form values for user
            genes_entered_str = ", ".join(genes_entered)
            cancer_types_entered_str = ", ".join(cancer_types_entered)
            st.write("**Signature Name**: ", signature_name, "  \n**Gene Names**: ", genes_entered_str, 
                     "  \n**Cancer Types**: ", cancer_types_entered_str, "  \n**Cut-point**: ", cut_point_entered)
            st.divider()

            # Create placeholders to hold the ssGSEA, KM plot and download button content
            # ssgsea_placeholder = st.empty()
            km_plot_placeholder = st.empty()
            download_results_placeholder = st.empty()
            # with ssgsea_placeholder:
                # Display ssGSEA output
                # st.dataframe(ssgsea_scores.head())
            with km_plot_placeholder:
                # Display the KM plot image created
                st.pyplot(fig=km_plot_figure, use_container_width=True)
            with download_results_placeholder:
                st.button(":arrow_down: Download Results", on_click=download_output, args=(ssgsea_scores, km_plot_figure,))
    


# ------------------------------------ RUN THE APP ------------------------------------
if __name__ == "__main__":
    main()