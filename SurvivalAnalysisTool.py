# Import statements
import streamlit as st
from streamlit_tags import st_tags
import pandas as pd
import kaplanmeier as km
import matplotlib.pyplot as plt # TESTING for image export of KM plots
from GSVA import gsva
from datetime import datetime # TESTING for filename saving
import numpy as np # TESTING for st.pyplot
import streamlit.components.v1 as components # TESTING for km plot page anchor
import time # TESTING page anchor scrolling
import gseapy as gp # for GSVA calculation
import threading # for GSVA calculation

# ------------------------------------ DATA ------------------------------------
# Cache the dataframe using st.cache_data decorator
@st.cache_data
# Function to read in the RNA & ID/Gene mapping tsv files
def load_data(rna_filename, gene_mapping_filename, survival_filename, phenotype_filename):
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
def handle_submit(form_validation_placeholder):
    # If the entire form is filled out, calculate GSVA, display the kaplan meier plot and submit
    if validate_form():
        # Calculate GSVA
        calculate_gsva()
        # Display the kaplan meier plot
        create_km_plot()
        # Mark the form as submitted
        st.session_state.form_submitted = True

# Function to check if all data fields were filled out
def validate_form():
    # Use state sessions to get form values
    genes_entered = st.session_state.get('genes_entered', '')
    cancer_types_entered = st.session_state.get('cancer_types_entered', '')
    cut_point_entered = st.session_state.get('cut_point_entered', '')
    
    # If all form fields filled out, return True, else False
    if genes_entered and cancer_types_entered and cut_point_entered:
        valid_form = True
        return True
    else:
        valid_form = False
        return False

# TODO: Function to calculate GSVA scores
def calculate_gsva():
    # Use state sessions to get form values
    genes_entered = st.session_state.get('genes_entered', '')
    cancer_types_entered = st.session_state.get('cancer_types_entered', '')
    cut_point_entered = st.session_state.get('cut_point_entered', '')
    
    # # Determine the number of threads to run the calculations on
    # n_threads=threading.active_count()-1
    
    # # Calculate the GSVA scores
    # ss = gp.ssgsea(data=merged_trimmed_df, gene_sets=signature, outdir=None, 
    #            sample_norm_method='rank', threads=n_threads, min_size=2, 
    #            verbose=True)
    
    print("calculating the gsva")

# TODO: Function to display the Kaplan Meier plot
def create_km_plot():
    # TESTING -- create and show KM plot
    arr = np.random.normal(1, 1, size=100)
    fig, ax = plt.subplots()
    ax.hist(arr, bins=20)
    return fig

# TODO: Function to download the GSVA data to CSV file, and KM plot to PNG
def download_output():
    # Get the current date and time for file naming
    today = datetime.now().strftime("%Y-%m-%d_%H:%M:%S")
    
    # GSVA export
    # Transform GSVA data into CSV format
    gsva = {'col1': [1, 2], 'col2': [3, 4]}
    gsva_df = pd.DataFrame(data=gsva)
    gsva_df_csv = gsva_df.to_csv(f'gsva_{today}.csv', index=False)

    # TESTING -- KM plot export
    plt.plot([1, 2, 3], [1, 4, 9])
    plt.savefig(f'km_plot_{today}.png')

# Function to block the form from submitting on Enter press with text_input (built-in streamlit functionality)
def block_form_submit():
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


# ------------------------------------ STYLING FUNCTIONS ------------------------------------
# Function to inject CSS and change the gene multiselect tag colours to green once entered
def customize_multiselect_colours() -> None:
    # Multiselect input outline CSS styling
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
        </style>
    """, unsafe_allow_html=True)

# Function to inject CSS and change the gene signature text input colours
def customize_text_input() -> None:
    # Text Input colour CSS styling    
    st.markdown("""
        <style>
        /* Change the outline styling for text_inputs */
        .stTextInput > div[class]:focus-within {
            border-bottom-color: #c4c1c1 !important;
            border-top-color: #c4c1c1 !important;
            border-left-color: #c4c1c1 !important;
            border-right-color: #c4c1c1 !important;
        }
        [data-testid="InputInstructions"] { 
            display: None;
        }
        /* Change the placeholder text opacity */
        .stTextInput > div > div > input::placeholder {
            opacity: 1.0;
        }
        /* Change the input field border radius */
        .stTextInput > div {
            border-radius: 4px;
        }
        </style>
    """, unsafe_allow_html=True)

# Function to customize the style of buttons
def customize_buttons() -> None:
    st.markdown("""
        <style>
        /* Initial button style */
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
    # Define JavaScript code for auto-scroll to the results section
    scroll_js = '''
        <script>
            console.log(window.parent.document.querySelector(".main"));
            window.parent.document.querySelector(".main").scrollTo({top: 500, behavior: 'smooth'});
        </script>
        '''

    # App title
    st.title("SURVIVAL ANALYSIS")
    st.divider()
    
    # Create a field for informational text
    # st.write("Describe KM plots, GSVA and what to do on this webpage...")

    # Call the load data method
    df, survival_df, phenotype_df = load_data('./data/GDC-PANCAN.htseq_fpkm-uq.parquet', 
                                              './data/gencode.v22.annotation.gene.probeMap',
                                              './data/GDC-PANCAN.survival.parquet', 
                                              './data/GDC-PANCAN.basic_phenotype.parquet')

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
        submit_button = st.form_submit_button(":chart_with_downwards_trend: Create KM Plot", on_click=handle_submit, args=(form_validation_placeholder,))
        
        # If the submit button was clicked, check that all fields are filled in
        if submit_button:
            # If not all fields are filled in
            if not signature_name or not genes_entered or not cancer_types_entered or not cut_point_entered:
                # Display red text form validation
                with form_validation_placeholder:
                    st.markdown("""
                        <p style='color:#cc0000; text-align:center;'>Please fill out all fields</p>
                    """, unsafe_allow_html=True)
    
    # TODO: IF GSVA TAKES TOO LONG TO PROCESS, CONSIDER USING st.status WHILE LOADING RESULTS

    # Block the form from submitting on Enter press with text_input (built-in streamlit functionality)
    block_form_submit()

    # Apply CSS for custom styling
    customize_multiselect_colours()
    customize_buttons()
    customize_text_input()

    # If the submit button was pressed and submitted successfully
    if st.session_state.get('form_submitted', False):
        # Define an empty element for scrolling to results
        scroll_temp = st.empty()
        
        # Display the results
        with st.container(border=True):
            st.subheader("Results")

            # Use state sessions to get form values
            signature_name = st.session_state.get('signature_name', '')
            genes_entered = st.session_state.get('genes_entered', '')
            cancer_types_entered = st.session_state.get('cancer_types_entered', '')
            cut_point_entered = st.session_state.get('cut_point_entered', '')
            # Display the entered form values
            genes_entered_str = ", ".join(genes_entered)
            cancer_types_entered_str = ", ".join(cancer_types_entered)
            st.write("**Signature Name**: ", signature_name, "  \n**Gene Names**: ", genes_entered_str, "  \n**Cancer Types**: ", cancer_types_entered_str, "  \n**Cut-point**: ", cut_point_entered)
            st.divider()

            # Create placeholders to hold the gsva and km plot content
            gsva_placeholder = st.empty()
            km_plot_placeholder = st.empty()
            download_results_placeholder = st.empty()
            with gsva_placeholder:
                # TODO: Display GSVA output
                st.write("GSVA OUTPUT HERE")
            with km_plot_placeholder:
                # TODO: Display KM plot
                fig = create_km_plot()
                st.pyplot(fig)
            with download_results_placeholder:
                st.button(":arrow_down: Download Results", on_click=download_output)
            # Add in auto-scroll behaviour
            with scroll_temp:
                st.components.v1.html(scroll_js, height=0)
                time.sleep(.5) # Ensure the script can execute before being deleted
            scroll_temp.empty()

    
# ------------------------------------ RUN THE APP ------------------------------------
if __name__ == "__main__":
    main()