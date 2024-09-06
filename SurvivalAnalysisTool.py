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
    merged_df = pd.merge(df, mapping, on='xena_sample', how='left')

    # Load in survival and phenotype dataframes
    survival_df = pd.read_parquet(survival_filename)
    phenotype_df = pd.read_parquet(phenotype_filename)

    # Create a list of all samples in each dataframe
    samples_rna = list(merged_df.columns)
    samples_pheno = list(phenotype_df['sample'].values)
    samples_survival = list(survival_df['sample'].values)
    
    # Find all common samples in all three lists
    common_samples = list(set(samples_rna) & set(samples_pheno) & set(samples_survival))
    
    # Subset and reorder all three datasets by common_samples
    merged_filtered_df = merged_df[common_samples] # Filter merged_df by columns in common_samples
    phenotype_filtered_df = phenotype_df[phenotype_df['sample'].isin(common_samples)]
    survival_filtered_df = survival_df[survival_df['sample'].isin(common_samples)]
    
    # Reorder phenotype and survival dataframes to match the columns of the rna matrix
    column_order = list(merged_filtered_df.columns)
    phenotype_filtered_ordered_df = phenotype_filtered_df.set_index('sample').loc[column_order].reset_index()
    survival_filtered_ordered_df = survival_filtered_df.set_index('sample').loc[column_order].reset_index()
    
    return merged_df, survival_filtered_ordered_df, phenotype_filtered_ordered_df


# ------------------------------------ HELPER FUNCTIONS ------------------------------------
def handle_submit():
    # If the entire form is filled out, calculate GSVA, display the kaplan meier plot and submit
    if validate_form():
        # Calculate GSVA
        calculate_gsva()
        # Display the kaplan meier plot
        create_km_plot()
        # Mark the form as submitted
        st.session_state.form_submitted = True
    else:
        with form_validation_placeholder:
            st.write("Please fill in all input fields")
        print("Form not complete - not submitted")

# Function to check if all data fields were filled out
def validate_form():
    # Use state sessions to get form values
    genes_entered = st.session_state.get('genes_entered', '')
    cancer_types_entered = st.session_state.get('cancer_types_entered', '')
    cut_point_entered = st.session_state.get('cut_point_entered', '')
    
    # If all form fields filled out, return True, else False
    if genes_entered and cancer_types_entered and cut_point_entered:
        return True
    else:
        return False

# TODO: Function to calculate GSVA scores
def calculate_gsva():
    # Use state sessions to get form values
    genes_entered = st.session_state.get('genes_entered', '')
    cancer_types_entered = st.session_state.get('cancer_types_entered', '')
    cut_point_entered = st.session_state.get('cut_point_entered', '')
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

def anchor():
    


# ------------------------------------ STYLING FUNCTIONS ------------------------------------
# Function to inject CSS and change the gene multiselect tag colours to green once entered
def change_multiselect_colours(tag_colour: str) -> None:
    # Multiselect Tag colour CSS styling
    tag_css = f"""
    <style>
    .stMultiSelect div[data-baseweb="select"] span[data-baseweb="tag"] {{
        background-color: {tag_colour} !important;
        color: white !important;
    }}
    </style>
    """
    st.markdown(tag_css, unsafe_allow_html=True)

    # Multiselect input outline CSS styling
    st.markdown("""
   <style>
    /* Initial border color */
    div[data-baseweb="select"] > div {
        border: 0em solid white !important;
        border-radius: 4px !important;
        box-shadow: 0 0 0 1px white !important; /* Maintain custom box-shadow */
    }
    /* Outline color when focused */
    div[data-baseweb="select"] > div:focus-within {
        border: 0em solid #c4c1c1;
        box-shadow: 0 0 0 1px #c4c1c1 !important;  
        border-color: #c4c1c1 !important;
    }
    </style>
    """, unsafe_allow_html=True)

# Function to customize the style of buttons
def customize_buttons(initial_bg_color: str, initial_outline_color: str, initial_text_color: str,
                            hover_bg_color: str, hover_outline_color: str,
                            active_bg_color: str, active_outline_color: str) -> None:
    css = f"""
    <style>
    /* Initial button style */
    div.stButton > button, div.stFormSubmitButton > button {{
        background-color: {initial_bg_color} !important;
        border-color: {initial_outline_color} !important;
        color: {initial_text_color} !important;
    }}
    
    /* Button hover effect */
    div.stButton > button:hover, div.stFormSubmitButton > button:hover {{
        background-color: {hover_bg_color} !important;
        border-color: {hover_outline_color} !important;
        color: white !important;
    }}

    /* Button active effect */
    div.stButton > button:active, div.stFormSubmitButton > button:active {{
        background-color: {active_bg_color} !important;
        border-color: {active_outline_color} !important;
        color: white !important;
    }}
    </style>
    """
    st.markdown(css, unsafe_allow_html=True)

    
# ------------------------------------ APP ------------------------------------
def main():    
    # App title
    st.title("SURVIVAL ANALYSIS")
    st.divider()
    
    # Create a field for informational text
    st.write("Describe KM plots, GSVA and what to do on this webpage...")

    # Apply the CSS to change the color of the multiselect tags etc.
    change_multiselect_colours("#4A9661")
    # Apply the CSS to change the color of the button and outline on hover and active
    customize_buttons(
        initial_bg_color="white", 
        initial_outline_color="#D5D6D8", 
        initial_text_color="#31333F",
        hover_bg_color="#1f77b4", 
        hover_outline_color="#1f77b4",
        active_bg_color="#5a9bd4", 
        active_outline_color="#5a9bd4"
    )

    # Call the load data method
    df, survival_df, phenotype_df = load_data('./data/GDC-PANCAN.htseq_fpkm-uq.parquet', 
                                              './data/gencode.v22.annotation.gene.probeMap',
                                              './data/GDC-PANCAN.survival.parquet', 
                                              './data/GDC-PANCAN.basic_phenotype.parquet')
    
    # Locate all gene names in a list
    gene_names = df['gene'].unique()

    # Create a form for data input
    with st.form("km_plot_form", clear_on_submit=False):
        # Create a placeholder for form validation error correction
        form_validation_placeholder = st.empty()
        
        # Multiselect element for gene selection
        genes_entered = st.multiselect(
            "Gene Names:",
            gene_names,
            placeholder="Enter gene names",
            key='genes_entered',
        )
        
        # Dropdown for cancer type
        cancer_types = phenotype_df['project_id'].unique()
        # TODO: what to do with "None" project_id's?
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

    with st.container():
        km_plot_anchor_placeholder = st.empty()
        gsva_placeholder = st.empty()
        km_plot_placeholder = st.empty()
        download_results_placeholder = st.empty()
    
        # Conditionally display the results and download button after form submission
        if st.session_state.get('form_submitted', False):
            with km_plot_anchor_placeholder:
                st.header('', anchor='km_plot_anchor')
                # Jump to the KM plot anchor using JavaScript
                js_code = """
                        <script>
                            document.getElementById('km_plot_anchor').scrollIntoView({behavior: 'smooth'});
                        </script>
                        """
                components.html(js_code)
            with gsva_placeholder:
                st.write("GSVA OUTPUT HERE")
            with km_plot_placeholder:
                st.write("KM PLOT OUTPUT HERE")
                fig = create_km_plot()
                st.pyplot(fig)
            with download_results_placeholder:
                st.button(":arrow_down: Download Results", on_click=download_output)


# ------------------------------------ RUN THE APP ------------------------------------
if __name__ == "__main__":
    main()