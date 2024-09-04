# Import statements
import streamlit as st
from streamlit_tags import st_tags
import pandas as pd
import kaplanmeier as km

# ------------------ DATA ------------------
# Cache the dataframe using st.cache_data decorator
@st.cache_data
# Function to read in the RNA & ID/Gene mapping tsv files
def load_data(filename, mapping_filename):
    # RNA matrix 
    df = pd.read_csv(filename, sep='\t')
    # ID/Gene Mapping file
    mapping = pd.read_csv(mapping_filename, sep='\t')
    mapping.rename(columns={'id': 'xena_sample'}, inplace=True) # rename column to merge on
    merged_df = pd.merge(df, mapping, on='xena_sample', how='left')
    return merged_df


# ------------------ HELPER FUNCTIONS ------------------
def handle_submit(genes_entered):
    # TODO: If there are genes entered, display the kaplan meier plot and submit
    if genes_entered:
        print("submitted")

# Function to inject CSS and change the gene multiselect tag colours to green once entered
def change_multiselect_colours(tag_colour: str, outline_colour: str) -> None:
    tag_css = f"""
    <style>
    .stMultiSelect div[data-baseweb="select"] span[data-baseweb="tag"] {{
        background-color: {tag_colour} !important;
        color: white !important;
    }}
    </style>
    """
    st.markdown(tag_css, unsafe_allow_html=True)

    outline_css = f"""
    <style>
    .stMultiSelect div[data-baseweb="select"]:focus-within > div {{
        border-color: {outline_colour} !important;
    }}
    </style>
    """
    st.markdown(outline_css, unsafe_allow_html=True)


def customize_submit_button(initial_bg_color: str, initial_outline_color: str, initial_text_color: str,
                            hover_bg_color: str, hover_outline_color: str,
                            active_bg_color: str, active_outline_color: str) -> None:
    css = f"""
    <style>
    /* Initial button style */
    div.stButton > button {{
        background-color: {initial_bg_color} !important;
        border-color: {initial_outline_color} !important;
        color: {initial_text_color} !important;
    }}
    
    /* Button hover effect */
    div.stButton > button:hover {{
        background-color: {hover_bg_color} !important;
        border-color: {hover_outline_color} !important;
        color: white !important;
    }}

    /* Button active effect */
    div.stButton > button:active {{
        background-color: {active_bg_color} !important;
        border-color: {active_outline_color} !important;
        color: white !important;
    }}
    </style>
    """
    st.markdown(css, unsafe_allow_html=True)


# ------------------ APP ------------------
def main():
    # App title
    st.title("SURVIVAL ANALYSIS")
    st.divider()
    
    # Create a field for informational text
    st.write("Informational text area")

    # Call the load data method
    df = load_data('./data/GDC-PANCAN.htseq_fpkm-uq.tsv', './data/gencode.v22.annotation.gene.probeMap')
    
    # Locate all gene names in a list
    gene_names = df['gene'].unique()
    
    # Multiselect element for gene selection
    genes_entered = st.multiselect(
        "Enter Genes:",
        gene_names,
    )

    # Apply the CSS to change the color of the multiselect tags etc.
    change_multiselect_colours("#4A9661", "gray")
    # Apply the CSS to change the color of the button and outline on hover and active
    customize_submit_button(
        initial_bg_color="white", 
        initial_outline_color="#D5D6D8", 
        initial_text_color="#31333F",
        hover_bg_color="#1f77b4", 
        hover_outline_color="#1f77b4",
        active_bg_color="#5a9bd4", 
        active_outline_color="#5a9bd4"
    )
    
    # Submit button
    st.button("Submit", on_click=handle_submit, args=(genes_entered,))


# ------------------ RUN THE APP ------------------
if __name__ == "__main__":
    main()