# ------------------------------------ IMPORTS ------------------------------------
import streamlit as st # for UI
import matplotlib.pyplot as plt # for KM plots
# Import external files
from data import *
from helpers import *
from styling import *


# ------------------------------------ MAIN APP ------------------------------------
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
    gene_names, cancer_types, phenotype_df, survival_df = load_data()

    # <%%%%%%%%%% TESTING
    gc.collect()
    memory_load = get_memory_usage()
    st.write(f"Memory usage after data load: {memory_load:.2f} MB")
    
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
        df = create_rna_dataframe(cancer_types_entered)

        # <%%%%%%%%%% TESTING
        gc.collect()
        memory_before = get_memory_usage()
        st.write(f"Memory usage before calculation: {memory_before:.2f} MB")
        
        # Calculate ssGSEA
        ssgsea_info = st.info('Calculating ssGSEA scores...', icon="🔄")
        ssgsea_scores = calculate_ssgsea(df, phenotype_df)
        ssgsea_info.empty()

        # <%%%%%%%%%% TESTING
        gc.collect()
        memory_after = get_memory_usage()
        st.write(f"Memory usage after calculation: {memory_after:.2f} MB")
        
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

            # Display the Kaplan Meier plot and download button
            st.pyplot(fig=km_plot_figure, use_container_width=True)
            st.button(":arrow_down: Download Results", on_click=download_output, args=(ssgsea_scores, km_plot_figure,))
    


# ------------------------------------ RUN THE APP ------------------------------------
if __name__ == "__main__":
    main()