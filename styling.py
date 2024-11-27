# ------------------------------------ IMPORTS ------------------------------------
import streamlit as st # for UI
import streamlit.components.v1 as components # for page scrolling



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
        div.stButton > button, div.stFormSubmitButton > button, div.stDownloadButton > button {
            background-color: white !important;
            border-color: #D5D6D8 !important;
            color: #31333F !important;
        }
        /* Button hover effect */
        div.stButton > button:hover, div.stFormSubmitButton > button:hover, div.stDownloadButton > button:hover {
            background-color: #1f77b4 !important;
            border-color: #1f77b4 !important;
            color: white !important;
        }
        /* Button active effect */
        div.stButton > button:active, div.stFormSubmitButton > button:active, div.stDownloadButton > button:active {
            background-color: #5a9bd4 !important;
            border-color: #5a9bd4 !important;
            color: white !important;
        }
        </style>
    """, unsafe_allow_html=True)
    
    
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