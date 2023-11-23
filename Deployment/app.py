from urllib.request import BaseHandler
from wsgiref.handlers import BaseCGIHandler
import streamlit as st 
import pandas as pd
import subprocess
import os
import base64
#import pickle
import sys
import joblib

padel_descriptor_path = r'C:\Users\DELL\Desktop\CDD\CDD_RUNX1_Bioactivity_Data\Deployment\bioactivity-prediction-app-main'

# Molecular descriptor calculator
def desc_calc(padel_descriptor_path):
    # Perform the descriptor calculation
    bashCommand = f"java -Xms2G -Xmx2G -Djava.awt.headless=true -jar {padel_descriptor_path}\\PaDEL-Descriptor.jar -removesalt -standardizenitro -fingerprints -descriptortypes {padel_descriptor_path}\\PubchemFingerprinter.xml -dir {padel_descriptor_path} -file descriptors_output.csv"
    process = subprocess.Popen(bashCommand, stdout=subprocess.PIPE, shell=True)
    output, error = process.communicate()
    if os.path.exists('molecule.smi'):
        os.remove('molecule.smi')

# File download
def filedownload(df):
    csv = df.to_csv(index=False)
    b64 = base64.b64encode(csv.encode()).decode()  # strings <-> bytes conversions
    href = f'<a href="data:file/csv;base64, {b64}" download= "prediction.csv">Download Predictions</a>'
    return href
    
# Model building
def build_model(load_data, desc_subset):
    # Reads in saved regression model
    load_model = joblib.load(r'C:\Users\DELL\Desktop\CDD\CDD_RUNX1_Bioactivity_Data\Deployment\RUNX1_model.joblib')

    # Apply model to make predictions
    prediction = load_model.predict(desc_subset)
    st.header('**Prediction output**')
    prediction_output = pd.Series(prediction, name='pIC50')
    molecule_name = pd.Series(load_data[1], name='molecule_name')
    df = pd.concat([molecule_name, prediction_output], axis=1)
    st.write(df)

# Logo image
# image = Image.open('logo.png')
# st.image(image, use_column_width=True)

# Page title
st.markdown("""
            
<style>
.stApp {
    background-image: url('https://assets-global.website-files.com/60f3bab37e9f03de6ca10f7f/62c5e62dd35f070a766d3cae_AdobeStock_173905768.jpeg');  # Replace with your background image URL
    background-size: cover;
    background-repeat: no-repeat;
}
</style>            
# Bioactivity Prediction App (RUNX1)
### This app allows you to predict the bioactivity towards inhibiting the 'RUNX1' transcription factor. 'RUNX1' is a target protein for Acute Myeloid Leukemia.
""",
unsafe_allow_html=True
)

with st.sidebar.header('1. Upload your CSV data'):
    uploaded_file = st.sidebar.file_uploader('Upload your input file', type=['csv'])

if st.sidebar.button('Predict'):
    if uploaded_file is not None:
        try:
            load_data = pd.read_csv(uploaded_file, header=None, sep='\t' if uploaded_file.name.endswith('.tsv') else ',', encoding='utf-8')
        except UnicodeDecodeError:
            st.warning("UTF-8 decoding failed. Trying 'latin-1' encoding.")
            try:
                load_data = pd.read_csv(uploaded_file, header=None, sep='\t' if uploaded_file.name.endswith('.tsv') else ',', encoding='latin-1')
            except UnicodeDecodeError:
                st.error("Both UTF-8 and 'latin-1' decoding failed. Please check the file encoding.")
                sys.exit()
            except pd.errors.EmptyDataError:
                st.error("The file is empty. Please upload a file with valid data.")
                sys.exit()
        except pd.errors.EmptyDataError:
            st.error("The file is empty. Please upload a file with valid data.")
            sys.exit()

        # Rest of your code...
        load_data.to_csv('molecule.smi', sep='\t', header=False, index=False)
        st.header('**Original input data**')
        st.write(load_data)

        with st.spinner("Calculating descriptors..."):
            desc_calc(padel_descriptor_path)

        # Read in calculated descriptors and display the dataframe
        st.header('**Calculated molecular descriptors**')
        desc = pd.read_csv(r'C:\Users\DELL\Desktop\CDD\CDD_RUNX1_Bioactivity_Data\Data collected in each section\Part-3, Descriptor_Dataset_Preparation\descriptors_output.csv')
        st.write(desc)
        st.write(desc.shape)

        # Read descriptor list used in previously built model
        st.header('**Subset of descriptors from previously built model**')
        Xlist = list(pd.read_csv(r'C:\Users\DELL\Desktop\CDD\CDD_RUNX1_Bioactivity_Data\Data collected in each section\Part-4, Model Building\descriptor_list.csv').columns)
        desc_subset = desc[Xlist]
        st.write(desc_subset)
        st.write(desc_subset.shape)

        # Apply trained model to make predictions on query compounds
        build_model(load_data, desc_subset)
    else:
        st.warning('Please upload a file before attempting to predict.')
