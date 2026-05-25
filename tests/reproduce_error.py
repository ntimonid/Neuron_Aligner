import sys
import os
sys.path.append(os.getcwd())

from src import utils

data_repository = 'data_repository'
try:
    acr2id, ancestorsById, neuriteLengthDistribution, acr_to_morpho_id = utils.load_useful_variables(data_repository)
    print("load_useful_variables successful")
except Exception as e:
    print(f"Error: {e}")
    # Print CWD for debugging
    print(f"Current working directory: {os.getcwd()}")
    # List files in the data_repository
    print(f"Files in {data_repository}: {os.listdir(data_repository)}")
