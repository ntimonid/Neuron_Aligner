import os
import json

PORT = 8000
# Assuming config.py is in src/, project root is one level up
PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_REPOSITORY = os.path.join(PROJECT_ROOT, 'data_repository')

# Load configuration from JSON
CONFIG_PATH = os.path.join(PROJECT_ROOT, 'config.json')
with open(CONFIG_PATH, 'r') as f:
    config = json.load(f)

MY_REGION = config.get('myRegion', 'VPM')
CPD_PARAMS = config.get('cpd_params', {'max_it': 3, 'flag_in' : [1,1,-1], 'tol' : 0.001, 'branch_constraint': False})
SOMA_THR = config.get('soma_thr', 30)
