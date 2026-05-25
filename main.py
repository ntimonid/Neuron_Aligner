import os
import threading
import time
import logging
from src import sba_api
from src.utils import *
import src.api_server
from src.api_server import run_server
from src import config

# Setup logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Global state variables
sba_interface = None
neuriteLengthDistribution = None
mySomaLocations = None
acr2id = None
ancestorsById = None

def initialize_application():
    global sba_interface, neuriteLengthDistribution, mySomaLocations, acr2id, ancestorsById
    
    data_repository = config.DATA_REPOSITORY
    myRegion = config.MY_REGION
    
    # Load useful variables
    logger.info("Loading useful variables...")
    acr2id, ancestorsById, neuriteLengthDistribution, acr_to_morpho_id = load_useful_variables(data_repository)
    logger.info("Useful variables loaded.")
    
    # Get the cell bodies
    mySomaLocations = get_soma_locations(myRegion, neuriteLengthDistribution, acr2id, ancestorsById)
    logger.info(f"Loaded {len(mySomaLocations.keys())} soma locations.")
    
    # Create SBA interface
    sba_interface = sba_api.SBA_interface()
    
    # Set parameters
    cpd_params = config.CPD_PARAMS
    soma_thr = config.SOMA_THR
    sba_interface.set_cpd_params(cpd_params, soma_thr)
    
    # Create objects corresponding to somata of neurons residing in region of interest (ex. VPM) to SBA
    sbaCommand_center = sba_interface.create_soma_object(mySomaLocations, myRegion)
    
    # Interface is ready
    logger.info("SBA Interface ready.")

if __name__ == "__main__":
    initialize_application()
    
    # Initialize API Server state
    src.api_server.set_api_state(sba_interface, neuriteLengthDistribution, mySomaLocations)

    # Start server in a separate thread
    server_thread = threading.Thread(target=run_server)
    server_thread.daemon = True
    server_thread.start()
    
    logger.info(f"Application running. API Server started on port {config.PORT}.")
    
    # Keep main thread alive
    while True:
        time.sleep(1)
