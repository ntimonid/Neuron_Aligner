import numpy as np
import requests
import json
import zlib
from src.cpd_registration import RigidRegistration, NonRigidRegistration
from src.neuron_morphology import NeuronMorphology

def load_neuron_from_web(neuronName, db='mouselight'):
    neuronPath = 'https://neuroinformatics.nl/HBP/neuronsreunited-viewer/{}_json_gz/{}.json.gz'.format(db,neuronName)
    print(f"Fetching {neuronPath}...")
    response = requests.get(neuronPath)
    if response.status_code != 200:
        raise ValueError(f"Failed to fetch neuron {neuronName}. Status code: {response.status_code}")
    
    file_content = response.content
    # Decompress using zlib as in sba_api.py
    source_neuron_dict = json.loads(zlib.decompress(file_content, 16+zlib.MAX_WBITS))
    source_neuron_dict['fname'] = neuronName
    
    return NeuronMorphology(neuronDict=source_neuron_dict)

def get_points(morphology):
    # Based on usage in sba_api.py
    lines, points = morphology.subsampledCopy([1, 2], minDistance=10.0)
    # Convert list of points to numpy array, taking only first 3 coordinates (x,y,z)
    points_array = np.array(points)
    # The list contains [0,0,0,0,0] as the first element which should be ignored as per comments in subsampledCopy
    return points_array[1:, 0:3]

def compute_mse(TY, X):
    # Mean squared error using nearest neighbor distances
    # Distance between each point in TY and all points in X
    dists = np.sqrt(np.sum((TY[:, None, :] - X[None, :, :]) ** 2, axis=2))
    # Min distance to any point in X for each point in TY
    min_dists = np.min(dists, axis=1)
    return np.mean(min_dists ** 2)

# Load two neurons
# IDs from data_repository\neuriteLengthDistribution(mouselight).json
neuron_id1 = 'AA0428'
neuron_id2 = 'AA0571'

morpho1 = load_neuron_from_web(neuron_id1, db='mouselight')
morpho2 = load_neuron_from_web(neuron_id2, db='mouselight')

X = get_points(morpho1)
Y = get_points(morpho2)

print(f"Target (X) shape: {X.shape}")
print(f"Source (Y) shape: {Y.shape}")

# Rigid Registration
print("\n--- Running Rigid Registration ---")
rigid_reg = RigidRegistration(X=X, Y=Y, max_iterations=100)
TY_rigid, params_rigid = rigid_reg.register()
# Compute MSE
mse_rigid = compute_mse(TY_rigid, X)
print(f"Rigid MSE: {mse_rigid:.4f}")

# Non-Rigid Registration
print("\n--- Running Non-Rigid Registration ---")
non_rigid_reg = NonRigidRegistration(X=X, Y=Y, max_iterations=100, beta=1.0, lambda_reg=0.1)
TY_non_rigid, params_non_rigid = non_rigid_reg.register()
# Compute MSE
mse_non_rigid = compute_mse(TY_non_rigid, X)
print(f"Non-Rigid MSE: {mse_non_rigid:.4f}")

print("\n--- Summary ---")
if mse_non_rigid < mse_rigid:
    print("Non-Rigid registration was more accurate.")
else:
    print("Rigid registration was more accurate.")
