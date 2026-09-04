import numpy as np
import requests
import json
import zlib
import sys
import os

# Add src to path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from src.neuron_morphology import NeuronMorphology
from src.synthetic_morphology import MorphologyInterpolator

def load_neuron_from_web(neuronName, db='mouselight'):
    neuronPath = 'https://neuroinformatics.nl/HBP/neuronsreunited-viewer/{}_json_gz/{}.json.gz'.format(db,neuronName)
    print(f"Fetching {neuronPath}...")
    response = requests.get(neuronPath)
    if response.status_code != 200:
        raise ValueError(f"Failed to fetch neuron {neuronName}. Status code: {response.status_code}")
    
    file_content = response.content
    # Decompress using zlib
    source_neuron_dict = json.loads(zlib.decompress(file_content, 16+zlib.MAX_WBITS))
    source_neuron_dict['fname'] = neuronName
    
    return NeuronMorphology(neuronDict=source_neuron_dict)

def compute_mse(pts1, pts2):
    # Simple MSE between corresponding points
    return np.mean(np.sum((pts1 - pts2)**2, axis=1))

def test_interpolation():
    print("Loading neurons...")
    # Using small neurons or just a few points if it's too slow
    # For testing purposes, we can also use subsampled versions
    neuron_id1 = 'AA0428'
    neuron_id2 = 'AA0571'
    
    morpho1 = load_neuron_from_web(neuron_id1)
    morpho2 = load_neuron_from_web(neuron_id2)
    
    # To speed up test, let's use a very small number of points
    # Actually, CPD with all points might be slow. 
    # But for a test, let's try it.
    
    print(f"Morpho 1 points: {len(morpho1.points)}")
    print(f"Morpho 2 points: {len(morpho2.points)}")
    
    interpolator = MorphologyInterpolator(morpho1, morpho2)
    
    weights = [0.0, 0.5, 1.0]
    results = {}
    
    for w in weights:
        print(f"Interpolating with weight {w}...")
        synth = interpolator.interpolate(w, max_iterations=20) # Low iterations for speed
        results[w] = synth
        print(f"Generated synthetic morphology for weight {w}")

    # Verification
    # Weight 0.0 should be exactly morpho2
    pts_0 = results[0.0].points[:, 0:3]
    pts_source = morpho2.points[:, 0:3]
    mse_0 = compute_mse(pts_0, pts_source)
    print(f"MSE(Weight=0, Source): {mse_0:.6f}")
    
    # Weight 1.0 should be the fully warped version
    pts_1 = results[1.0].points[:, 0:3]
    # We can't easily compare to morpho1 because they have different number of points
    # and CPD doesn't make them identical, just maps Y to X's space.
    
    # Check if 0.5 is indeed halfway between 0.0 and 1.0
    pts_05 = results[0.5].points[:, 0:3]
    expected_pts_05 = pts_source + 0.5 * (pts_1 - pts_source)
    mse_half = compute_mse(pts_05, expected_pts_05)
    print(f"MSE(Weight=0.5, LinearInterpolation): {mse_half:.6f}")
    
    assert mse_0 < 1e-6, "Weight 0.0 should match source exactly"
    assert mse_half < 1e-6, "Weight 0.5 should be linearly interpolated displacement"
    print("Interpolation test passed!")

if __name__ == "__main__":
    test_interpolation()
