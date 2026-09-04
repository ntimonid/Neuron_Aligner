import numpy as np
import copy
from src.cpd_registration import NonRigidRegistration
from src.neuron_morphology import NeuronMorphology

class MorphologyInterpolator:
    """
    Class to create a synthetic morphology that is a weighted interpolation between two morphologies.
    """

    def __init__(self, target_morphology, source_morphology):
        """
        Initialize the interpolator with two NeuronMorphology objects.
        
        Parameters:
        target_morphology: NeuronMorphology (morphology 1)
        source_morphology: NeuronMorphology (morphology 2)
        """
        self.morpho1 = target_morphology
        self.morpho2 = source_morphology
        
        # Extract points for registration
        # We use all points for a faithful representation of the morphology
        self.X = self.morpho1.points[:, 0:3]
        self.Y = self.morpho2.points[:, 0:3]

    def interpolate(self, weight, **cpd_kwargs):
        """
        Create a synthetic morphology.
        
        Parameters:
        weight: float between 0 and 1. 
                0.0 = morphology 2 (source)
                1.0 = fully warped to morphology 1 (target)
                0.6 = 60% close to morphology 1, 40% close to morphology 2.
                
        **cpd_kwargs: arguments for NonRigidRegistration (max_iterations, beta, lambda_reg, etc.)
        
        Returns:
        synthetic_morphology: NeuronMorphology
        """
        if not (0 <= weight <= 1):
            raise ValueError("Weight must be between 0 and 1.")

        # 1. Run Non-Rigid Registration from Y (morpho2) to X (morpho1)
        # Default CPD parameters if not provided
        cpd_params = {
            'max_iterations': 100,
            'beta': 2.0,
            'lambda_reg': 2.0,
            'tolerance': 0.001
        }
        cpd_params.update(cpd_kwargs)
        
        reg = NonRigidRegistration(X=self.X, Y=self.Y, **cpd_params)
        reg.register()
        
        # 2. Get the displacement field (GW)
        # In NonRigidRegistration.transform_point_cloud: TY = Y + G*W
        # The full displacement is G*W.
        displacement = np.dot(reg.G, reg.W)
        
        # 3. Apply weighted displacement to source points
        # If weight=0, points remain Y. If weight=1, points become Y + displacement (fully warped).
        interpolated_points_xyz = self.Y + weight * displacement
        
        # 4. Create new NeuronMorphology object
        # Clone the source morphology structure
        synth_dict = self.morpho2.asDict()
        
        # Update points (keeping original columns like radius/type if present)
        # points array in NeuronMorphology has multiple columns. Usually [x, y, z, r, type, ...]
        new_points = np.array(synth_dict['treePoints']['data'])
        new_points[:, 0:3] = interpolated_points_xyz
        synth_dict['treePoints']['data'] = new_points.tolist()
        
        # Create the new morphology
        synthetic_morpho = NeuronMorphology(neuronDict=synth_dict)
        
        return synthetic_morpho
