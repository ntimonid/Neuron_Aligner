import numpy as np
from src.cpd_registration import NonRigidRegistration

# Create a simple non-linear warping test
X = np.array([[0, 0], [1, 0], [2, 0]])
Y = np.array([[0, 0], [1, 1], [2, 0]]) # Bended upwards

# The registration should bring Y back to X (straight line)

print("Initializing registration...")
reg = NonRigidRegistration(X=X, Y=Y, max_iterations=100)
print("Running registration...")
TY, params = reg.register()

print("Registration finished.")
print("TY:\n", TY)
