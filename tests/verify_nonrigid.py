import numpy as np
from src.cpd_registration import NonRigidRegistration

# Simple 2D example
X = np.array([[0, 0], [1, 1], [2, 2]])
Y = np.array([[0, 0], [1, 1], [2, 2]])

# Introduce some noise to Y
Y[1, 0] += 0.1

print("Initializing registration...")
reg = NonRigidRegistration(X=X, Y=Y, max_iterations=10)
print("Running registration...")
TY, params = reg.register()

print("Registration finished.")
print("TY:\n", TY)
