import numpy as np
from src.cpd_registration import RigidRegistration, NonRigidRegistration

def generate_synthetic_tree(n_points=50):
    # Generate a simple branch (a line in 3D)
    points = np.linspace(0, 10, n_points)
    return np.stack([points, np.zeros(n_points), np.zeros(n_points)], axis=1)

def deform_tree(points):
    # Apply a non-linear deformation (bending)
    deformed = points.copy()
    # Bend it along the Y-axis
    deformed[:, 1] = np.sin(points[:, 0] / 2)
    return deformed

def compute_mse(TY, X):
    # Mean squared error between registered points TY and target points X
    # Assuming corresponding points for simplicity (as in synthetic data)
    return np.mean(np.sum((TY - X) ** 2, axis=1))

# Generate data
n_points = 50
X = generate_synthetic_tree(n_points)
Y = deform_tree(X)

print(f"Generating {n_points} points.")
print(f"Target (X) shape: {X.shape}")
print(f"Source (Y) shape: {Y.shape}")

# Rigid Registration
print("\n--- Running Rigid Registration ---")
rigid_reg = RigidRegistration(X=X, Y=Y, max_iterations=100)
TY_rigid, params_rigid = rigid_reg.register()
mse_rigid = compute_mse(TY_rigid, X)
print(f"Rigid MSE: {mse_rigid:.4f}")

# Non-Rigid Registration
print("\n--- Running Non-Rigid Registration ---")
non_rigid_reg = NonRigidRegistration(X=X, Y=Y, max_iterations=100, beta=1.0, lambda_reg=0.1)
TY_non_rigid, params_non_rigid = non_rigid_reg.register()
mse_non_rigid = compute_mse(TY_non_rigid, X)
print(f"Non-Rigid MSE: {mse_non_rigid:.4f}")

print("\n--- Summary ---")
if mse_non_rigid < mse_rigid:
    print("Non-Rigid registration was more accurate.")
else:
    print("Rigid registration was more accurate.")
