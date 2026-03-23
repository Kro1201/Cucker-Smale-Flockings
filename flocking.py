import numpy as np
from numpy import *
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

# The phi function we used for our simulations
def phi1(s, a, b):
    return (a/(s**b))
def phi2(s, a, b):
    return (a/((1+s**2)**(b/2)))


# Euler function to simulate the evolution of the flock
def euler_optimise(x, y, vx, vy, a, b, T, M):
    """
    Simulate the evolution of a flock using the Euler method.
    Return 
    """
    N = len(x)
    h = T/M
    lamb = 1.0

    # Initialize arrays to store positions and velocities
    X_hist = np.zeros((M+1, N))
    Y_hist = np.zeros((M+1, N))
    VX_hist = np.zeros((M+1, N))
    VY_hist = np.zeros((M+1, N))

    # Initial conditions
    X_hist[0] = x
    Y_hist[0] = y
    VX_hist[0] = vx
    VY_hist[0] = vy

    # Temporary arrays to store the new positions and velocities
    X1, X2 = x.copy(), y.copy()
    V1, V2 = vx.copy(), vy.copy()

    for t_step in range(1, M+1):
        # dX[i, j] corresponds to X_j - X_i
        dX = X1 - X1[:, np.newaxis]
        dY = X2 - X2[:, np.newaxis]
        dV1 = V1 - V1[:, np.newaxis]
        dV2 = V2 - V2[:, np.newaxis]

        # Compute the distances matrix
        dist = np.sqrt(dX**2 + dY**2)

        # Compute the weights matrix using phi2
        W = phi2(dist, a, b)

        # Weighted sum on the neighbouring particles axis (axis=1)
        S1 = np.sum(W * dV1, axis=1)
        S2 = np.sum(W * dV2, axis=1)

        # Updated velocities
        V1_new = V1 + (lamb * h * S1) / N
        V2_new = V2 + (lamb * h * S2) / N

        # Updated positions
        X1 = X1 + h * V1_new
        X2 = X2 + h * V2_new

        # Store the new positions and velocities
        V1, V2 = V1_new, V2_new
        X_hist[t_step] = X1
        Y_hist[t_step] = X2
        VX_hist[t_step] = V1
        VY_hist[t_step] = V2

    return X_hist, Y_hist, VX_hist, VY_hist

# Display functions

def tracerfleche(x, y, vx, vy, title="Flocks: Positions and Velocities"):
    """
    Display the particles' positions and velocities using a quiver plot.
    """
    plt.figure(figsize=(8, 6))
    plt.scatter(x, y, color='blue', label="Positions")
    plt.quiver(x, y, vx, vy, color='red', label="Velocities")

    # Display the mean velocity as a green arrow
    plt.scatter(np.mean(x), np.mean(y), color='green', marker='X', s=100, label="Gravity center")
    plt.quiver(np.mean(x), np.mean(y), np.mean(vx), np.mean(vy), color='green', width=0.005)

    plt.title(title)
    plt.legend()
    plt.grid(True)
    plt.show()

def calc_diff_moy_vectorise(V1_hist, V2_hist):
    """
    Calculate the maximum distance of each particle's velocity from the mean velocity at each time step.
    """
    Vm1 = np.mean(V1_hist, axis=1)
    Vm2 = np.mean(V2_hist, axis=1)

    # Calculate the distance of each particle's velocity from the mean velocity
    dist = np.sqrt((V1_hist - Vm1[:, np.newaxis])**2 + (V2_hist - Vm2[:, np.newaxis])**2)

    # Return the maximum distance for each time step
    return np.max(dist, axis=1)


