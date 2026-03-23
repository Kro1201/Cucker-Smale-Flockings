import numpy as np
import matplotlib.pyplot as plt
from flocking import euler_optimise, tracerfleche, calc_diff_moy_vectorise

if __name__ == "__main__":
    # Parameters
    N = 20  # Number of particles
    T = 10.0  # Total time
    M = 1000  # Number of time steps
    a = 1.0  # Parameter for phi function
    b = 2  # Parameter for phi function

    # Initial conditions: random positions and velocities
    np.random.seed(42)
    x = np.random.normal(0, 1, N)
    y = np.random.normal(0, 1, N)
    vx = np.random.normal(0, 1, N)
    vy = np.random.normal(0, 1, N)

    # Display the initial positions and velocities
    tracerfleche(x, y, vx, vy, title="Initial Positions and Velocities")

    # Simulate the evolution of the flock
    X_hist, Y_hist, V1_hist, V2_hist = euler_optimise(x, y, vx, vy, a, b, T, M)

    # Display the final positions and velocities
    tracerfleche(X_hist[-1], Y_hist[-1], V1_hist[-1], V2_hist[-1], title="Final Positions and Velocities")

    # Evolution of the velocities alignment over time
    diff_vitesses = calc_diff_moy_vectorise(V1_hist, V2_hist)

    plt.figure(figsize=(8, 6))
    plt.semilogy(np.linspace(0, T, M+1), diff_vitesses)
    plt.xlabel("Time")
    plt.ylabel("Max(|v_i - v_moy|)")
    plt.title("Evolution of Velocity Alignment")
    plt.grid(True)
    plt.show()

    # Evolution of the positions over time
    diff_positions = calc_diff_moy_vectorise(X_hist, Y_hist)

    plt.figure(figsize=(8, 6))
    plt.plot(np.linspace(0, T, M+1), diff_positions)
    plt.xlabel("Time")
    plt.ylabel("Max(|x_i - x_moy|)")
    plt.title("Evolution of spatial dispersion")
    plt.grid(True)
    plt.show()