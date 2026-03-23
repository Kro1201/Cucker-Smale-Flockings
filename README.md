# Flocking Models: Dynamical Systems and Cucker-Smale Model

This repository contains the research report and the numerical simulation code developed for my Master 1 thesis at Paris-Dauphine PSL University.

The project aims to understand the collective behavior (or flocking) of interacting agents, focusing primarily on the Cucker-Smale mathematical model.

📂 Repository Structure

    Flockings_Mémoire_M1_révisé.pdf: The thesis paper detailing the mathematical proofs, the reduction of the dynamical system, and the analysis of the simulations.

    flocking.py: The optimised Python script simulating the model. It uses a vectorized Euler method via NumPy for high-performance computation of positions and velocities over time.
    

📝 Project Description

The report analytically and numerically explores the necessary conditions for a group of agents to reach a velocity consensus (alignment) while remaining grouped (spatial confinement).

Key topics covered include:

    The reduction of the Cucker-Smale dynamical system.

    The fundamental Flocking Theorem.

    Analysis of singular and regular communication rates.

In our simulations, the regular interaction function between particles is parameterized as follows:
ψ(s)=(1+s2)β/2α​

The codebase allows users to test the emergence of unconditional or conditional flocking by tweaking the α and β parameters.

🛠️ Prerequisites and Installation

To run the simulation code, you need Python 3 installed along with standard scientific libraries.

```
pip install numpy matplotlib
```

🚀 Usage

Run the main script to execute the optimized simulation:
```
python main.py
```

Configuration:
You can directly modify the variables within main.py to explore different system behaviors:

    N: Number of particles.

    a (α) and b (β): Influence parameters dictating the strength of the interaction relative to distance.

    np.random.seed(): Fixes the initial state to ensure reproducible results.

At the end of the execution, the script generates three plots:

    A visual representation of the particles' positions and velocity vectors at the final state.

    The evolution of velocity alignment over time on a logarithmic scale (convergence toward 0).

    The evolution of spatial dispersion (positions must remain bounded to confirm flocking).
