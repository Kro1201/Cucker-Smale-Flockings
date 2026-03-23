from numpy import *
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

# The phi function we used for our simulations
def phi1(s, a, b):
    return (a/(s**b))
def phi2(s, a, b):
    return (a/((1+s**2)**(b/2)))

# The time and the time-step
M=1000
T=10
h=T/M

# The Explicit Euler function to evolve the flock
def euler(x, y, vx, vy, a, b):
    X =[x]
    Y =[y]
    V1 =[vx]
    V2 =[vy]
    lamb = 1
    X1 = x.copy()
    X2 = y.copy()
    V1temp = vx.copy()
    V2temp = vy.copy()
    for t in np.linspace(0, T, M):
        S1 = np.zeros(N)
        S2 = np.zeros(N)
    for i in range(N):
        for j in range(N):
            vv1 = phi1(sqrt((X1[i]-X1[j])**2+(X2[i]-X2[j])**2), a, b)*(V1temp[j]-V1temp[i])
            vv2 = phi1(sqrt((X1[i]-X1[j])**2+(X2[i]-X2[j])**2), a, b)*(V2temp[j]-V2temp[i])
            S1[i] += vv1
            S2[i] += vv2
    V1temp_new = V1temp + (lamb*h*S1)/N
    V2temp_new = V2temp + (lamb*h*S2)/N
    for i in range(N):
        X1[i] += h*V1temp_new[i]
        X2[i] += h*V2temp_new[i]
    V1temp = V1temp_new
    V2temp = V2temp_new
    X.append(X1.copy())
    Y.append(X2.copy())
    V1.append(V1temp.copy())
    V2.append(V2temp.copy())
    return X, Y, V1, V2

# The function used to display the flock
def tracerfleche(x1, x2, x3, x4):
    plt.scatter(np.mean(x1), np.mean(x2), color='green')
    plt.quiver(np.mean(x1), np.mean(x2), np.mean(x3), np.mean(x4), color='green')
    plt.scatter(x1, x2, color='blue', label="Positions")
    #plt.xlim(np.mean(x1)-5, np.mean(x1)+5)
    #plt.ylim(np.mean(x2)-5, np.mean(x2)+5)
    plt.quiver(x1, x2, x3, x4, color='red', label="Velocities")
    plt.show()

# The initial conditions
N = 20 # Number of particles
# A position and a velocity in R2 are assigned to each particle
x = np.random.normal(0, 1, N)
y = np.random.normal(0, 1, N)
vx = np.random.normal(0, 1, N)
vy = np.random.normal(0, 1, N)

# The parameters of the phi function
a = 1
b = 1

# We use the Euler function to simulate the evolution of our random flock
X, Y, V1, V2 = euler(x, y, vx, vy, a, b)
tracerfleche(x, y, vc, vy)
tracerfleche(X[M], Y[M], V1[M], V2[M])

# This part plots the convergence of the velocities and positions
def diffmoy(Z1, Z2): 
    Y1 = []
    for t in range(len(Z1)):
        B1 = []
        Vm = [np.mean(Z1[t]), np.mean(Z2[t])]
        for j in range(N):
            u = sqrt((Z1[t][j]-Vm[0])**2+(Z2[t][j]-Vm[1])**2)
            B1.append(u)
            m = max(B1)
        Y1.append(m)
    return Y1

Y1 = diffmoy(V1, V2)
plt.semilogy(np.linspace(0, T, len(Y1)), Y1)
plt.xlabel("Time")
plt.ylabel("Max (vi-vm)")
plt.title("Speed-alignment")
plt.grid(True)
plt.show() 

xminorme = diffmoy(X, Y)
plt.plot(range(len(xminorme)), xminorme)

plt.xlabel("Time")
plt.ylabel("Max (xi-xm)")
plt.title("The positions are bounded")
plt.grid(True)
plt.show()

# Algorithm to find the critical value of b for which the flock converges

# Define the 2n order ODE system
def equation(t, x):
    dxdt = x[1]
    d2xdt2 = -dxdt / np.sqrt(1+2*x[0])
    return [dxdt, d2xdt2]

# Function to simulate the system for given initial conditions
def simulate_system(x0, v0, t_span = (0, 100)):
    sol = solve_ivp(equation, t_span, [x0, v0], t_eval=np.linspace(t_span[0], t_span[1], 1000))
    return sol

# Initial velocity
v0 = -1.0


