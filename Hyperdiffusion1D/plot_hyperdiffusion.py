
import numpy as np
import matplotlib.pyplot as plt

Nx = 1000
L = 1
T = 0.000025
dt = 0.000001

dx = L/Nx
x_vec = np.linspace(0, L, Nx+1, endpoint = True)
Nt = int(np.floor(T/dt))
t_vec = np.linspace(0, T, Nt+1, endpoint = True)

# mesh grid
X, T = np.meshgrid(x_vec, t_vec)

plot_hyperdiffusion = np.loadtxt('plot_hyperdiffusion.txt', delimiter=',')

plt.figure()
CS = plt.contourf(X, T, plot_hyperdiffusion, 50)
plt.xlabel('X')
plt.ylabel('T')
plt.title('hyperdiffusion')
plt.colorbar(CS, shrink=0.8, extend='both')
plt.show()
