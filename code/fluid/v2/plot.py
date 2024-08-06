import numpy as np
import matplotlib.pyplot as plt

Nt = 50
Nx = 16
Ny = 16


plt.figure(figsize=(10, 10))
plt.ion()

def velocity(u,v):
    return([[np.sqrt(u[i][j]**2 + v[i][j]**2)for j in range(Ny)]for i in range(Nx)])

# Read the data
for time in range(Nt):
    plt.clf()
    filename = 'res/res_' + str(time) + '.dat'
    data = np.loadtxt(filename)

    # Reshape the data to 2D arrays
    u_t = data[:, 0].reshape((Nx, Ny))
    v_t = data[:, 1].reshape((Nx, Ny))

    # Create a grid for the plot
    x = np.linspace(0, Nx-1, Nx)
    y = np.linspace(0, Ny-1, Ny)
    X, Y = np.meshgrid(x, y)

    V = velocity(u_t, v_t)

    # Plot the vector field
    plt.streamplot(X, Y, u_t, v_t, density=0.8, linewidth=1, arrowsize=1.5, arrowstyle='->', color='white')
    #plt.contourf(X, Y, V, alpha=0.5, cmap='plasma')
    plt.imshow(V, alpha=0.5, cmap='plasma')
    plt.colorbar(label='Velocity')
    plt.title('Fluid Simulation Vector Field, t='+str(time))
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.scatter([int(3*Nx/7),int(3*Nx/7)],[int(3*Nx/7),int(4*Nx/7)])
    plt.scatter([int(3*Nx/7),int(4*Nx/7)],[int(4*Nx/7),int(4*Nx/7)])
    plt.scatter([int(4*Nx/7),int(4*Nx/7)],[int(4*Nx/7),int(3*Nx/7)])
    plt.scatter([int(4*Nx/7),int(3*Nx/7)],[int(3*Nx/7),int(3*Nx/7)])
    plt.savefig("res/fig_"+str(time)+".png")
    plt.pause(0.05)
