from vpython import *
import numpy as np

N = 100
R, lamda = 1.0, 500E-9
d = 100E-6

dx, dy = d/N, d/N
scene1 = canvas(align='left', height=600, width=600,
                center=vector(N*dx/2, N*dy/2, 0))
scene2 = canvas(align='right', x=600, height=600,
                width=600, center=vector(N*dx/2, N*dy/2, 0))
scene1.lights, scene2.lights = [], []
scene1.ambient, scene2.ambient = color.gray(0.99), color.gray(0.99)
side = np.linspace(-0.01*pi, 0.01*pi, N)
x, y = np.meshgrid(side, side)


# Calculate E_field
l_side = np.linspace(-N/2*(10**(-6)), N/2*(10**(-6)), N)
l_x, l_y = np.meshgrid(l_side, l_side)
E_field = np.zeros((N, N))
k = 2*pi/lamda
k_x = k*x/R
k_y = k*y/R
for i in range(N):
    j = 50
    for j in range(N):
        X = l_x[i][j]
        Y = l_y[i][j]
        if X**2 + Y**2 > (d/2)**2:
            continue
        E_field += 1/R * (dx*dy) * np.cos(k_x*X+k_y*Y)

# Get radius of first dark ring
# theta = -1
for i in range(50, 0, -1):
    if abs(E_field[i][50]) < abs(E_field[i+1][50]) \
          and abs(E_field[i][50]) < abs(E_field[i-1][50]):
        theta = sqrt(x[i][50]**2 + y[i][50]**2)/R
        rayleigh = 1.22*lamda/d
        print(f"Simulated theta             = {theta:.3E}")
        print(f"Rayleigh criterion          = {rayleigh:.3E}")
        print(f"Relative error to Rayleigh  = {abs(theta-rayleigh)/rayleigh:.2%}")
        break




Inte = abs(E_field) ** 2
maxI = np.amax(Inte)
for i in range(N):
    for j in range(N):
        box(canvas=scene1, pos=vector(i*dx, j*dy, 0), length=dx, height=dy, width=dx,
            color=vector(Inte[i, j]/maxI, Inte[i, j]/maxI, Inte[i, j]/maxI))

Inte = abs(E_field)
maxI = np.amax(Inte)
for i in range(N):
    for j in range(N):
        box(canvas=scene2, pos=vector(i*dx, j*dy, 0), length=dx, height=dy, width=dx,
            color=vector(Inte[i, j]/maxI, Inte[i, j]/maxI, Inte[i, j]/maxI))
