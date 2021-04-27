from vpython import *
import numpy as np

prob = 0.005
N, L = 400, 7E-9/2.0
E = 1000000
q, m, size = 1.6E-19, 1E-6/6E23, 0.1E-9  # artificial charge particle
t, dt, vrms = 0, 1E-16, 10000.0
atoms, atoms_v = [], []

# initialization
scene = canvas(width=600, height=600, align='left',
               background=vector(0.2, 0.2, 0))
scenev = canvas(width=600, height=600, align='left',
                range=4E4, background=vector(0.2, 0.2, 0))
container = box(canvas=scene, length=2*L, height=2*L,
                width=2*L, opacity=0.2, color=color.yellow)

pos_array = -L + 2*L*np.random.rand(N, 3)
X, Y, Z = np.random.normal(0, vrms, N), np.random.normal(
    0, vrms, N), np.random.normal(0, vrms, N)
v_array = np.transpose([X, Y, Z])


def a_to_v(a):  # array to vector
    return vector(a[0], a[1], a[2])


for i in range(N):
    atom = sphere(canvas=scene, pos=a_to_v(
        pos_array[i]), radius=size, color=a_to_v(np.random.rand(3, 1)))
    atoms.append(atom)
    atoms_v.append(sphere(canvas=scenev, pos=a_to_v(
        v_array[i]), radius=vrms/30, color=a_to_v(np.random.rand(3, 1))))

# the average velocity and two axes in velocity space
vd_ball = sphere(canvas=scenev, pos=vec(0, 0, 0),
                 radius=vrms/15, color=color.red)
x_axis = curve(canvas=scenev, pos=[
               vector(-2*vrms, 0, 0), vector(2*vrms, 0, 0)], radius=vrms/100)
y_axis = curve(canvas=scenev, pos=[vector(
    0, -2*vrms, 0), vector(0, 2*vrms, 0)], radius=vrms/100)
vv = vector(0, 0, 0)    # for calculating the average velocity
total_c = 0             # the total number of collisions

rng = np.random.default_rng()
coli_cnt = 0

while True:
    t += dt
    rate(10000)

    v_array[:, 0] += q*E/m*dt
    pos_array += v_array*dt  # calculate new positions for all atoms
    outside = abs(pos_array) >= L
    pos_array[outside] = - pos_array[outside]

    # handle collision here

    v_sq = np.sum(v_array**2, axis=1)
    v_rms = np.sqrt(np.mean(v_sq))
    # exit()
    for i in range(N):
        if rng.random() <= prob:
            coli_cnt += 1
            
            hori_angle = rng.random()*2*np.pi
            elev_angle = rng.random()*2*np.pi - pi
            v_array[i][0] = v_rms*np.cos(elev_angle)*np.cos(hori_angle)
            v_array[i][1] = v_rms*np.cos(elev_angle)*np.sin(hori_angle)
            v_array[i][2] = v_rms*np.sin(elev_angle)
        else:
            continue


    vv += a_to_v(np.sum(v_array, axis=0)/N)

    if int(t/dt) % 2000 == 0:
        # tau = 0  # need to be modified
        tau = (t*N/coli_cnt)
        vd = mag(vv/(t/dt))
        print(tau, vd, q*E*tau/m)
    vd_ball.pos = vv/(t/dt)

    for i in range(N):
        atoms_v[i].pos, atoms[i].pos = a_to_v(v_array[i]), a_to_v(pos_array[i])
