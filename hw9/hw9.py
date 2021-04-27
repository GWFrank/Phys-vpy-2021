import numpy as np
import time
st = time.time()

m = 2*10**3 # disk partition
n = 2*10**3 # loop partition

mu_0 = 1.256637 * 10**(-6)
R = 12 * 10**(-2)
r = 6 * 10**(-2)
H = 10 * 10**(-2)

# (1)
large_loop = np.empty((n, 3), dtype=np.float64)
for i in range(n):
    large_loop[i][0] = R*np.cos(i/n*2*np.pi)
    large_loop[i][1] = R*np.sin(i/n*2*np.pi)
    large_loop[i][2] = 0
current = 1

ds = np.empty((n, 3), dtype=np.float64)
for i in range(n):
    ds[i][0] = -(2*np.pi*R/n)*np.sin(i/n*2*np.pi)
    ds[i][1] = (2*np.pi*R/n)*np.cos(i/n*2*np.pi)
    ds[i][2] = 0

Phi_B = np.empty(m, dtype=np.float64)
for i in range(1, m+1):
    r_i = i/m * r
    p = np.array([r_i, 0, H], dtype=np.float64)
    r_vec = p-large_loop
    r_len = np.sqrt(np.sum(r_vec**2, axis=1))
    dB = mu_0/(4*np.pi) * current * np.cross(ds, r_vec) / ((r_len**3)[:, np.newaxis])
    B_ring = np.sum(dB, axis=0)[2]
    Phi_B[i-1] = B_ring * np.pi * (r_i**2 - ((i-1)/m*r)**2)

total_Phi_1 = np.sum(Phi_B)

# (2)
small_loop = np.empty((n, 3), dtype=np.float64)
for i in range(n):
    small_loop[i][0] = r*np.cos(i/n*2*np.pi)
    small_loop[i][1] = r*np.sin(i/n*2*np.pi)
    small_loop[i][2] = 0
current = 1

ds = np.empty((n, 3), dtype=np.float64)
for i in range(n):
    ds[i][0] = -(2*np.pi*r/n)*np.sin(i/n*2*np.pi)
    ds[i][1] = (2*np.pi*r/n)*np.cos(i/n*2*np.pi)
    ds[i][2] = 0

Phi_B = np.empty(m, dtype=np.float64)
for i in range(1, m+1):
    R_i = i/m * R
    p = np.array([R_i, 0, H], dtype=np.float64)
    r_vec = p-small_loop
    r_len = np.sqrt(np.sum(r_vec**2, axis=1))
    dB = mu_0/(4*np.pi) * current * np.cross(ds, r_vec) / ((r_len**3)[:, np.newaxis])
    B_ring = np.sum(dB, axis=0)[2]
    Phi_B[i-1] = B_ring * np.pi * (R_i**2 - ((i-1)/m*R)**2)

total_Phi_2 = np.sum(Phi_B)
ed = time.time()

print("="*35)
print(f"calculation time    : {ed-st:.3f} sec")
print(f"ring partition number {m=}")
print(f"loop partition number {n=}")
print(f"magnetic flux by (1): {total_Phi_1:.6E}")
print(f"magnetic flux by (2): {total_Phi_2:.6E}")
print("="*35)
