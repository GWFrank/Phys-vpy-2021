from vpython import*
import math
import time
from vpython.rate_control import INTERACT_PERIOD

from vpython.vpython import print_to_string

fd = 120 # 120Hz
R = 30
L = 200
C = 20
L_val = L*10**(-3)
C_val = C*10**(-6)
T = 1/fd

t = 0
dt = 1.0/(fd * 5000) # 5000 simulation points per cycle

scene1 = graph(align = 'left', xtitle='t', ytitle='i (A) blue, v (100V) red,', background=vector(0.2, 0.6, 0.2))
scene2 = graph(align = 'left', xtitle='t', ytitle='Energy (J)', background=vector(0.2, 0.6, 0.2))

i_t = gcurve(color=color.blue, graph = scene1)
v_t = gcurve(color=color.red, graph = scene1)
E_t = gcurve(color=color.red, graph = scene2)


def voltage(t):
    # return 36
    if t < 0 or t >= 12*T:
        return 0
    else:
        return 36*math.sin(2*math.pi*fd*t)

pre_pre_q = 0
pre_q = 0
pre_pre_I = 0
pre_I = 0

I_m = 48763
t_m = 48763

E_off = 48763
t_tenth = 48763

while t <= 20*T:
    t += dt
    V = voltage(t)
    
    coef1 = (R*C_val) / (dt)
    coef2 = (L_val*C_val) / (dt**2)
    q = ( C_val*V + coef1*(pre_q) + coef2*(2*pre_q-pre_pre_q) ) / (1 + coef1 + coef2)
    I = (q - pre_q)/dt
    I_prime = (I - pre_I)/dt
    E_C = q**2 / (2*C_val)
    E_L = (L_val * I**2)/2

    v_t.plot(pos = (t, V/100))
    i_t.plot(pos = (t, I))
    E_t.plot(pos = (t, E_C+E_L))
    
    if t >= 9*T and t <= 10*T:
        if pre_I == max(I, pre_I, pre_pre_I):
            I_m = pre_I
            t_m = t
    
    if abs((t-12*T) / T) <= 10**(-4):
        E_off = E_C+E_L
    
    if t > 12*T:
        if abs(((E_L+E_C)-0.1*E_off) / E_off) <= 10**(-5):
            t_tenth = t

    pre_pre_q = pre_q
    pre_q = q
    pre_pre_I = pre_I
    pre_I = I

phi = (2*math.pi)*(t_m/T) - (9.25*2*math.pi)
X_L = 2*math.pi*fd*L_val
X_C = 1 / (2*math.pi*fd*C_val)
I_theory = 36/sqrt(R**2 + (X_L-X_C)**2)
phi_theory = math.atan((X_L-X_C)/R)

print("="*20)
print("Simulated:")
print(f"I_m = {I_m:.5f}")
print(f"phi = {phi:.5f}")
print("="*20)
print("Theory:")
print(f"I_m = {I_theory:.5f}")
print(f"phi = {phi_theory:.5f}")
print("="*20)
print(f"10% energy t: {t_tenth}")
print("="*20)
