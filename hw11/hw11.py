from vpython import *
import math

scene = canvas(background=vec(0.8, 0.8, 0.8), width=1200,
               height=300, center=vec(3, 0, 10), fov=0.004)

lens_surface1 = shapes.arc(radius=0.15, angle1=0, angle2=pi)
circle1 = paths.arc(pos=vec(0, 0, 0), radius=0.0000001,
                    angle2=2*pi, up=vec(1, 0, 0))
lens_surface2 = shapes.arc(radius=0.15, angle1=-pi, angle2=0)
circle2 = paths.arc(pos=vec(0, 0, 0), radius=0.0000001,
                    angle2=2*pi, up=vec(1, 0, 0))
extrusion(path=circle1, shape=lens_surface1, color=color.yellow, opacity=0.6)
extrusion(path=circle2, shape=lens_surface2, color=color.yellow, opacity=0.6)
curve(pos=[vec(-7, 0, 0), vec(13, 0, 0)], color=color.red, radius=0.02)
arrow(pos=vec(-6, 0, 0), axis=vec(0, 0.5, 0), shaftwidth=0.1)
arrow(pos=vec(12, 0, 0), axis=vec(0, -1, 0), shaftwidth=0.1)


def refraction_vector(n1, n2, v_in, normal_v):
    # find the unit vector of velocity of the outgoing ray
    theta_1 = math.acos(dot(norm(v_in), norm(normal_v)))
    if v_in.y < 0:
        theta_1 = -theta_1
    theta_2 = math.asin(n1*math.sin(theta_1)/n2)
    delta = theta_2-theta_1
    basis_1 = vec(math.cos(delta), math.sin(delta), 0)
    basis_2 = vec(-math.sin(delta), math.cos(delta), 0)
    v_out = v_in.x * basis_1 + v_in.y * basis_2
    return v_out

R = 4.0
thickness = 0.3
g1center = vec(-R + thickness/2, 0, 0)
g2center = vec(R - thickness/2, 0, 0)
nair = 1
nglass = 1.5

for angle in range(-7, 2):
    ray = sphere(pos=vec(-6, 0.5, 0), color=color.blue,
             radius=0.01, make_trail=True)
    ray.v = vector(cos(angle/40.0), sin(angle/40.0), 0)
    dt = 0.002
    in_air = True
    while True:
        rate(1000)
        ray.pos = ray.pos + ray.v*dt
        # your code here
        if ray.pos.x < 0:
            if in_air and abs(mag(ray.pos-g2center)-R)/R <= 10**(-3.4):
                in_air = not in_air
                normal_v = norm(g2center-ray.pos)
                ray.v = refraction_vector(nair, nglass, ray.v, normal_v)
        else:
            if not in_air and abs(mag(ray.pos-g1center)-R)/R <= 10**(-3.4):
                in_air = not in_air
                normal_v = norm(ray.pos-g1center)
                ray.v = refraction_vector(nglass, nair, ray.v, normal_v)
            

        if ray.pos.x >= 12:
            print(ray.pos.y)
            break
