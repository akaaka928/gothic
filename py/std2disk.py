import numpy as np

# orbit in Fardal et al. (2007)
x0, y0, z0 = -34.75, 19.37, -13.99# unit is kpc
vx0, vy0, vz0 = 67.34, -26.12, 13.50# unit is km/s

# set Euler angle
inc = (180.0 - 77.0) * np.pi / 180.0
pa  =          37.0  * np.pi / 180.0

sini, sinp = np.sin(inc), np.sin(pa)
cosi, cosp = np.cos(inc), np.cos(pa)

# reference is line 69-71 of dm:~/work/HA-PACS/130616survey/ml/src/m31/coordinate.c
rot = np.array([[sinp, cosp, 0.0], [-cosi * cosp, cosi * sinp, sini], [sini * cosp, -sini * sinp, cosi]])
print("rotation matrix")
print(rot)

# rot: observed frame (x, y, z) ==>> disk orthogonal frame (X, Y, Z)
# inv: disk orthogonal frame (X, Y, Z) ==>> observed frame (x, y, z)
inv = np.linalg.inv(rot)

pos = np.array([x0, y0, z0])
vel = np.array([vx0, vy0, vz0])

print("initial position:")
print(np.dot(rot, pos))

print("initial velocity:")
print(np.dot(rot, vel))


ini = np.array([0.0, 0.0, 1.0])# rotation axis of M31's disk in (X, Y, Z)-coordinate
fin = np.dot(inv, ini)# rotation axis of M31's disk in (x, y, z)-coordinate
print("rotation axis of M31's disk in observed frame:")
print(fin)

axis = np.cross(ini, fin)
sintheta = np.sqrt(np.dot(axis, axis))
costheta = np.dot(ini, fin)
axis = axis / np.linalg.norm(axis)
print("rotation axis of the Euler angle:")
print(axis)
print("sine and cosine of the rotation from disk frame to observed frame")
print(sintheta, costheta)
print("rotation angle (disk frame to observed frame):")
print(np.arcsin(sintheta), np.arccos(costheta))
print("sin(acos(costheta)) and cos(asin(sintheta)):")
print(np.sin(np.arccos(costheta)), np.cos(np.arcsin(sintheta)))
print("sin(asin(sintheta)) and cos(acos(costheta)):")
print(np.sin(np.arcsin(sintheta)), np.cos(np.arccos(costheta)))
