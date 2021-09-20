import numpy as np

import m31 as m31
# print(m31.vx(), m31.vy(), m31.vz())


# rot: observed frame (x, y, z) ==>> M31 disk frame (X, Y, Z)
# inv: M31 disk frame (X, Y, Z) ==>> observed frame (x, y, z)
rot, inv = m31.set_rotation_matrix()

# orbit in Fardal et al. (2007)
pos = np.array([-34.75,  19.37, -13.99])# unit is kpc
vel = np.array([ 67.34, -26.12,  13.50])# unit is km/s

print("initial position in M31 disk frame:")
print(np.dot(rot, pos))

print("initial velocity in M31 disk frame:")
print(np.dot(rot, vel))


# rotation matrix
print("rotation matrix: observed frame (x, y, z) ==>> M31 disk frame (X, Y, Z)")
print(rot)

print("rotation matrix: M31 disk frame (X, Y, Z) ==>> observed frame (x, y, z)")
print(inv)


# orbital angular momentum
orm = np.cross(pos, vel)# set in observed frame
print ("orbital angular momentum in observed frame:")
print(orm)
print ("orbital angular momentum in M31 disk frame:")
print(np.dot(rot, orm))

print ("direction of orbital angular momentum in observed frame:")
print(orm / np.linalg.norm(orm))
print ("direction of orbital angular momentum in M31 disk frame:")
print(np.dot(rot, orm) / np.linalg.norm(np.dot(rot, orm)))


# spin angular momentum
spin = np.array([0.0, 0.0, 1.0])# set in M31 disk frame
print ("direction of spin angular momentum in observed frame:")
print(np.dot(inv, spin))
print ("direction of spin angular momentum in M31 disk frame:")
print(spin)


# colliding point of NW Stream and a DM subhalo
# xi, eta, D = -5.7083, 5.2083, 870.17
# vrot = 147.133755043495863996
# xi, eta, D = 1.1404, 8.2456, 908.22
# vrot = 137.6723563
# xi, eta, D = 1.41605839416058394161, 8.53291536050156739812, 913.013698630136986301
# vrot = 136.003236111249302739
xi, eta, D = -5.84952978056426332288, 8.0, 923.972602739726027397
vrot = 131.309100747738482496
x0, y0, z0 = m31.cartesian_coordinate(xi, eta, D)
obs = np.array([x0, y0, z0])
print("collision point in observed frame:")
print(obs, np.sqrt(x0**2 + y0**2 + z0**2))
print("collision point in M31 disk frame:")
print(np.dot(rot, obs))

vx0 = vrot / np.sqrt(1.0 + (y0 / x0)**2 + ((x0**2 + y0**2) / (z0 * x0))**2)
vy0 = vx0 * y0 / x0
vz0 = -vx0 * (x0**2 + y0**2) / (z0 * x0)
print("collision velocity 0 in observed frame:")
vel = np.array([vx0, vy0, vz0])
print(vel, np.sqrt(vx0**2 + vy0**2 + vz0**2))
print("collision velocity 0 in M31 disk frame:")
print(np.dot(rot, vel))

A = -(x0 * vz0 - z0 * vx0) / (y0 * vz0 - z0 * vy0)
B = -(x0 + A * y0) / z0
vx1 = vrot / np.sqrt(1.0 + A * A + B * B)
vy1 = A * vx1
vz1 = B * vx1
print("collision velocity 1 in observed frame:")
vel = np.array([vx1, vy1, vz1])
print(vel, np.sqrt(vx1**2 + vy1**2 + vz1**2))
print("collision velocity 1 in M31 disk frame:")
print(np.dot(rot, vel))


# inflowing angle of the GSS progenitor
disk = np.array([0.0, 0.0, 1.0])
view = np.dot(inv, disk)
print("rotation axis of the M31 in observed frame:")
print(view)
flow = np.array([394.0, -353.0, 316.0]) # in units of km/s
vvec = flow / np.linalg.norm(flow)
cosv = np.dot(view, vvec)
out = np.cross(view, vvec)
print("cosine:")
print(cosv)
print("sine:")
print(np.linalg.norm(out))
# print(cosv**2 + np.linalg.norm(out) **2)
print("arcsin:", np.arcsin(np.linalg.norm(out)), np.rad2deg(np.arcsin(np.linalg.norm(out))))
print("arccos:", np.arccos(cosv), np.rad2deg(np.arccos(cosv)))
