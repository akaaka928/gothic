import numpy as np

import m31 as m31

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
xi, eta, D = -5.7083, 5.2083, 870.17
x0, y0, z0 = m31.cartesian_coordinate(xi, eta, D)
obs = np.array([x0, y0, z0])
print("collision point in observed frame:")
print(obs, np.sqrt(x0**2 + y0**2 + z0**2))
print("collision point in M31 disk frame:")
print(np.dot(rot, obs))

vx0 = 147.133755043495863996 / np.sqrt(1.0 + (y0 / x0)**2 + ((x0**2 + y0**2) / (z0 * x0))**2)
vy0 = vx0 * y0 / x0
vz0 = -vx0 * (x0**2 + y0**2) / (z0 * x0)
print("collision velocity 0 in observed frame:")
vel = np.array([vx0, vy0, vz0])
print(vel, np.sqrt(vx0**2 + vy0**2 + vz0**2))
print("collision velocity 0 in M31 disk frame:")
print(np.dot(rot, vel))

A = -(x0 * vz0 - z0 * vx0) / (y0 * vz0 - z0 * vy0)
B = -(x0 + A * y0) / z0
vx1 = 147.133755043495863996 / np.sqrt(1.0 + A * A + B * B)
vy1 = A * vx1
vz1 = B * vx1
print("collision velocity 1 in observed frame:")
vel = np.array([vx1, vy1, vz1])
print(vel, np.sqrt(vx1**2 + vy1**2 + vz1**2))
print("collision velocity 1 in M31 disk frame:")
print(np.dot(rot, vel))
