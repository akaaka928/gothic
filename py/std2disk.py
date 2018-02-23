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
