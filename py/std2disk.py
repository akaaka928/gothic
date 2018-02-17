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

print("rotation matrix: observed frame (x, y, z) ==>> M31 disk frame (X, Y, Z)")
print(rot)

print("rotation matrix: M31 disk frame (X, Y, Z) ==>> observed frame (x, y, z)")
print(inv)
