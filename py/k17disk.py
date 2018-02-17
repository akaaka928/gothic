# set rotation axis of the infalling satellite in the M31's disk coordinate
import numpy as np

phi   = -15.0
theta =  30.0

# implementation A
rotY = -theta * np.pi / 180.0
rotZ = (phi - 180.0) * np.pi / 180.0

sinY, cosY = np.sin(rotY), np.cos(rotY)
sinZ, cosZ = np.sin(rotZ), np.cos(rotZ)

rot = np.array([
    [cosZ * cosY, -sinZ, cosZ * sinY],
    [sinZ * cosY,  cosZ, sinZ * sinY],
    [      -sinY,   0.0,        cosY]
])
print("rotation matrix:")
print(rot)

spin = np.dot(rot, np.array([0.0, 0.0, 1.0]))
print("rotation axis of the satellite in the M31 disk frame:")
print(spin)


# implementation B
import m31 as m31
axis = np.array([1.0, 0.0, 0.0])
angle = -0.5 * np.pi
rot1, inv1 = m31.set_rodrigues_matrix(axis, np.sin(angle), np.cos(angle))

axis = np.dot(rot1, np.array([0.0, 0.0, 1.0]))
angle = (180.0 - theta) * np.pi / 180.0
rot2, inv2 = m31.set_rodrigues_matrix(axis, np.sin(angle), np.cos(angle))

axis = np.dot(rot1, np.array([0.0, 1.0, 0.0]))
angle = (180.0 - phi) * np.pi / 180.0
rot3, inv3 = m31.set_rodrigues_matrix(axis, np.sin(angle), np.cos(angle))

rot = np.dot(rot3, np.dot(rot2, rot1))
print("rotation matrix:")
print(rot)

spin = np.dot(rot, np.array([0.0, 0.0, 1.0]))
print("rotation axis of the satellite in the M31 disk frame:")
print(spin)
