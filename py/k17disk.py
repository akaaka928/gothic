import numpy as np

phi = -15.0
theta = 30.0

phip   = (180.0 -   phi) * np.pi / 180.0
thetap = (180.0 - theta) * np.pi / 180.0

sinp, sint = np.sin(phip), np.sin(thetap)
cosp, cost = np.cos(phip), np.cos(thetap)

# rot = np.array([[cosp, -sinp * cost, sinp * sint], [sinp, cosp * cost, -cosp * sinp], [0.0, sint, cost]])
rot = np.array([[cosp, sinp, 0.0], [-sinp * cost, cosp * cost, sint], [sinp * sint, -cosp * sint, cost]])
print("rotation matrix")
print(rot)

# rot: disk rotation axis ==>> inclination in Kirihara et al. (2017a)
ini = np.array([0.0, 0.0, 1.0])# rotation axis of M31's disk in (X, Y, Z)-coordinate
fin = np.dot(rot, ini)# inclined disk (of the satellite)
print("rotation axis of the satellite in the disk orthogonal frame")
print(fin)

axis = np.cross(ini, fin)
sintheta = np.sqrt(np.dot(axis, axis))
costheta = np.dot(ini, fin)
axis = axis / np.linalg.norm(axis)
print("rotation axis of the rotation:")
print(axis)
print("sine and cosine of the rotation from disk frame to inclined disk")
print(sintheta, costheta)
print("rotation angle (disk frame to inclined disk):")
print(np.arcsin(sintheta), np.arccos(costheta))
print("sin(acos(costheta)) and cos(asin(sintheta)):")
print(np.sin(np.arccos(costheta)), np.cos(np.arcsin(sintheta)))
print("sin(asin(sintheta)) and cos(acos(costheta)):")
print(np.sin(np.arcsin(sintheta)), np.cos(np.arccos(costheta)))
