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


# implementation C
# l2s and s2l are provided by Kirihara-kun
l2s = np.array([
    [ 0.27136106, 0.61941386, 0.73667473],
    [ 0.27499914, 0.68358978, -0.67607727],
    [-0.92235495, 0.38604596, 0.01516137]
])
s2l = np.array([
    [ 0.27136106, 0.27499914, -0.92235495],
    [ 0.61941386, 0.68358978, 0.38604596],
    [ 0.73667473, -0.67607727, 0.01516137]
])
# orbit in Fardal et al. (2007)
pos = np.array([-34.75,  19.37, -13.99])# unit is kpc
vel = np.array([ 67.34, -26.12,  13.50])# unit is km/s
vec = np.cross(pos, vel)# L in observed frame
print ("orbital angular momentum in observed frame:")
print(vec)
orb = vec / np.linalg.norm(vec)
spn = np.dot(l2s, orb)
# spn = np.dot(s2l, orb)
print("rotation axis of the satellite in the observed frame:")
print(spn)
# rot4: observed frame (x, y, z) ==>> M31 disk frame (X, Y, Z)
# inv4: M31 disk frame (X, Y, Z) ==>> observed frame (x, y, z)
rot4, inv4 = m31.set_rotation_matrix()
fin = np.dot(rot4, spn)# spin axis in M31 disk frame
print("rotation axis of the satellite in the M31 disk frame:")
print(fin)
print("rotation axis of the satellite in the observed frame:")
print(np.dot(inv4, fin))
ini = np.array([0.0, 0.0, 1.0])
out = np.cross(ini, fin)
sint = np.sqrt(np.dot(out, out))
cost = np.dot(ini, fin)
axis = out / np.linalg.norm(out)
rot5, inv5 = m31.set_rodrigues_matrix(axis, sint, cost)
print("rotation matrix [(0, 0, 1) in M31 disk frame to Kirihara's rotation in M31 disk frame]:")
print(rot5)

print("spin axis in the observed frame:")
print(np.dot(inv4, np.dot(rot5, np.array([0.0, 0.0, 1.0]))))


# implementation D
l2s_obs = np.array([
    [-0.01358629, -0.74585316, -0.66597183],
    [0.79317924, 0.39750387, -0.46136469],
    [0.60883669, -0.53450326, 0.58619464]
])
# orbit in Fardal et al. (2007)
pos = np.array([-34.75,  19.37, -13.99])# unit is kpc
vel = np.array([ 67.34, -26.12,  13.50])# unit is km/s
vec = np.cross(pos, vel)# L in observed frame
print ("orbital angular momentum in observed frame:")
print(vec)
spin_obs = np.dot(l2s_obs, vec)
spin_obs_dir = spin_obs / np.linalg.norm(spin_obs)
print ("spin angular momentum in observed frame:")
print(spin_obs_dir)
obs2disk, disk2obs = m31.set_rotation_matrix()
spin_disk = np.dot(obs2disk, spin_obs_dir)
ini = np.array([0.0, 0.0, 1.0])
out = np.cross(ini, spin_disk)
sint = np.sqrt(np.dot(out, out))
cost = np.dot(ini, spin_disk)
axis = out / np.linalg.norm(out)
rot6, inv6 = m31.set_rodrigues_matrix(axis, sint, cost)
print("rotation matrix [(0, 0, 1) in M31 disk frame to Kirihara's rotation in M31 disk frame]:")
print(rot6)
print("spin angular momentum in observed frame:")
print(np.dot(disk2obs, np.dot(rot6, np.array([0.0, 0.0, 1.0]))))
