# set initial location and velocity of the orbitin satellite in the M31's disk coordinate
import numpy as np
import m31 as m31


# print information on M31
print("distance to M31 is", m31.distance(), "kpc")
print("vM31_x is", m31.vx(), "km/s")
print("vM31_y is", m31.vy(), "km/s")
print("vM31_z is", m31.vz(), "km/s")

# prepare rotation matrix
rot, inv = m31.set_rotation_matrix()
print("rotation matrix [observed frame to M31 disk frame]:")
print(rot)
print("rotation matrix [M31 disk frame to observed frame]:")
print(inv)
# Fardal orbit in M31 disk frame
print("Fardal orbit: initial location in observed frame (kpc):", np.dot(rot, np.array([-34.75, 19.37, -13.99])))
print("Fardal orbit: initial velocity in M31 disk frame (km/s):", np.dot(rot, np.array([67.34, -26.12, 13.50])))


# # initial location
# xi, eta, D = -2.218220e+01, -6.034168e-01, 8.600729e+02# given by Kirihara-kun (private communication, original value in Komiyama et al. 2018)
# xx, yy, zz = m31.cartesian_coordinate(xi, eta, D)
# xi0, eta0, D0 = m31.standard_coordinate(xx, yy, zz)
# print("initial location in standard coordinate (given by Kirihara-kun)  :", xi, eta, D)
# print("initial location in standard coordinate (std -> cartesian -> std):", xi0, eta0, D0)
# print("initial location in observed frame (kpc):", xx, yy, zz)
# print("initial location in M31 disk frame (kpc):", np.dot(rot, np.array([xx, yy, zz])))


# # initial velocity
# vxi, veta, vlos = -9.877951e+00, -2.336741e+01, -3.021508e+02# given by Kirihara-kun (private communication, original value in Komiyama et al. 2018)
# vx, vy, vz = m31.cartesian_velocity(vxi, veta, vlos, xx, yy, zz)
# vxi0, veta0, vlos0 = m31.observed_velocity(xx, yy, zz, vx, vy, vz)
# print("initial velocity in standard coordinate (given by Kirihara-kun)  :", vxi, veta, vlos)
# print("initial velocity in standard coordinate (std -> cartesian -> std):", vxi0, veta0, vlos0)
# print("initial velocity in observed frame (km/s):", vx, vy, vz)
# print("initial velocity in M31 disk frame (km/s):", np.dot(rot, np.array([vx, vy, vz])))


# # initial location in Kirihara coordinate
# rot_k2m = np.array([
#     [0.0, 1.0, 0.0],
#     [1.0, 0.0, 0.0],
#     [0.0, 0.0, -1.0]
# ])
# XYZ_kiri = np.array([1.895850e+02, -2.556563e+01, 2.378530e+02])
# XYZ_disk = np.dot(rot_k2m, XYZ_kiri)
# print("trial version: initial location in M31 disk frame (kpc)", XYZ_disk)
# print("trial version: initial location in observed frame (kpc)", np.dot(inv, XYZ_disk))
# tmp = np.dot(inv, XYZ_disk)
# print("trial version: initial location in standard coordinate:", m31.standard_coordinate(tmp[0], tmp[1], tmp[2]))


# initial condition in test simulation
POS = np.array([2.556563e+1, -1.895850e+2, -2.378530e+2])
VEL = np.array([2.015646e+1, -2.068789e+1, 1.864773e+1])
print("initial location in M31 disk frame (kpc):", POS)
pos = np.dot(inv, POS)
print("initial location in observed frame (kpc):", pos)
xi, eta, D = m31.standard_coordinate(pos[0], pos[1], pos[2])
print("initial location in standard coordinate (deg., deg., kpc):", xi, eta, D)

print("initial velocity in M31 disk frame (km/s):", VEL)
vel = np.dot(inv, VEL)
print("initial velocity in observed frame (km/s):", vel)
vxi, veta, vlos = m31.observed_velocity(pos[0], pos[1], pos[2], vel[0], vel[1], vel[2])
print("initial velocity in standard coordinate (km/s):", vxi, veta, vlos)
