import numpy as np

dist = 137.0 # kpc

num = 5
name = ["GC1", "GC2", "GC3", "GC4", "GC5"]
mass = [4.57, 5.26, 5.56, 5.12, 5.25] # log_10 M
r0 = [14.4, 8.7, 1.4, 1.0, 4.2] # arcsec
rt = [91, 114, 95, 66, 76] # arcsec

Nhalo = 4
halo = ["cusp", "core1", "core2", "core3"]
Menc = [4.8882e+8, 1.7268e+9, 9.2246e+8, 4.9936e+8] # Msun (r < 1 kpc)

for ii in range(num):
    print(name[ii] + ":")
    mtot = 10.0 ** mass[ii]
    print("mass = {:.6e} Msun".format(mtot))

    print(" core radius = {:.6e} (pc)".format(1.0e+3 * dist * np.tan(np.deg2rad(r0[ii] / 3600))))
    print("tidal radius = {:.6e} (pc)".format(1.0e+3 * dist * np.tan(np.deg2rad(rt[ii] / 3600))))
    print("concentration parameter = {:.6f}".format(np.log10(np.tan(np.deg2rad(rt[ii] / 3600)) / np.tan(np.deg2rad(r0[ii] / 3600)))))

    for jj in range(Nhalo):
        rHill = (mtot / (3.0 * Menc[jj])) ** (1.0 / 3.0) # * 1 kpc
        print("Hill radius for {:<} is {:.6e} (pc)".format(halo[jj], 1.0e+3 * rHill))


print("GC0:")
print("mass = {:.6e} Msun".format(1.82e+5))
print("core radius = {:.6e} (pc)".format(1.0e+3 * 0.01051))
for jj in range(Nhalo):
    rHill = (1.82e+5 / (3.0 * Menc[jj])) ** (1.0 / 3.0) # * 1 kpc
    print("Hill radius for {:<} is {:.6e} (pc)".format(halo[jj], 1.0e+3 * rHill))
