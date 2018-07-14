# RaDec --> M31 standard coordinate: https://github.com/astropy/astropy/blob/master/docs/coordinates/matchsep.rst "Sky Offset" Frames
import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.coordinates import angles

# M31 location (taken from NED, reference is 2010ApJS..189...37E)
m31 = SkyCoord(ra = "00h42m44.3503s", dec = "+41d16m08.634s", frame = "icrs")
# m31 = SkyCoord(ra = 10.684793 * u.deg, dec = 41.269065 * u.deg, frame = "icrs")
# m31 = SkyCoord.from_name("M31")
# print(m31)
# print(m31.galactic)

# comparison: galaxies taken from NED
galaxy = {
    "NGC 147": SkyCoord(ra = "00h33m12.12s", dec = "+48d30m31.5s", frame = "icrs"),
    "M 33   ": SkyCoord(ra = "01h33m50.89s", dec = "+30d39m36.8s", frame = "icrs"),
    "And VII": SkyCoord(ra = "23h26m31.74s", dec = "+50d40m32.6s", frame = "icrs")
}
for data_label, data in galaxy.items():
    # inm31 = data.transform_to(m31.skyoffset_frame())
    # xi, eta = inm31.lon.degree, inm31.lat.degree
    # print(data_label, ":", xi, eta)
    sep = data.separation(m31).degree
    pos = m31.position_angle(data)
    xi  = sep * np.sin(pos)
    eta = sep * np.cos(pos)
    print(data_label, ":", xi, eta)

# comparison: Table 1 of Martin et al. (2013)
martin13_tab1 = {
    "Detection 24": SkyCoord(ra = "01h27m40s", dec = "+28d05m53s", frame = "icrs"),
    "Detection 62": SkyCoord(ra = "01h31m38s", dec = "+29d18m15s", frame = "icrs"),
    "Detection  6": SkyCoord(ra = "23h54m43s", dec = "+42d30m28s", frame = "icrs"),
    "Detection  9": SkyCoord(ra = "00h36m21s", dec = "+49d36m40s", frame = "icrs")
}
# martin13_tab1_d24 = SkyCoord(ra = "01h27m40s", dec = "+28d05m53s", frame = "icrs")# X0 = 10.2 deg., Y0 = -12.9 deg. # And XXII
# martin13_tab1_d62 = SkyCoord(ra = "01h31m38s", dec = "+29d18m15s", frame = "icrs")# X0 = 11.0 deg., Y0 = -11.6 deg. # Overlap with M33 stellar debris
# martin13_tab1_d06 = SkyCoord(ra = "23h54m43s", dec = "+42d30m28s", frame = "icrs")# X0 = -8.9 deg., Y0 = 1.9 deg. # And XXI
# martin13_tab1_d09 = SkyCoord(ra = "00h36m21s", dec = "+49d36m40s", frame = "icrs")# X0 = -1.0 deg., Y0 = 8.4 deg. # Cas II / And XXX
# # martin13_tab1_d24 = SkyCoord(ra = "01h27m40s", dec = "+28d05m25s", frame = "icrs")# X0 = 10.2 deg., Y0 = -12.9 deg. # And XXII (RaDec is taken from NED)
# # martin13_tab1_d06 = SkyCoord(ra = "23h54m47.7s", dec = "+42d28m15s", frame = "icrs")# X0 = -8.9 deg., Y0 = 1.9 deg. # And XXI (RaDec is taken from NED)
# # martin13_tab1_d09 = SkyCoord(ra = "00h36m34.90s", dec = "+49d38m48.0s", frame = "icrs")# X0 = -1.0 deg., Y0 = 8.4 deg. # Cas II / And XXX (RaDec is taken from NED)
for data_label, data in martin13_tab1.items():
    # inm31 = data.transform_to(m31.skyoffset_frame())
    # xi, eta = inm31.lon.degree, inm31.lat.degree
    # print(data_label, ":", xi, eta)
    sep = data.separation(m31).degree
    pos = m31.position_angle(data)
    xi  = sep * np.sin(pos)
    eta = sep * np.cos(pos)
    print(data_label, ":", xi, eta)


# observation fields in Komiyama et al. (2018): taken from Table 1
komiyama18_tab1 = {
    "f003": SkyCoord(ra = "00h16m31.7s", dec = "+44d43m30s", frame = "icrs"),
    "f004": SkyCoord(ra = "00h10m04.7s", dec = "+45d27m47s", frame = "icrs"),
    "f009": SkyCoord(ra = "00h10m24.5s", dec = "+46d49m07s", frame = "icrs"),
    "f022": SkyCoord(ra = "00h16m04.2s", dec = "+43d22m15s", frame = "icrs"),
    "f023": SkyCoord(ra = "00h09m45.8s", dec = "+44d06m28s", frame = "icrs")
}
for data_label, data in komiyama18_tab1.items():
    # inm31 = data.transform_to(m31.skyoffset_frame())
    # xi, eta = inm31.lon.degree, inm31.lat.degree
    # print(data_label, ":", xi, eta)
    sep = data.separation(m31).degree
    pos = m31.position_angle(data)
    xi  = sep * np.sin(pos)
    eta = sep * np.cos(pos)
    print(data_label, ":", xi, eta)


# distances to the NW Stream in Komiyama et al. (2018): taken from Table 4
komiyama18_tab4 = {
    "Stream 1": SkyCoord(ra = 2.06 * u.deg, dec = 46.90 * u.deg, frame = "icrs"),
    "Stream 2": SkyCoord(ra = 3.00 * u.deg, dec = 45.45 * u.deg, frame = "icrs"),
    "Stream 3": SkyCoord(ra = 3.78 * u.deg, dec = 44.25 * u.deg, frame = "icrs"),
    "Stream 4": SkyCoord(ra = 4.52 * u.deg, dec = 43.10 * u.deg, frame = "icrs")
}
for data_label, data in komiyama18_tab4.items():
    # inm31 = data.transform_to(m31.skyoffset_frame())
    # xi, eta = inm31.lon.degree, inm31.lat.degree
    # print(data_label, ":", xi, eta)
    sep = data.separation(m31).degree
    pos = m31.position_angle(data)
    xi  = sep * np.sin(pos)
    eta = sep * np.cos(pos)
    print(data_label, ":", xi, eta)


# locations of globular clusters in Veljanoski et al. (2014): taken from Table 4
# PAndAS-04: V_hel = -397 \pm  7 km/s
# PAndAS-09: V_hel = -444 \pm 21 km/s
# PAndAS-10: V_hel = -435 \pm 10 km/s
# PAndAS-11: V_hel = -447 \pm 13 km/s
# PAndAS-12: V_hel = -472 \pm  5 km/s
veljanoski14_tab4 = {
    "PAndAS-04": SkyCoord(ra = "00h04m42.9s", dec = "+47d21m42s", frame = "icrs"),
    "PAndAS-09": SkyCoord(ra = "00h12m54.6s", dec = "+45d05m55s", frame = "icrs"),
    "PAndAS-10": SkyCoord(ra = "00h13m38.6s", dec = "+45d11m11s", frame = "icrs"),
    "PAndAS-11": SkyCoord(ra = "00h14m55.6s", dec = "+44d37m16s", frame = "icrs"),
    "PAndAS-12": SkyCoord(ra = "00h17m40.0s", dec = "+43d18m39s", frame = "icrs")
}
for data_label, data in veljanoski14_tab4.items():
    sep = data.separation(m31).degree
    pos = m31.position_angle(data)
    xi  = sep * np.sin(pos)
    eta = sep * np.cos(pos)
    print(data_label, ":", xi, eta)
