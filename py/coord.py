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


# distances to the NW Stream in Komiyama et al. (2018): taken from Table 4
komiyama18_tab4 = {
    "Stream 1": SkyCoord(ra = 2.06 * u.deg, dec = 46.90 * u.deg, frame = "icrs"),
    "Stream 2": SkyCoord(ra = 3.00 * u.deg, dec = 45.45 * u.deg, frame = "icrs"),
    "Stream 3": SkyCoord(ra = 3.78 * u.deg, dec = 44.25 * u.deg, frame = "icrs"),
    "Stream 4": SkyCoord(ra = 4.52 * u.deg, dec = 43.10 * u.deg, frame = "icrs")
}
for data_label, data in komiyama18_tab4.items():
    sep = data.separation(m31).degree
    pos = m31.position_angle(data)
    xi  = sep * np.sin(pos)
    eta = sep * np.cos(pos)
    print(data_label, ":", xi, eta)


# Stream fields in Komiyama et al. (2018): given by Komiyama-san
komiyama18_field = {
    "NE0": SkyCoord(ra = 1.3  * u.deg, dec = 48.8  * u.deg, frame = "icrs"),
    "NE1": SkyCoord(ra = 5.3  * u.deg, dec = 42.8  * u.deg, frame = "icrs"),
    "SW0": SkyCoord(ra = 1.3  * u.deg, dec = 47.32 * u.deg, frame = "icrs"),
    "SW1": SkyCoord(ra = 5.3  * u.deg, dec = 40.92 * u.deg, frame = "icrs"),
    "s12": SkyCoord(ra = 2.58125 * u.deg, dec = 46.1  * u.deg, frame = "icrs"),
    "s23": SkyCoord(ra = 3.420833333 * u.deg, dec = 44.8  * u.deg, frame = "icrs"),
    "s34": SkyCoord(ra = 4.13125 * u.deg, dec = 43.7  * u.deg, frame = "icrs")
}
for data_label, data in komiyama18_field.items():
    sep = data.separation(m31).degree
    pos = m31.position_angle(data)
    xi  = sep * np.sin(pos)
    eta = sep * np.cos(pos)
    print(data_label, ":", xi, eta)


# edge of Stream fields
NE0 = SkyCoord(ra = 1.3  * u.deg, dec = 48.8  * u.deg, frame = "icrs")
NE1 = SkyCoord(ra = 5.3  * u.deg, dec = 42.8  * u.deg, frame = "icrs")
SW0 = SkyCoord(ra = 1.3  * u.deg, dec = 47.32 * u.deg, frame = "icrs")
SW1 = SkyCoord(ra = 5.3  * u.deg, dec = 40.92 * u.deg, frame = "icrs")
NE0_sep = NE0.separation(m31).degree
NE0_pos = m31.position_angle(NE0)
NE0_xi  = NE0_sep * np.sin(NE0_pos)
NE0_eta = NE0_sep * np.cos(NE0_pos)
NE1_sep = NE1.separation(m31).degree
NE1_pos = m31.position_angle(NE1)
NE1_xi  = NE1_sep * np.sin(NE1_pos)
NE1_eta = NE1_sep * np.cos(NE1_pos)
SW0_sep = SW0.separation(m31).degree
SW0_pos = m31.position_angle(SW0)
SW0_xi  = SW0_sep * np.sin(SW0_pos)
SW0_eta = SW0_sep * np.cos(SW0_pos)
SW1_sep = SW1.separation(m31).degree
SW1_pos = m31.position_angle(SW1)
SW1_xi  = SW1_sep * np.sin(SW1_pos)
SW1_eta = SW1_sep * np.cos(SW1_pos)
# print(NE0_xi, NE0_eta)
# print(NE1_xi, NE1_eta)
# print(SW0_xi, SW0_eta)
# print(SW1_xi, SW1_eta)
NE_slope     = (NE1_eta          - NE0_eta         ) / (NE1_xi - NE0_xi)
NE_intercept = (NE0_eta * NE1_xi - NE1_eta * NE0_xi) / (NE1_xi - NE0_xi)
SW_slope     = (SW1_eta          - SW0_eta         ) / (SW1_xi - SW0_xi)
SW_intercept = (SW0_eta * SW1_xi - SW1_eta * SW0_xi) / (SW1_xi - SW0_xi)
print("NE edge: ", "slope = ", NE_slope, ", intercept = ", NE_intercept)
print("SW edge: ", "slope = ", SW_slope, ", intercept = ", SW_intercept)
sep12 = SkyCoord(ra = 2.58125 * u.deg, dec = 46.1  * u.deg, frame = "icrs")
sep23 = SkyCoord(ra = 3.420833333 * u.deg, dec = 44.8  * u.deg, frame = "icrs")
sep34 = SkyCoord(ra = 4.13125 * u.deg, dec = 43.7  * u.deg, frame = "icrs")
sep12_sep = sep12.separation(m31).degree
sep12_pos = m31.position_angle(sep12)
sep12_xi  = sep12_sep * np.sin(sep12_pos)
sep12_eta = sep12_sep * np.cos(sep12_pos)
sep23_sep = sep23.separation(m31).degree
sep23_pos = m31.position_angle(sep23)
sep23_xi  = sep23_sep * np.sin(sep23_pos)
sep23_eta = sep23_sep * np.cos(sep23_pos)
sep34_sep = sep34.separation(m31).degree
sep34_pos = m31.position_angle(sep34)
sep34_xi  = sep34_sep * np.sin(sep34_pos)
sep34_eta = sep34_sep * np.cos(sep34_pos)
print("edges of boundary between Stream fields 1 & 2:", (sep12_eta - NE_intercept) / NE_slope, sep12_eta, ", and", (sep12_eta - SW_intercept) / SW_slope, sep12_eta)
print("edges of boundary between Stream fields 2 & 3:", (sep23_eta - NE_intercept) / NE_slope, sep23_eta, ", and", (sep23_eta - SW_intercept) / SW_slope, sep23_eta)
print("edges of boundary between Stream fields 3 & 4:", (sep34_eta - NE_intercept) / NE_slope, sep34_eta, ", and", (sep34_eta - SW_intercept) / SW_slope, sep34_eta)


# observation fields in Komiyama et al. (2018): taken from Table 1
komiyama18_tab1 = {
    "f003": SkyCoord(ra = "00h16m31.7s", dec = "+44d43m30s", frame = "icrs"),
    "f004": SkyCoord(ra = "00h10m04.7s", dec = "+45d27m47s", frame = "icrs"),
    "f009": SkyCoord(ra = "00h10m24.5s", dec = "+46d49m07s", frame = "icrs"),
    "f022": SkyCoord(ra = "00h16m04.2s", dec = "+43d22m15s", frame = "icrs"),
    "f023": SkyCoord(ra = "00h09m45.8s", dec = "+44d06m28s", frame = "icrs")
}
rad = 0.75# degree
for data_label, data in komiyama18_tab1.items():
    sep = data.separation(m31).degree
    pos = m31.position_angle(data)
    xi  = sep * np.sin(pos)
    eta = sep * np.cos(pos)
    print(data_label, ":", xi, eta)

    # find points of intersection with Subaru/HSC fields and NW Stream fields
    # 1. NE edge
    bb = NE_intercept + NE_slope * xi - eta
    common = (1.0 + NE_slope ** 2)
    discriminant = common * (rad ** 2) - (bb ** 2)
    if discriminant >= 0:
        val = np.sqrt(discriminant)
        print("NE edge with", data_label, ":", xi - (NE_slope * bb - val) / common, eta + (bb + NE_slope * val) / common)
        print("NE edge with", data_label, ":", xi - (NE_slope * bb + val) / common, eta + (bb - NE_slope * val) / common)

    # 2. SW edge
    bb = SW_intercept + SW_slope * xi - eta
    common = (1.0 + SW_slope ** 2)
    discriminant = common * (rad ** 2) - (bb ** 2)
    if discriminant >= 0:
        val = np.sqrt(discriminant)
        print("SW edge with", data_label, ":", xi - (SW_slope * bb - val) / common, eta + (bb + SW_slope * val) / common)
        print("SW edge with", data_label, ":", xi - (SW_slope * bb + val) / common, eta + (bb - SW_slope * val) / common)
