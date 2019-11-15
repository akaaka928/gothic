# M31 standard coordinate --> RaDec (inverse of py/coord.py)
import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.coordinates import angles


# M31 location (taken from NED, reference is 2010ApJS..189...37E)
m31 = SkyCoord(ra = '00h42m44.3503s', dec = '+41d16m08.634s', frame = 'icrs')


def std2radec(xi, eta):
    sep = np.sqrt(xi ** 2 + eta ** 2)
    pos = np.arctan2(xi, eta)
    # print(sep, pos)
    gal = m31.directional_offset_by(pos, sep * u.degree)
    # print(gal) # print RaDec in units of degree
    # print(gal.to_string('hmsdms')) # print RaDec in hms and dms format
    return gal


def report(name, gal):
    print()
    print(name, ':', gal.to_string('hmsdms'))
    print(gal)
    print(gal.galactic)


# test (NGC 147)
# 'NGC 147': SkyCoord(ra = '00h33m12.12s', dec = '+48d30m31.5s', frame = 'icrs')
# xi_ngc147, eta_ngc147 = -1.58359005307475, 7.262465927867442
# std2radec(xi_ngc147, eta_ngc147)
# ngc147 = SkyCoord(ra = '00h33m12.12s', dec = '+48d30m31.5s', frame = 'icrs')
# sep = m31.separation(ngc147).degree
# pos = m31.position_angle(ngc147)
# print(ngc147)
# print(sep, pos)
# print(sep * np.sin(pos), sep * np.cos(pos))


# # show M31 info
# print('M31:', m31, m31.galactic, m31.to_string('decimal'), m31.to_string('hmsdms'))
report('M31', m31)


# observation fields in Komiyama et al. (2018): taken from Table 1
komiyama18_tab1 = {
    'f003': SkyCoord(ra = '00h16m31.7s', dec = '+44d43m30s', frame = 'icrs'),
    'f004': SkyCoord(ra = '00h10m04.7s', dec = '+45d27m47s', frame = 'icrs'),
    'f009': SkyCoord(ra = '00h10m24.5s', dec = '+46d49m07s', frame = 'icrs'),
    'f022': SkyCoord(ra = '00h16m04.2s', dec = '+43d22m15s', frame = 'icrs'),
    'f023': SkyCoord(ra = '00h09m45.8s', dec = '+44d06m28s', frame = 'icrs')
}
for data_label, data in komiyama18_tab1.items():
    # sep = m31.separation(data).degree
    # pos = m31.position_angle(data)
    # xi  = sep * np.sin(pos)
    # eta = sep * np.cos(pos)
    report(data_label, data)


# candidates of additional HSC observations (Northern fields)
# north_xi[0] = -6.1; north_eta[0] = 6.9; north_rad[0] = 0.75 # 0.75 degree is the Subaru/HSC FoV
# north_xi[1] = -6.3; north_eta[1] = 8.2; north_rad[1] = 0.75 # 0.75 degree is the Subaru/HSC FoV
report('N0', std2radec(-6.1, 6.9))
report('N1', std2radec(-6.3, 8.2))


# candidates of additional HSC/PFS observations (Southern fields)
# south_xi[0] = -4.05; south_eta[0] = 1.25; south_rad[0] = 0.75 # 0.75 degree is the Subaru/HSC FoV
# south_xi[1] = -3.2; south_eta[1] = 0.2; south_rad[1] = 0.75 # 0.75 degree is the Subaru/HSC FoV
report('S0', std2radec(-4.05, 1.25))
report('S1', std2radec(-3.2, 0.2))


# location of the progenitor core
