import numpy as np

zm31 = 776.0# same with Komiyama et al. (2018); median of NED
# zm31 = 773.0# M31 distance in Conn et al. (2016)

year = 3.15576e+7# year in units of sec
kpc2km = 3.0856775975e+21 * 1.0e-5

# vm31x = zm31 * kpc2km * np.tan(np.deg2rad( 49.0 * 1.0e-6 / 3600.0)) / year# Gaia DR2 + HST (van der Marel et al. 2018)
# vm31y = zm31 * kpc2km * np.tan(np.deg2rad(-38.0 * 1.0e-6 / 3600.0)) / year# Gaia DR2 + HST (van der Marel et al. 2018)
vm31x = 0.0
vm31y = 0.0
vm31z = -300.0

def distance():
    return zm31

def vx():
    return vm31x

def vy():
    return vm31y

def vz():
    return vm31z

def set_rodrigues_matrix(axis, sinx, cosx):
    rot = np.array([
        [(1.0 - cosx) * axis[0] * axis[0] + cosx          , (1.0 - cosx) * axis[0] * axis[1] - sinx * axis[2], (1.0 - cosx) * axis[2] * axis[0] + sinx * axis[1]],
        [(1.0 - cosx) * axis[0] * axis[1] + sinx * axis[2], (1.0 - cosx) * axis[1] * axis[1] + cosx          , (1.0 - cosx) * axis[1] * axis[2] - sinx * axis[0]],
        [(1.0 - cosx) * axis[2] * axis[0] - sinx * axis[1], (1.0 - cosx) * axis[1] * axis[2] + sinx * axis[0], (1.0 - cosx) * axis[2] * axis[2] + cosx          ]
    ])
    inv = rot.T

    return rot, inv


# M31 disk frame: Z-axis is rotation axis, Y-axis is major axis (points north-east), X-axis is minor axis (points north-west)
# observed frame: z-axis is line-of-sight, y-axis points north, x-axis points east
# rot: observed frame (x, y, z) ==>> M31 disk frame (X, Y, Z)
# inv: M31 disk frame (X, Y, Z) ==>> observed frame (x, y, z)
def set_rotation_matrix():
    theta = -37.0 * np.pi / 180.0# position angle
    phi   = -77.0 * np.pi / 180.0# inclination

    cost, cosp = np.cos(theta), np.cos(phi)
    sint, sinp = np.sin(theta), np.sin(phi)

    # 1st rotation: rotation axis is z-axis, rotation angle is theta
    rot1, inv1 = set_rodrigues_matrix(np.array([0.0, 0.0, 1.0]), sint, cost)

    # 2nd rotation: rotation axis is rotated y-axis (rot1 * (0, 1, 0)), rotation angle is phi
    # print(np.dot(rot1, np.array([0.0, 1.0, 0.0])))
    rot2, inv2 = set_rodrigues_matrix(np.dot(rot1, np.array([0.0, 1.0, 0.0])), sinp, cosp)

    # get rotation axis of M31's disk in observed frame
    spin = np.dot(rot2, np.dot(rot1, np.array([0.0, 0.0, -1.0])))
    # print("rotation axis of M31's disk in observed frame:")
    # print(spin)

    # rotate spin axis of M31's disk
    ini = spin
    fin = np.array([0.0, 0.0, 1.0])
    out = np.cross(ini, fin)
    sint = np.sqrt(np.dot(out, out))
    cost = np.dot(ini, fin)
    axis = out / np.linalg.norm(out)
    # print("rotation axis is:")
    # print(axis)
    # print("sintheta, costheta are:")
    # print(sint, cost)
    rot3, inv3 = set_rodrigues_matrix(axis, sint, cost)

    # major axis in M31 disk frame
    major_disk = np.array([0.0, 1.0, 0.0])
    major_obs = np.dot(rot2, np.dot(rot1, major_disk))
    major_tmp = np.dot(inv3, major_disk)

    rot4, inv4 = set_rodrigues_matrix(spin, -np.sqrt(np.dot(np.cross(major_tmp, major_obs), np.cross(major_tmp, major_obs))), np.dot(major_tmp, major_obs))

    inv = np.dot(rot4, inv3)
    rot = inv.T

    # print("X-axis in observed frame:")
    # print(np.dot(inv, np.array([1.0, 0.0, 0.0])))
    # print("Y-axis in observed frame:")
    # print(np.dot(inv, np.array([0.0, 1.0, 0.0])))
    # print("Z-axis in observed frame:")
    # print(np.dot(inv, np.array([0.0, 0.0, 1.0])))

    return rot, inv


def kpc2degree(kpc):
    degree = np.rad2deg(np.arctan(kpc / zm31))

    return degree


def standard_coordinate(xx, yy, zz):
    xi = np.rad2deg(np.arctan(xx / (zm31 + zz)))
    eta = np.rad2deg(np.arctan(yy / (zm31 + zz)))
    D = np.sqrt(xx ** 2 + yy ** 2 + (zm31 + zz) ** 2)

    return xi, eta, D


def observed_velocity(xx, yy, zz, vx, vy, vz):
    zz += zm31
    D = np.sqrt(xx ** 2 + yy ** 2 + zz ** 2)
    vx += vm31x
    vy += vm31y
    vz += vm31z
    vlos = (xx * vx + yy * vy + zz * vz) / D
    vxi  = (zz * vx - xx * vz) / D
    veta = (zz * vy - yy * vz) / D

    return vxi, veta, vlos


def cartesian_coordinate(xi, eta, D):
    tan_xi  = np.tan(np.deg2rad(xi))
    tan_eta = np.tan(np.deg2rad(eta))
    dist = D / np.sqrt(1.0 + tan_xi ** 2 + tan_eta ** 2)
    xx = dist * tan_xi
    yy = dist * tan_eta
    zz = dist - zm31

    return xx, yy, zz


def cartesian_velocity(vxi, veta, vlos, xx, yy, zz):
    zp = zz + zm31
    dist = np.sqrt(xx ** 2 + yy ** 2 + zp ** 2)
    vz = (zp * vlos - xx * vxi - yy * veta) / dist
    vx = (dist * vxi  + xx * vz) / zp
    vy = (dist * veta + yy * vz) / zp
    vx -= vm31x
    vy -= vm31y
    vz -= vm31z

    return vx, vy, vz


def disk_ellipse():
    num = 128 + 1
    xi  = [0] * num
    eta = [0] * num
    D   = [0] * num

    rot, inv = set_rotation_matrix()

    rad = 30.0
    dtheta = 2.0 * np.pi / (num - 1)
    for ii in range(num):
        theta = dtheta * ii
        obs = np.dot(inv, np.array([rad * np.cos(theta), rad * np.sin(theta), 0.0]))

        xi[ii], eta[ii], D[ii] = standard_coordinate(obs[0], obs[1], obs[2])
        # D[ii] -= zm31

    return num, xi, eta, D


def Fardal2007_shell():
    Neast = 8
    Eshell_xi  = [0] * Neast
    Eshell_eta = [0] * Neast

    Eshell_xi[0] = -0.56;    Eshell_eta[0] =  1.61
    Eshell_xi[1] = -1.03;    Eshell_eta[1] =  1.49
    Eshell_xi[2] = -1.45;    Eshell_eta[2] =  1.16
    Eshell_xi[3] = -1.58;    Eshell_eta[3] =  0.97
    Eshell_xi[4] = -1.78;    Eshell_eta[4] =  0.47
    Eshell_xi[5] = -1.87;    Eshell_eta[5] =  0.17
    Eshell_xi[6] = -1.92;    Eshell_eta[6] = -0.31
    Eshell_xi[7] = -1.95;    Eshell_eta[7] = -0.56

    Nwest = 7
    Wshell_xi  = [0] * Nwest
    Wshell_eta = [0] * Nwest

    Wshell_xi[0] = 1.95;    Wshell_eta[0] =  1.46
    Wshell_xi[1] = 1.96;    Wshell_eta[1] =  0.93
    Wshell_xi[2] = 1.91;    Wshell_eta[2] =  0.59
    Wshell_xi[3] = 1.86;    Wshell_eta[3] =  0.41
    Wshell_xi[4] = 1.78;    Wshell_eta[4] =  0.31
    Wshell_xi[5] = 1.70;    Wshell_eta[5] =  0.17
    Wshell_xi[6] = 1.45;    Wshell_eta[6] = -0.01

    return Neast, Eshell_xi, Eshell_eta, Nwest, Wshell_xi, Wshell_eta


# convert pixeld ID to (xi, eta)
def Conn2016_coord(xx, yy):
    # x = -15: xi = 6 degree
    # x = 403: xi = 0 degree
    # y = 131: eta = -6 degree
    # y = 689: eta =  2 degree
    # y = 550: eta = 0 degree
    xi0 = 403
    eta0 = 550
    pix2xi  = -0.014354066# (0 - 6) / (403 - (-15))
    pix2eta =  0.014336917# (2 - (-6)) / (689 - 131)

    xi  = (xx -  xi0) * pix2xi
    eta = (yy - eta0) * pix2eta

    return xi, eta


def GSS_obs_field():
    # # Table 1 in Font et al. (2006)
    # Nfield = 11
    # field_xi  = [0] * Nfield
    # field_eta = [0] * Nfield

    # field_xi[ 0] =  1.077;    field_eta[ 0] = -2.021# a3
    # field_xi[ 1] =  2.015;    field_eta[ 1] = -3.965# 1
    # field_xi[ 2] =  1.745;    field_eta[ 2] = -3.525# 2
    # field_xi[ 3] =  1.483;    field_eta[ 3] = -3.087# 3
    # field_xi[ 4] =  1.226;    field_eta[ 4] = -2.653# 4
    # field_xi[ 5] =  0.969;    field_eta[ 5] = -2.264# 5
    # field_xi[ 6] =  0.717;    field_eta[ 6] = -1.768# 6
    # field_xi[ 7] =  0.467;    field_eta[ 7] = -1.327# 7
    # field_xi[ 8] =  0.219;    field_eta[ 8] = -0.886# 8
    # field_xi[ 9] = -0.731;    field_eta[ 9] =  0.891# 12
    # field_xi[10] = -0.963;    field_eta[10] =  1.342# 13

    # read from Figure 3 of Conn et al. (2016)
    Nfield = 8
    field_xi  = np.zeros((Nfield, 4 + 1))
    field_eta = np.zeros((Nfield, 4 + 1))
    # M1
    field_xi[0][0], field_eta[0][0] = Conn2016_coord(286, 289)
    field_xi[0][1], field_eta[0][1] = Conn2016_coord(286, 257)
    field_xi[0][2], field_eta[0][2] = Conn2016_coord(238, 257)
    field_xi[0][3], field_eta[0][3] = Conn2016_coord(238, 289)
    # M2
    field_xi[1][0], field_eta[1][0] = Conn2016_coord(305, 320)
    field_xi[1][1], field_eta[1][1] = Conn2016_coord(305, 287)
    field_xi[1][2], field_eta[1][2] = Conn2016_coord(257, 287)
    field_xi[1][3], field_eta[1][3] = Conn2016_coord(257, 320)
    # M3
    field_xi[2][0], field_eta[2][0] = Conn2016_coord(324, 350)
    field_xi[2][1], field_eta[2][1] = Conn2016_coord(324, 318)
    field_xi[2][2], field_eta[2][2] = Conn2016_coord(275, 318)
    field_xi[2][3], field_eta[2][3] = Conn2016_coord(275, 350)
    # M4
    field_xi[3][0], field_eta[3][0] = Conn2016_coord(342, 381)
    field_xi[3][1], field_eta[3][1] = Conn2016_coord(342, 348)
    field_xi[3][2], field_eta[3][2] = Conn2016_coord(293, 348)
    field_xi[3][3], field_eta[3][3] = Conn2016_coord(293, 381)
    # M5
    field_xi[4][0], field_eta[4][0] = Conn2016_coord(359, 408)
    field_xi[4][1], field_eta[4][1] = Conn2016_coord(359, 375)
    field_xi[4][2], field_eta[4][2] = Conn2016_coord(311, 375)
    field_xi[4][3], field_eta[4][3] = Conn2016_coord(311, 408)
    # M6
    field_xi[5][0], field_eta[5][0] = Conn2016_coord(377, 442)
    field_xi[5][1], field_eta[5][1] = Conn2016_coord(377, 410)
    field_xi[5][2], field_eta[5][2] = Conn2016_coord(328, 410)
    field_xi[5][3], field_eta[5][3] = Conn2016_coord(328, 442)
    # M7
    field_xi[6][0], field_eta[6][0] = Conn2016_coord(395, 473)
    field_xi[6][1], field_eta[6][1] = Conn2016_coord(395, 441)
    field_xi[6][2], field_eta[6][2] = Conn2016_coord(346, 441)
    field_xi[6][3], field_eta[6][3] = Conn2016_coord(346, 473)
    # M8
    field_xi[7][0], field_eta[7][0] = Conn2016_coord(412, 504)
    field_xi[7][1], field_eta[7][1] = Conn2016_coord(412, 472)
    field_xi[7][2], field_eta[7][2] = Conn2016_coord(363, 472)
    field_xi[7][3], field_eta[7][3] = Conn2016_coord(363, 504)
    for ii in range(Nfield):
        field_xi[ii][4]  = field_xi[ii][0]
        field_eta[ii][4] = field_eta[ii][0]

    return Nfield, field_xi, field_eta


def GSS_distance():
    # NOTE
    # subfields GSSx.5 are skipped
    # GSS1 is skipped to remove overlap with M31 disk

    # Table 1 in Conn et al. (2016)
    Ngss = 9
    gss_xi  = [0] * Ngss
    gss_eta = [0] * Ngss
    gss_D   = [0] * Ngss
    gss_Dep = [0] * Ngss
    gss_Dem = [0] * Ngss

    gss_xi[0] = -0.047;    gss_eta[0] = -1.462;    gss_D[0] = 760.0;    gss_Dep[0] = 10.0;    gss_Dem[0] =  6.0# GSS2
    gss_xi[1] =  0.297;    gss_eta[1] = -1.937;    gss_D[1] = 787.0;    gss_Dep[1] =  5.0;    gss_Dem[1] =  7.0# GSS3
    gss_xi[2] =  0.641;    gss_eta[2] = -2.411;    gss_D[2] = 800.0;    gss_Dep[2] =  8.0;    gss_Dem[2] =  7.0# GSS4
    gss_xi[3] =  0.984;    gss_eta[3] = -2.885;    gss_D[3] = 821.0;    gss_Dep[3] =  7.0;    gss_Dem[3] =  9.0# GSS5
    gss_xi[4] =  1.328;    gss_eta[4] = -3.358;    gss_D[4] = 859.0;    gss_Dep[4] =  7.0;    gss_Dem[4] = 21.0# GSS6
    gss_xi[5] =  1.671;    gss_eta[5] = -3.830;    gss_D[5] = 826.0;    gss_Dep[5] = 18.0;    gss_Dem[5] =  8.0# GSS7
    gss_xi[6] =  2.015;    gss_eta[6] = -4.302;    gss_D[6] = 871.0;    gss_Dep[6] =  6.0;    gss_Dem[6] =  7.0# GSS8
    gss_xi[7] =  2.358;    gss_eta[7] = -4.772;    gss_D[7] = 844.0;    gss_Dep[7] =  8.0;    gss_Dem[7] =  7.0# GSS9
    gss_xi[8] =  2.701;    gss_eta[8] = -5.242;    gss_D[8] = 870.0;    gss_Dep[8] = 25.0;    gss_Dem[8] = 41.0# GSS10

    # for ii in range(Ngss):
    #     gss_D[ii] -= zm31

    gss_field_xi  = np.zeros((Ngss, 4 + 1))
    gss_field_eta = np.zeros((Ngss, 4 + 1))
    # GSS2
    gss_field_xi[0][0], gss_field_eta[0][0] = Conn2016_coord(479, 388)
    gss_field_xi[0][1], gss_field_eta[0][1] = Conn2016_coord(455, 358)
    gss_field_xi[0][2], gss_field_eta[0][2] = Conn2016_coord(333, 462)
    gss_field_xi[0][3], gss_field_eta[0][3] = Conn2016_coord(357, 493)
    # GSS3
    gss_field_xi[1][0], gss_field_eta[1][0] = Conn2016_coord(455, 358)
    gss_field_xi[1][1], gss_field_eta[1][1] = Conn2016_coord(431, 327)
    gss_field_xi[1][2], gss_field_eta[1][2] = Conn2016_coord(309, 432)
    gss_field_xi[1][3], gss_field_eta[1][3] = Conn2016_coord(333, 462)
    # GSS4
    gss_field_xi[2][0], gss_field_eta[2][0] = Conn2016_coord(431, 327)
    gss_field_xi[2][1], gss_field_eta[2][1] = Conn2016_coord(407, 297)
    gss_field_xi[2][2], gss_field_eta[2][2] = Conn2016_coord(285, 401)
    gss_field_xi[2][3], gss_field_eta[2][3] = Conn2016_coord(309, 432)
    # GSS5
    gss_field_xi[3][0], gss_field_eta[3][0] = Conn2016_coord(407, 297)
    gss_field_xi[3][1], gss_field_eta[3][1] = Conn2016_coord(383, 266)
    gss_field_xi[3][2], gss_field_eta[3][2] = Conn2016_coord(261, 371)
    gss_field_xi[3][3], gss_field_eta[3][3] = Conn2016_coord(285, 401)
    # GSS6
    gss_field_xi[4][0], gss_field_eta[4][0] = Conn2016_coord(383, 266)
    gss_field_xi[4][1], gss_field_eta[4][1] = Conn2016_coord(359, 236)
    gss_field_xi[4][2], gss_field_eta[4][2] = Conn2016_coord(237, 340)
    gss_field_xi[4][3], gss_field_eta[4][3] = Conn2016_coord(261, 371)
    # GSS7
    gss_field_xi[5][0], gss_field_eta[5][0] = Conn2016_coord(359, 236)
    gss_field_xi[5][1], gss_field_eta[5][1] = Conn2016_coord(335, 205)
    gss_field_xi[5][2], gss_field_eta[5][2] = Conn2016_coord(213, 310)
    gss_field_xi[5][3], gss_field_eta[5][3] = Conn2016_coord(237, 340)
    # GSS8
    gss_field_xi[6][0], gss_field_eta[6][0] = Conn2016_coord(335, 205)
    gss_field_xi[6][1], gss_field_eta[6][1] = Conn2016_coord(311, 175)
    gss_field_xi[6][2], gss_field_eta[6][2] = Conn2016_coord(189, 279)
    gss_field_xi[6][3], gss_field_eta[6][3] = Conn2016_coord(213, 310)
    # GSS9
    gss_field_xi[7][0], gss_field_eta[7][0] = Conn2016_coord(311, 175)
    gss_field_xi[7][1], gss_field_eta[7][1] = Conn2016_coord(287, 144)
    gss_field_xi[7][2], gss_field_eta[7][2] = Conn2016_coord(165, 249)
    gss_field_xi[7][3], gss_field_eta[7][3] = Conn2016_coord(189, 279)
    # GSS10
    gss_field_xi[8][0], gss_field_eta[8][0] = Conn2016_coord(287, 144)
    gss_field_xi[8][1], gss_field_eta[8][1] = Conn2016_coord(263, 114)
    gss_field_xi[8][2], gss_field_eta[8][2] = Conn2016_coord(141, 218)
    gss_field_xi[8][3], gss_field_eta[8][3] = Conn2016_coord(165, 249)
    for ii in range(Ngss):
        gss_field_xi[ii][4]  = gss_field_xi[ii][0]
        gss_field_eta[ii][4] = gss_field_eta[ii][0]

    return Ngss, gss_xi, gss_eta, gss_D, gss_Dep, gss_Dem, gss_field_xi, gss_field_eta


def StreamC_distance():
    # Table 2 in Conn et al. (2016)
    NstreamC = 4
    streamC_xi  = [0] * NstreamC
    streamC_eta = [0] * NstreamC
    streamC_D   = [0] * NstreamC
    streamC_Dep = [0] * NstreamC
    streamC_Dem = [0] * NstreamC

    streamC_xi[0] = 2.558;    streamC_eta[0] = -3.676;    streamC_D[0] = 809.0;    streamC_Dep[0] =  9.0;    streamC_Dem[0] = 15.0# C1
    streamC_xi[1] = 3.182;    streamC_eta[1] = -3.023;    streamC_D[1] = 854.0;    streamC_Dep[1] =  8.0;    streamC_Dem[1] = 28.0# C2
    streamC_xi[2] = 3.580;    streamC_eta[2] = -1.896;    streamC_D[2] = 819.0;    streamC_Dep[2] =  9.0;    streamC_Dem[2] = 56.0# C3
    streamC_xi[3] = 3.715;    streamC_eta[3] = -0.499;    streamC_D[3] = 831.0;    streamC_Dep[3] =  8.0;    streamC_Dem[3] = 20.0# C4

    # for ii in range(NstreamC):
    #     streamC_D[ii] -= zm31

    streamC_field_xi  = np.zeros((NstreamC, 4 + 1))
    streamC_field_eta = np.zeros((NstreamC, 4 + 1))
    # C1
    streamC_field_xi[0][0], streamC_field_eta[0][0] = Conn2016_coord(251, 316)
    streamC_field_xi[0][1], streamC_field_eta[0][1] = Conn2016_coord(235, 254)
    streamC_field_xi[0][2], streamC_field_eta[0][2] = Conn2016_coord(174, 285)
    streamC_field_xi[0][3], streamC_field_eta[0][3] = Conn2016_coord(224, 331)
    # C2
    streamC_field_xi[1][0], streamC_field_eta[1][0] = Conn2016_coord(224, 331)
    streamC_field_xi[1][1], streamC_field_eta[1][1] = Conn2016_coord(174, 285)
    streamC_field_xi[1][2], streamC_field_eta[1][2] = Conn2016_coord(138, 369)
    streamC_field_xi[1][3], streamC_field_eta[1][3] = Conn2016_coord(187, 391)
    # C3
    streamC_field_xi[2][0], streamC_field_eta[2][0] = Conn2016_coord(187, 391)
    streamC_field_xi[2][1], streamC_field_eta[2][1] = Conn2016_coord(138, 369)
    streamC_field_xi[2][2], streamC_field_eta[2][2] = Conn2016_coord(120, 472)
    streamC_field_xi[2][3], streamC_field_eta[2][3] = Conn2016_coord(167, 466)
    # C4
    streamC_field_xi[3][0], streamC_field_eta[3][0] = Conn2016_coord(167, 466)
    streamC_field_xi[3][1], streamC_field_eta[3][1] = Conn2016_coord(120, 472)
    streamC_field_xi[3][2], streamC_field_eta[3][2] = Conn2016_coord(122, 555)
    streamC_field_xi[3][3], streamC_field_eta[3][3] = Conn2016_coord(160, 625)
    for ii in range(NstreamC):
        streamC_field_xi[ii][4]  = streamC_field_xi[ii][0]
        streamC_field_eta[ii][4] = streamC_field_eta[ii][0]

    return NstreamC, streamC_xi, streamC_eta, streamC_D, streamC_Dep, streamC_Dem, streamC_field_xi, streamC_field_eta


def StreamD_distance():
    # Table 2 in Conn et al. (2016)
    NstreamD = 5
    streamD_xi  = [0] * NstreamD
    streamD_eta = [0] * NstreamD
    streamD_D   = [0] * NstreamD
    streamD_Dep = [0] * NstreamD
    streamD_Dem = [0] * NstreamD

    streamD_xi[0] = 2.174;    streamD_eta[0] = -2.142;    streamD_D[0] = 779.0;    streamD_Dep[0] =  7.0;    streamD_Dem[0] = 21.0# D1
    streamD_xi[1] = 2.728;    streamD_eta[1] = -1.423;    streamD_D[1] = 782.0;    streamD_Dep[1] = 32.0;    streamD_Dem[1] = 26.0# D2
    streamD_xi[2] = 2.947;    streamD_eta[2] = -0.579;    streamD_D[2] = 781.0;    streamD_Dep[2] = 26.0;    streamD_Dem[2] = 13.0# D3
    streamD_xi[3] = 3.097;    streamD_eta[3] =  0.469;    streamD_D[3] = 818.0;    streamD_Dep[3] = 50.0;    streamD_Dem[3] = 15.0# D4
    streamD_xi[4] = 2.932;    streamD_eta[4] =  1.198;    streamD_D[4] = 783.0;    streamD_Dep[4] = 13.0;    streamD_Dem[4] = 17.0# D5

    # for ii in range(NstreamD):
    #     streamD_D[ii] -= zm31

    streamD_field_xi  = np.zeros((NstreamD, 4 + 1))
    streamD_field_eta = np.zeros((NstreamD, 4 + 1))
    # D1
    streamD_field_xi[0][0], streamD_field_eta[0][0] = Conn2016_coord(290, 402)
    streamD_field_xi[0][1], streamD_field_eta[0][1] = Conn2016_coord(248, 360)
    streamD_field_xi[0][2], streamD_field_eta[0][2] = Conn2016_coord(202, 414)
    streamD_field_xi[0][3], streamD_field_eta[0][3] = Conn2016_coord(253, 440)
    # D2
    streamD_field_xi[1][0], streamD_field_eta[1][0] = Conn2016_coord(253, 440)
    streamD_field_xi[1][1], streamD_field_eta[1][1] = Conn2016_coord(202, 414)
    streamD_field_xi[1][2], streamD_field_eta[1][2] = Conn2016_coord(171, 478)
    streamD_field_xi[1][3], streamD_field_eta[1][3] = Conn2016_coord(233, 486)
    # D3
    streamD_field_xi[2][0], streamD_field_eta[2][0] = Conn2016_coord(233, 486)
    streamD_field_xi[2][1], streamD_field_eta[2][1] = Conn2016_coord(171, 478)
    streamD_field_xi[2][2], streamD_field_eta[2][2] = Conn2016_coord(160, 545)
    streamD_field_xi[2][3], streamD_field_eta[2][3] = Conn2016_coord(223, 541)
    # D4
    streamD_field_xi[3][0], streamD_field_eta[3][0] = Conn2016_coord(223, 541)
    streamD_field_xi[3][1], streamD_field_eta[3][1] = Conn2016_coord(160, 545)
    streamD_field_xi[3][2], streamD_field_eta[3][2] = Conn2016_coord(158, 624)
    streamD_field_xi[3][3], streamD_field_eta[3][3] = Conn2016_coord(213, 584)
    # D5
    streamD_field_xi[4][0], streamD_field_eta[4][0] = Conn2016_coord(213, 584)
    streamD_field_xi[4][1], streamD_field_eta[4][1] = Conn2016_coord(158, 624)
    streamD_field_xi[4][2], streamD_field_eta[4][2] = Conn2016_coord(224, 682)
    streamD_field_xi[4][3], streamD_field_eta[4][3] = Conn2016_coord(237, 612)
    for ii in range(NstreamD):
        streamD_field_xi[ii][4]  = streamD_field_xi[ii][0]
        streamD_field_eta[ii][4] = streamD_field_eta[ii][0]

    return NstreamD, streamD_xi, streamD_eta, streamD_D, streamD_Dep, streamD_Dem, streamD_field_xi, streamD_field_eta


def NWS_Komiyama2018():
    # Subaru/HSC files: Table 1 in Komiyama et al. (2018)
    Nhsc = 5
    hsc_xi  = [0] * Nhsc
    hsc_eta = [0] * Nhsc
    hsc_rad = [0] * Nhsc
    hsc_fov = 0.75# degree

    hsc_xi[0] = -4.653744156858361 ;    hsc_eta[0] = 3.635683362752261;    hsc_rad[0] = hsc_fov# f003
    hsc_xi[1] = -5.7227710294896745;    hsc_eta[1] = 4.47091006650441 ;    hsc_rad[1] = hsc_fov# f004
    hsc_xi[2] = -5.530826768687439 ;    hsc_eta[2] = 5.816784796173433;    hsc_rad[2] = hsc_fov# f009
    hsc_xi[3] = -4.842739589150247 ;    hsc_eta[3] = 2.290423176281251;    hsc_rad[3] = hsc_fov# f022
    hsc_xi[4] = -5.912437203950868 ;    hsc_eta[4] = 3.125050968157858;    hsc_rad[4] = hsc_fov# f023

    # Distances to NW Stream measured by RC stars: Table 4 in Komiyama et al. (2018)
    Nnws = 4
    nws_xi  = [0] * Nnws
    nws_eta = [0] * Nnws
    nws_D   = [0] * Nnws
    nws_Dep = [0] * Nnws
    nws_Dem = [0] * Nnws
    nws_mod = [0] * Nnws
    nws_ran = [0] * Nnws
    nws_sys = [0] * Nnws

    nws_xi[0] = -5.891761379732909;    nws_eta[0] = 5.934934874532188 ;    nws_mod[0] = 24.64;    nws_ran[0] = 0.196;    nws_sys[0] = 0.057# Stream 1
    nws_xi[1] = -5.388266383737322;    nws_eta[1] = 4.426258888399156 ;    nws_mod[1] = 24.62;    nws_ran[1] = 0.186;    nws_sys[1] = 0.057# Stream 2
    nws_xi[2] = -4.942619286035076;    nws_eta[2] = 3.18150303805077  ;    nws_mod[2] = 24.59;    nws_ran[2] = 0.183;    nws_sys[2] = 0.057# Stream 3
    nws_xi[3] = -4.49814447496522 ;    nws_eta[3] = 1.9926453169853198;    nws_mod[3] = 24.58;    nws_ran[3] = 0.185;    nws_sys[3] = 0.057# Stream 4

    for ii in range(Nnws):
        nws_D[ii] = 10.0 ** (nws_mod[ii] / 5.0 - 2.0)
        nws_Dep[ii] = 10.0 ** ((nws_mod[ii] + nws_ran[ii] + nws_sys[ii]) / 5.0 - 2.0) - nws_D[ii]
        nws_Dem[ii] = nws_D[ii] - 10.0 ** ((nws_mod[ii] - nws_ran[ii] - nws_sys[ii]) / 5.0 - 2.0)

    # Stream fields in Komiyama et al. (2018): Table 4 in Komiyama et al. (2018) and private communications with Komiyama-san
    # Stream 1
    Nnws_s1 = 5 + 64
    nws_s1_xi  = [0] * Nnws_s1
    nws_s1_eta = [0] * Nnws_s1
    nws_s1_xi[0], nws_s1_eta[0] = -5.7055741988006705, 6.546142963104346
    nws_s1_xi[1], nws_s1_eta[1] = -5.186841492207897 , 5.101772155015601
    nws_s1_xi[2], nws_s1_eta[2] = -5.907773581780884 , 5.101772155015601
    nws_s1_xi[3], nws_s1_eta[3] = -6.239719954692477 , 6.061673446044403
    nws_s1_min = np.pi - np.arcsin((nws_s1_eta[3] - hsc_eta[2]) / hsc_fov)
    nws_s1_max =         np.arccos((nws_s1_xi[0]  - hsc_xi[2] ) / hsc_fov)
    dtheta = (nws_s1_max - nws_s1_min) / 64
    for ii in range(64):
        nws_s1_xi[4 + ii]  = hsc_xi[2]  + hsc_fov * np.cos(nws_s1_min + float(ii) * dtheta)
        nws_s1_eta[4 + ii] = hsc_eta[2] + hsc_fov * np.sin(nws_s1_min + float(ii) * dtheta)
    # Stream 2
    Nnws_s2 = 5
    nws_s2_xi  = [0] * Nnws_s2
    nws_s2_eta = [0] * Nnws_s2
    nws_s2_xi[0], nws_s2_eta[0] = -5.186841492207897 , 5.101772155015601
    nws_s2_xi[1], nws_s2_eta[1] = -4.701958106244509 , 3.7516521430870844
    nws_s2_xi[2], nws_s2_eta[2] = -5.4408844910550105, 3.7516521430870844
    nws_s2_xi[3], nws_s2_eta[3] = -5.907773581780884 , 5.101772155015601
    # Stream 3
    Nnws_s3 = 5
    nws_s3_xi  = [0] * Nnws_s3
    nws_s3_eta = [0] * Nnws_s3
    nws_s3_xi[0], nws_s3_eta[0] = -4.701958106244509 , 3.7516521430870844
    nws_s3_xi[1], nws_s3_eta[1] = -4.292884147430157 , 2.612617574895511
    nws_s3_xi[2], nws_s3_eta[2] = -5.046991496988008 , 2.612617574895511
    nws_s3_xi[3], nws_s3_eta[3] = -5.4408844910550105, 3.7516521430870844
    # Stream 4
    Nnws_s4 = 5 + 64
    nws_s4_xi  = [0] * Nnws_s4
    nws_s4_eta = [0] * Nnws_s4
    nws_s4_xi[0], nws_s4_eta[0] = -4.6822224216441475, 1.5578016922010893
    nws_s4_xi[1], nws_s4_eta[1] = -5.046991496988008 , 2.612617574895511
    nws_s4_xi[2], nws_s4_eta[2] = -4.292884147430157 , 2.612617574895511
    nws_s4_xi[3], nws_s4_eta[3] = -4.113795100125332 , 2.1139580751416602
    nws_s4_min = 2.0 * np.pi + np.arcsin((nws_s4_eta[3] - hsc_eta[3]) / hsc_fov)
    nws_s4_max = 2.0 * np.pi - np.arccos((nws_s4_xi[0]  - hsc_xi[3] ) / hsc_fov)
    dtheta = (nws_s4_max - nws_s4_min) / 64
    for ii in range(64):
        nws_s4_xi[4 + ii]  = hsc_xi[3]  + hsc_fov * np.cos(nws_s4_min + float(ii) * dtheta)
        nws_s4_eta[4 + ii] = hsc_eta[3] + hsc_fov * np.sin(nws_s4_min + float(ii) * dtheta)

    nws_s1_xi[Nnws_s1 - 1] = nws_s1_xi[0];    nws_s1_eta[Nnws_s1 - 1] = nws_s1_eta[0]
    nws_s2_xi[Nnws_s2 - 1] = nws_s2_xi[0];    nws_s2_eta[Nnws_s2 - 1] = nws_s2_eta[0]
    nws_s3_xi[Nnws_s3 - 1] = nws_s3_xi[0];    nws_s3_eta[Nnws_s3 - 1] = nws_s3_eta[0]
    nws_s4_xi[Nnws_s4 - 1] = nws_s4_xi[0];    nws_s4_eta[Nnws_s4 - 1] = nws_s4_eta[0]

    return Nhsc, hsc_xi, hsc_eta, hsc_rad, Nnws, nws_xi, nws_eta, nws_D, nws_Dep, nws_Dem, nws_s1_xi, nws_s1_eta, nws_s2_xi, nws_s2_eta, nws_s3_xi, nws_s3_eta, nws_s4_xi, nws_s4_eta


def NWS_obs_field():
    # Table 1 in Komiyama et al. (2018)
    Nhsc = 5
    hsc_xi  = [0] * Nhsc
    hsc_eta = [0] * Nhsc
    hsc_rad = [0] * Nhsc

    hsc_xi[0] = -4.653744156858361 ;    hsc_eta[0] = 3.635683362752261;    hsc_rad[0] = 0.75# f003
    hsc_xi[1] = -5.7227710294896745;    hsc_eta[1] = 4.47091006650441 ;    hsc_rad[1] = 0.75# f004
    hsc_xi[2] = -5.530826768687439 ;    hsc_eta[2] = 5.816784796173433;    hsc_rad[2] = 0.75# f009
    hsc_xi[3] = -4.842739589150247 ;    hsc_eta[3] = 2.290423176281251;    hsc_rad[3] = 0.75# f022
    hsc_xi[4] = -5.912437203950868 ;    hsc_eta[4] = 3.125050968157858;    hsc_rad[4] = 0.75# f023

    return Nhsc, hsc_xi, hsc_eta, hsc_rad


def NWS_anl_field():
    # Table 4 in Komiyama et al. (2018), details are given by Komiyama-san via private communication
    Nfield = 4
    field_xi  = np.zeros((Nfield, 4 + 1))
    field_eta = np.zeros((Nfield, 4 + 1))
    # Stream 1
    field_xi[0][0], field_eta[0][0] = -5.7055741988006705, 6.546142963104346
    field_xi[0][1], field_eta[0][1] = -5.186841492207897 , 5.101772155015601
    field_xi[0][2], field_eta[0][2] = -5.907773581780884 , 5.101772155015601
    field_xi[0][3], field_eta[0][3] = -6.239719954692477 , 6.061673446044403
    # Stream 2
    field_xi[1][0], field_eta[1][0] = -5.186841492207897 , 5.101772155015601
    field_xi[1][1], field_eta[1][1] = -4.701958106244509 , 3.7516521430870844
    field_xi[1][2], field_eta[1][2] = -5.4408844910550105, 3.7516521430870844
    field_xi[1][3], field_eta[1][3] = -5.907773581780884 , 5.101772155015601
    # Stream 3
    field_xi[2][0], field_eta[2][0] = -4.701958106244509 , 3.7516521430870844
    field_xi[2][1], field_eta[2][1] = -4.292884147430157 , 2.612617574895511
    field_xi[2][2], field_eta[2][2] = -5.046991496988008 , 2.612617574895511
    field_xi[2][3], field_eta[2][3] = -5.4408844910550105, 3.7516521430870844
    # Stream 4
    field_xi[3][0], field_eta[3][0] = -4.292884147430157 , 2.612617574895511
    field_xi[3][1], field_eta[3][1] = -4.113795100125332 , 2.1139580751416602
    field_xi[3][2], field_eta[3][2] = -4.6822224216441475, 1.5578016922010893
    field_xi[3][3], field_eta[3][3] = -5.046991496988008 , 2.612617574895511

    for ii in range(Nfield):
        field_xi[ii][4]  = field_xi[ii][0]
        field_eta[ii][4] = field_eta[ii][0]

    return Nfield, field_xi, field_eta



def NWS_distance():
    # Table 4 in Komiyama et al. (2018)
    Nnws = 4
    nws_xi  = [0] * Nnws
    nws_eta = [0] * Nnws
    nws_D   = [0] * Nnws
    nws_Dep = [0] * Nnws
    nws_Dem = [0] * Nnws
    nws_mod = [0] * Nnws
    nws_ran = [0] * Nnws
    nws_sys = [0] * Nnws

    nws_xi[0] = -5.891761379732909;    nws_eta[0] = 5.934934874532188 ;    nws_mod[0] = 24.64;    nws_ran[0] = 0.196;    nws_sys[0] = 0.057# Stream 1
    nws_xi[1] = -5.388266383737322;    nws_eta[1] = 4.426258888399156 ;    nws_mod[1] = 24.62;    nws_ran[1] = 0.186;    nws_sys[1] = 0.057# Stream 2
    nws_xi[2] = -4.942619286035076;    nws_eta[2] = 3.18150303805077  ;    nws_mod[2] = 24.59;    nws_ran[2] = 0.183;    nws_sys[2] = 0.057# Stream 3
    nws_xi[3] = -4.49814447496522 ;    nws_eta[3] = 1.9926453169853198;    nws_mod[3] = 24.58;    nws_ran[3] = 0.185;    nws_sys[3] = 0.057# Stream 4

    for ii in range(Nnws):
        nws_D[ii] = 10.0 ** (nws_mod[ii] / 5.0 - 2.0)
        nws_Dep[ii] = 10.0 ** ((nws_mod[ii] + nws_ran[ii] + nws_sys[ii]) / 5.0 - 2.0) - nws_D[ii]
        nws_Dem[ii] = nws_D[ii] - 10.0 ** ((nws_mod[ii] - nws_ran[ii] - nws_sys[ii]) / 5.0 - 2.0)
        # print(nws_xi[ii], nws_eta[ii], nws_D[ii], nws_Dep[ii], nws_Dem[ii])

    return Nnws, nws_xi, nws_eta, nws_D, nws_Dep, nws_Dem


def NWS_velocity():
    # Table 4 in Veljanoski et al. (2014)
    Ngcs = 5
    nws_gc_xi  = [0] * Ngcs
    nws_gc_eta = [0] * Ngcs
    nws_gc_vel = [0] * Ngcs
    nws_gc_vep = [0] * Ngcs
    nws_gc_vem = [0] * Ngcs

    nws_gc_xi[0] = -6.436733608263167;    nws_gc_eta[0] = 6.4598886185030535;    nws_gc_vel[0] = -397;    nws_gc_vep[0] =  7;    nws_gc_vem[0] =  7;# PAndAS-04
    nws_gc_xi[1] = -5.260959997771581;    nws_gc_eta[1] = 4.061434922115992 ;    nws_gc_vel[1] = -444;    nws_gc_vep[1] = 21;    nws_gc_vem[1] = 21;# PAndAS-09
    nws_gc_xi[2] = -5.124223955350093;    nws_gc_eta[2] = 4.137743013677189 ;    nws_gc_vel[2] = -435;    nws_gc_vep[2] = 10;    nws_gc_vem[2] = 10;# PAndAS-10
    nws_gc_xi[3] = -4.946180555295868;    nws_gc_eta[3] = 3.5546551505544537;    nws_gc_vel[3] = -447;    nws_gc_vep[3] = 13;    nws_gc_vem[3] = 13;# PAndAS-11
    nws_gc_xi[4] = -4.557799679219338;    nws_gc_eta[4] = 2.2086054244275655;    nws_gc_vel[4] = -472;    nws_gc_vep[4] =  5;    nws_gc_vem[4] =  5;# PAndAS-12

    return Ngcs, nws_gc_xi, nws_gc_eta, nws_gc_vel, nws_gc_vep, nws_gc_vem
