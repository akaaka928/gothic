import numpy as np

# zm31 = 776.0# same with Komiyama et al. (2018); median of NED
zm31 = 773.0# M31 distance in Conn et al. (2016)

vm31x = 0.0
vm31y = 0.0
vm31z = -300.0

def distance():
    return zm31

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

def standard_coordinate(xx, yy, zz):
    tmp = 180.0 / (np.pi * (zm31 + zz))
    xi = xx * tmp
    eta = yy * tmp
    D = np.sqrt(xx ** 2 + yy ** 2 + (zm31 + zz) ** 2)

    return xi, eta, D


def observed_velocity(xx, yy, zz, vx, vy, vz):
    xi, eta, D = standard_coordinate(xx, yy, zz)

    vlos = (xx * (vx + vm31x) + yy * (vy + vm31y) + (zm31 + zz) * (vz + vm31z)) / D
    vxi = ((zm31 + zz) * vx - xx * vz) / D
    veta = ((zm31 + zz) * vy - yy * vz) / D

    return vxi, veta, vlos


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
        D[ii] -= zm31

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


def GSS_obs_field():
    # Table 1 in Font et al. (2006)
    Nfield = 11
    field_xi  = [0] * Nfield
    field_eta = [0] * Nfield

    field_xi[ 0] =  1.077;    field_eta[ 0] = -2.021# a3
    field_xi[ 1] =  2.015;    field_eta[ 1] = -3.965# 1
    field_xi[ 2] =  1.745;    field_eta[ 2] = -3.525# 2
    field_xi[ 3] =  1.483;    field_eta[ 3] = -3.087# 3
    field_xi[ 4] =  1.226;    field_eta[ 4] = -2.653# 4
    field_xi[ 5] =  0.969;    field_eta[ 5] = -2.264# 5
    field_xi[ 6] =  0.717;    field_eta[ 6] = -1.768# 6
    field_xi[ 7] =  0.467;    field_eta[ 7] = -1.327# 7
    field_xi[ 8] =  0.219;    field_eta[ 8] = -0.886# 8
    field_xi[ 9] = -0.731;    field_eta[ 9] =  0.891# 12
    field_xi[10] = -0.963;    field_eta[10] =  1.342# 13

    return Nfield, field_xi, field_eta


def GSS_distance():
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

    for ii in range(Ngss):
        gss_D[ii] -= zm31

    # NOTE
    # subfields GSSx.5 are skipped
    # GSS1 is skipped to remove overlap with M31 disk

    return Ngss, gss_xi, gss_eta, gss_D, gss_Dep, gss_Dem

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

    for ii in range(NstreamC):
        streamC_D[ii] -= zm31

    return NstreamC, streamC_xi, streamC_eta, streamC_D, streamC_Dep, streamC_Dem


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

    for ii in range(NstreamD):
        streamD_D[ii] -= zm31

    return NstreamD, streamD_xi, streamD_eta, streamD_D, streamD_Dep, streamD_Dem
