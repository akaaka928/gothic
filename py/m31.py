import numpy as np

zm31 = 776.0# same with Komiyama et al. (2018); median of NED

vm31x = 0.0
vm31y = 0.0
vm31z = -300.0


# rot: observed frame (x, y, z) ==>> disk orthogonal frame (X, Y, Z)
# inv: disk orthogonal frame (X, Y, Z) ==>> observed frame (x, y, z)
def Euler_angle():
    inc = (180.0 - 77.0) * np.pi / 180.0
    pa  =          37.0  * np.pi / 180.0

    sini, sinp = np.sin(inc), np.sin(pa)
    cosi, cosp = np.cos(inc), np.cos(pa)

    rot = np.array([[sinp, cosp, 0.0], [-cosi * cosp, cosi * sinp, sini], [sini * cosp, -sini * sinp, cosi]])
    inv = np.linalg.inv(rot)

    return rot, inv


def standard_coordinate(xx, yy, zz):
    xi = xx / (zm31 + zz)
    eta = yy / (zm31 + zz)
    D = np.sqrt(xx ** 2 + yy ** 2 + (zm31 + zz) ** 2)

    return xi, eta, D


def observed_velocity(xx, yy, zz, vx, vy, vz):
    xi, eta, D = standard_coordinate(xx, yy, zz)

    vlos = (xx * (vx + vm31x) + yy * (vy + vm31y) + (zm31 + zz) * (vz + vm31z)) / D
    vxi = ((zm31 + zz) * vx - xx * vz) / D
    veta = ((zm31 + zz) * vy - yy * vz) / D

    return vxi, veta, vlos


def disk_ellipse():
    num = 128
    xi  = [0] * num
    eta = [0] * num
    D   = [0] * num

    rot, inv = Euler_angle()

    rad = 30.0
    dtheta = 2.0 * np.pi / num
    for ii in range(num):
        theta = dtheta * ii
        pos = np.array([rad * np.cos(theta), rad * np.sin(theta), 0.0])
        obs = np.dot(rot, pos)

        xi[ii] = obs[0]
        eta[ii] = obs[1]
        D[ii] = obs[2]

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
