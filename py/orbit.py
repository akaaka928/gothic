import numpy as np


root = 'cfg/'
series_dir = 'nws/pickup/'

ini_sh = 'sh/slurm/submit_backward-ini.sh'
run_sh = 'sh/slurm/submit_backward-run.sh'


def write_config_file(folder, series, Nkind, test, mass, subhalo, orbit, backward, pos, vel):
    filename = series
    if backward:
        filename = series + '-backward'
    cfg = folder + filename + '.cfg'
    ini = folder + filename + '.ini'

    if backward:
        with open(ini_sh, mode = 'a') as fp:
            fp.write('sbatch --dependency=singleton' + ' ' + '--export=' + 'FILE=' + filename + ',' +  'CFG=' + cfg + ' ' + 'sh/slurm/backward-ini.sh\n')
        with open(run_sh, mode = 'a') as fp:
            fp.write('sbatch --dependency=singleton' + ' ' + '--export=' + 'FILE=' + filename                       + ' ' + 'sh/slurm/backward-run.sh\n')


    with open(root + cfg, mode = 'w') as fp:
        fp.write(str(1) + '\n') # unit system
        fp.write(str(Nkind) + '\n') # number of component(s)
        if Nkind > 1:
            fp.write('k18nws nws/continue.param\n')
        if test:
            fp.write('throw_BH_particle' + ' ' + ini + '\n')
        else:
            ini = folder + orbit + '.ini'
            fp.write(subhalo + ' ' + ini + '\n')

    with open(root + ini, mode = 'w') as fp:
        fp.write(str(0) + '\n')
        fp.write(str(1.0) + ' ' + str(0.0) + ' ' + str(0.0) + '\n')
        fp.write(str(0.0) + ' ' + str(1.0) + ' ' + str(0.0) + '\n')
        fp.write(str(0.0) + ' ' + str(0.0) + ' ' + str(1.0) + '\n')
        if backward:
            fp.write(str(pos[0]) + ' ' + str(pos[1]) + ' ' + str(pos[2]) + '\n')
            fp.write(str(vel[0]) + ' ' + str(vel[1]) + ' ' + str(vel[2]) + '\n')
        else:
            fp.write('\n\n')
        if test:
            fp.write(str(-1) + '\n')
            fp.write(str(mass) + '\n')
        else:
            fp.write(str(0) + '\n')


def prepare_config_files(series):
    for power in [7.0, 7.5, 8.0, 8.5, 9.0, 9.5]:
        mass = np.power(10.0, power)
        filename = series + '_mass' + str(int(power * 10.0))
        write_config_file(series_dir + 'test/', filename + '-test', 2, True, mass, None, None, False, None, None)

        write_config_file(series_dir + 'live/', filename + '-prada'   , 2, False, None, 'prada_mass'    + str(int(power * 10.0)), series, False, None, None)
        write_config_file(series_dir + 'live/', filename + '-ishiyama', 2, False, None, 'ishiyama_mass' + str(int(power * 10.0)), series, False, None, None)
        write_config_file(series_dir + 'live/', filename + '-gilman'  , 2, False, None, 'gilman_mass'   + str(int(power * 10.0)), series, False, None, None)


def generate_config_files(pos, vel, idx, vnorm, tag):
    series = 'snap' + str(idx) + '_vel' + str(int(vnorm)) + '_orbit' + str(tag)

    write_config_file(series_dir + 'back/', series, 1, True, 1.0e+9, None, None, True, pos, vel)
    prepare_config_files(series)


# read summary of the orbit-plane
summary = np.genfromtxt('dat/k18nws-chosen_group.csv', delimiter = ',')
time = summary[:, 1] * 0.1
X0 = summary[:, 2]
Y0 = summary[:, 3]
Z0 = summary[:, 4]
aa = summary[:, 5]
bb = summary[:, 6]
cc = summary[:, 7]

# normalize the vector (a, b, c)
inv = 1.0 / np.sqrt(aa ** 2 + bb ** 2 + cc ** 2)
aa *= inv
bb *= inv
cc *= inv

num = len(time)
# print(num)
idx = num

with open(ini_sh, mode = 'w') as fp:
    fp.write('#!/bin/sh\n')
with open(run_sh, mode = 'w') as fp:
    fp.write('#!/bin/sh\n')


while idx >= 0:
    idx = int(input("input the index of the target snapshot ([0:Nfile - 1]): "))
    # print(idx)

    if idx >= num:
        print("index must be in [0:", num - 1, "], retry")
    else:
        if idx >= 0:
            # determine the collision point
            pos = np.array([X0[idx], Y0[idx], Z0[idx]])
            # print(pos)

            # determine the normal vector
            nml = np.array([aa[idx], bb[idx], cc[idx]])
            # print(nml)

            # calculate
            ncr = np.cross(nml, pos) # vector product (\vb*{n} \times \vb*{r})
            rdn = np.dot(pos, nml)   # scalar product (\vb*{r} \cdot \vb*{n})
            # print(ncr)
            # print(rdn)

            vel = float(input("input the velocity of the DM sub-halo: "))
            # print(vel)

            # orbit 0
            v0 = vel * ncr / np.sqrt(np.dot(ncr, ncr))

            # orbit 1
            v1 = -v0

            # orbit 2
            vec = pos - rdn * nml
            v2 = vel * (vec) / np.sqrt(np.dot(vec, vec))

            # orbit 3
            v3 = -v2

            # print("v0:", v0)
            # print("v1:", v1)
            # print("v2:", v2)
            # print("v3:", v3)
            generate_config_files(pos, v0, idx, vel, 0)
            generate_config_files(pos, v1, idx, vel, 1)
            generate_config_files(pos, v2, idx, vel, 2)
            generate_config_files(pos, v3, idx, vel, 3)
