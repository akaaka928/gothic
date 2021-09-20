import numpy as np
import glob

# Ndir = 1
# pos = ["./"]
Ndir = 3
pos = ["original/", "cuda/", "compile/"]

Ndat = 15
num = [1024, 2048, 4096, 8192, 16384, 32768, 65536, 131072, 262144, 524288, 1048576, 2097152, 4194304, 8388608, 16777216]
tag = ["N001k", "N002k", "N004k", "N008k", "N016k", "N032k", "N064k", "N128k", "N256k", "N512k", "N001M", "N002M", "N004M", "N008M", "N016M"]


for hh in range(Ndir):
    output = pos[hh] + "summary.txt"
    fout = open(output, "w")
    for ii in range(Ndat):
        files = sorted(glob.glob(pos[hh] + tag[ii] + "/bench/hernquist.time????????.mean.dat"))

        init = np.genfromtxt(files[0], comments = "#", delimiter = "\t")
        kind = init.shape[init.ndim - 1]

        avg = [0] * kind

        for filename in files:
            tmp = np.genfromtxt(filename, comments = "#", delimiter = "\t")
            avg += tmp

        avg /= len(files)
        print("%d\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e" % (num[ii], avg[0], avg[1], avg[2], avg[3], avg[4], avg[5], avg[6], avg[7], avg[8], avg[9], avg[10], avg[11], avg[12]), file = fout)
