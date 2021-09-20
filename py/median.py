import numpy as np

logdir = "log/"
series = "cb17"
append = "magi.breakdown.txt"

output = logdir + series + ".median.txt"
fout = open(output, "w")

for Ntot in [131072, 262144, 524288, 1048576, 2097152, 4194304, 8388608, 16777216, 33554432, 67108864, 134217728, 268435456, 536870912]:
    filename = logdir + series + "." + "N" + str(Ntot) + "." + append
    measured = np.genfromtxt(filename,  comments = '#', delimiter = '\t')# for tsv

    ttot = measured[:, 0]
    bodyAlloc = measured[:, 1]
    fileWrite = measured[:, 2]
    spheAlloc = measured[:, 3]
    spheProf = measured[:, 4]
    spheDist = measured[:, 5]
    spheInfo = measured[:, 6]
    eddington = measured[:, 7]
    diskAlloc = measured[:, 8]
    diskProf = measured[:, 9]
    diskDist = measured[:, 10]
    diskInfo = measured[:, 11]
    diskTbl = measured[:, 12]
    diskVel = measured[:, 13]
    observe = measured[:, 14]
    vdisp   = measured[:, 15]
    column  = measured[:, 16]

    ttot_median = np.median(ttot)
    bodyAlloc_median = np.median(bodyAlloc)
    fileWrite_median = np.median(fileWrite)
    spheAlloc_median = np.median(spheAlloc)
    spheProf_median = np.median(spheProf)
    spheDist_median = np.median(spheDist)
    spheInfo_median = np.median(spheInfo)
    eddington_median = np.median(eddington)
    diskAlloc_median = np.median(diskAlloc)
    diskProf_median = np.median(diskProf)
    diskDist_median = np.median(diskDist)
    diskInfo_median = np.median(diskInfo)
    diskTbl_median = np.median(diskTbl)
    diskVel_median = np.median(diskVel)
    observe_median = np.median(observe)
    vdisp_median = np.median(vdisp)
    column_median = np.median(column)

    print("%d\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e" % (Ntot, ttot_median, bodyAlloc_median, fileWrite_median, spheAlloc_median, spheProf_median, spheDist_median, spheInfo_median, eddington_median, diskAlloc_median, diskProf_median, diskDist_median, diskInfo_median, diskTbl_median, diskVel_median, observe_median, vdisp_median, column_median), file = fout)

fout.close()
