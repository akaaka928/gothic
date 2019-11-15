#!/bin/sh
#$ -cwd
#$ -l q_core=1
#$ -l h_rt=0:05:00
#$ -N editor
#$ -hold_jid magi
###############################################################
# NCORES_PER_NODE=28 # for f_node
# NCORES_PER_NODE=14 # for h_node
# NCORES_PER_NODE=7 # for q_node
NCORES_PER_NODE=4 # for q_core
###############################################################
export OMP_NUM_THREADS=$NCORES_PER_NODE
###############################################################


###############################################################
# global configurations
###############################################################
EXEC=bin/editor
###############################################################
# dump file generation interval (in units of minute)
if [ -z "$SAVE" ]; then
    SAVE=55.0
fi
###############################################################
if [ -z "$EPS" ]; then
    EPS=1.5625e-2
fi
if [ -z "$ETA" ]; then
    ETA=0.5
fi
if [ -z "$FINISH" ]; then
    FINISH=14000.0
fi
if [ -z "$INTERVAL" ]; then
    INTERVAL=25.0
fi
###############################################################
# problem ID
if [ -z "$PROBLEM" ]; then
    # reproduction of Komiyama et al. (2018)
    # PROBLEM=13

    # prepare for throwing DM sub-halo into NW stream
    PROBLEM=30
    # PROBLEM=31
    # PROBLEM=32
    # PROBLEM=33
    # PROBLEM=34
    # PROBLEM=35
    # PROBLEM=36

    # continue N-body simulation to compare results with DM sub-halo free simulation
    # PROBLEM=40

    # collision between DM sub-halo (test-particle) and NW stream (t = 4.5 Gyr)
    # PROBLEM=100
    # PROBLEM=101
    # PROBLEM=102
    # PROBLEM=103
    # PROBLEM=104
    # PROBLEM=105
    # PROBLEM=106
    # PROBLEM=107

    # PROBLEM=110
    # PROBLEM=111
    # PROBLEM=112
    # PROBLEM=113
    # PROBLEM=114
    # PROBLEM=115
    # PROBLEM=116
    # PROBLEM=117


    # collision between DM sub-halo (Prada relation) and NW stream (t = 4.5 Gyr)
    # PROBLEM=120
    # PROBLEM=121
    # PROBLEM=122
    # PROBLEM=123
    # PROBLEM=124
    # PROBLEM=125
    # PROBLEM=126
    # PROBLEM=127

    # PROBLEM=130
    # PROBLEM=131
    # PROBLEM=132
    # PROBLEM=133
    # PROBLEM=134
    # PROBLEM=135
    # PROBLEM=136
    # PROBLEM=137


    # collision between DM sub-halo (Ishiyama & Ando relation) and NW stream (t = 4.5 Gyr)
    # PROBLEM=140
    # PROBLEM=141
    # PROBLEM=142
    # PROBLEM=143
    # PROBLEM=144
    # PROBLEM=145
    # PROBLEM=146
    # PROBLEM=147

    # PROBLEM=150
    # PROBLEM=151
    # PROBLEM=152
    # PROBLEM=153
    # PROBLEM=154
    # PROBLEM=155
    # PROBLEM=156
    # PROBLEM=157


    # collision between DM sub-halo (Gilman relation) and NW stream (t = 4.5 Gyr)
    # PROBLEM=160
    # PROBLEM=161
    # PROBLEM=162
    # PROBLEM=163
    # PROBLEM=164
    # PROBLEM=165
    # PROBLEM=166
    # PROBLEM=167

    # PROBLEM=170
    # PROBLEM=171
    # PROBLEM=172
    # PROBLEM=173
    # PROBLEM=174
    # PROBLEM=175
    # PROBLEM=176
    # PROBLEM=177


    # collision between DM sub-halo (test-particle) and NW stream (t = 5.0 Gyr)
    # PROBLEM=200
    # PROBLEM=201
    # PROBLEM=202
    # PROBLEM=203
    # PROBLEM=204
    # PROBLEM=205
    # PROBLEM=206
    # PROBLEM=207

    # PROBLEM=210
    # PROBLEM=211
    # PROBLEM=212
    # PROBLEM=213
    # PROBLEM=214
    # PROBLEM=215
    # PROBLEM=216
    # PROBLEM=217

    # collision between DM sub-halo (Prada relation) and NW stream (t = 5.0 Gyr)
    # PROBLEM=220
    # PROBLEM=221
    # PROBLEM=222
    # PROBLEM=223
    # PROBLEM=224
    # PROBLEM=225
    # PROBLEM=226
    # PROBLEM=227

    # PROBLEM=230
    # PROBLEM=231
    # PROBLEM=232
    # PROBLEM=233
    # PROBLEM=234
    # PROBLEM=235
    # PROBLEM=236
    # PROBLEM=237


    # collision between DM sub-halo (Ishiyama & Ando relation) and NW stream (t = 5.0 Gyr)
    # PROBLEM=240
    # PROBLEM=241
    # PROBLEM=242
    # PROBLEM=243
    # PROBLEM=244
    # PROBLEM=245
    # PROBLEM=246
    # PROBLEM=247

    # PROBLEM=250
    # PROBLEM=251
    # PROBLEM=252
    # PROBLEM=253
    # PROBLEM=254
    # PROBLEM=255
    # PROBLEM=256
    # PROBLEM=257


    # collision between DM sub-halo (Gilman relation) and NW stream (t = 5.0 Gyr)
    # PROBLEM=260
    # PROBLEM=261
    # PROBLEM=262
    # PROBLEM=263
    # PROBLEM=264
    # PROBLEM=265
    # PROBLEM=266
    # PROBLEM=267

    # PROBLEM=270
    # PROBLEM=271
    # PROBLEM=272
    # PROBLEM=273
    # PROBLEM=274
    # PROBLEM=275
    # PROBLEM=276
    # PROBLEM=277


    # collision between DM sub-halo (test-particle) and NW stream (t = 5.5 Gyr)
    # PROBLEM=300
    # PROBLEM=301
    # PROBLEM=302
    # PROBLEM=303
    # PROBLEM=304
    # PROBLEM=305
    # PROBLEM=306
    # PROBLEM=307

    # PROBLEM=310
    # PROBLEM=311
    # PROBLEM=312
    # PROBLEM=313
    # PROBLEM=314
    # PROBLEM=315
    # PROBLEM=316
    # PROBLEM=317


    # collision between DM sub-halo (Prada relation) and NW stream (t = 5.5 Gyr)
    # PROBLEM=320
    # PROBLEM=321
    # PROBLEM=322
    # PROBLEM=323
    # PROBLEM=324
    # PROBLEM=325
    # PROBLEM=326
    # PROBLEM=327

    # PROBLEM=330
    # PROBLEM=331
    # PROBLEM=332
    # PROBLEM=333
    # PROBLEM=334
    # PROBLEM=335
    # PROBLEM=336
    # PROBLEM=337

    # collision between DM sub-halo (Ishiyama & Ando relation) and NW stream (t = 5.5 Gyr)
    # PROBLEM=340
    # PROBLEM=341
    # PROBLEM=342
    # PROBLEM=343
    # PROBLEM=344
    # PROBLEM=345
    # PROBLEM=346
    # PROBLEM=347

    # PROBLEM=350
    # PROBLEM=351
    # PROBLEM=352
    # PROBLEM=353
    # PROBLEM=354
    # PROBLEM=355
    # PROBLEM=356
    # PROBLEM=357


    # collision between DM sub-halo (Gilman relation) and NW stream (t = 5.5 Gyr)
    # PROBLEM=360
    # PROBLEM=361
    # PROBLEM=362
    # PROBLEM=363
    # PROBLEM=364
    # PROBLEM=365
    # PROBLEM=366
    # PROBLEM=367

    # PROBLEM=370
    # PROBLEM=371
    # PROBLEM=372
    # PROBLEM=373
    # PROBLEM=374
    # PROBLEM=375
    # PROBLEM=376
    # PROBLEM=377


    # collision between DM sub-halo (test-particle) and NW stream (t = 6.0 Gyr)
    # PROBLEM=400
    # PROBLEM=401
    # PROBLEM=402
    # PROBLEM=403
    # PROBLEM=404
    # PROBLEM=405
    # PROBLEM=406
    # PROBLEM=407

    # PROBLEM=410
    # PROBLEM=411
    # PROBLEM=412
    # PROBLEM=413
    # PROBLEM=414
    # PROBLEM=415
    # PROBLEM=416
    # PROBLEM=417

    # collision between DM sub-halo (Prada relation) and NW stream (t = 6.0 Gyr)
    # PROBLEM=420
    # PROBLEM=421
    # PROBLEM=422
    # PROBLEM=423
    # PROBLEM=424
    # PROBLEM=425
    # PROBLEM=426
    # PROBLEM=427

    # PROBLEM=430
    # PROBLEM=431
    # PROBLEM=432
    # PROBLEM=433
    # PROBLEM=434
    # PROBLEM=435
    # PROBLEM=436
    # PROBLEM=437


    # collision between DM sub-halo (Ishiyama & Ando relation) and NW stream (t = 6.0 Gyr)
    # PROBLEM=440
    # PROBLEM=441
    # PROBLEM=442
    # PROBLEM=443
    # PROBLEM=444
    # PROBLEM=445
    # PROBLEM=446
    # PROBLEM=447

    # PROBLEM=450
    # PROBLEM=451
    # PROBLEM=452
    # PROBLEM=453
    # PROBLEM=454
    # PROBLEM=455
    # PROBLEM=456
    # PROBLEM=457

    # collision between DM sub-halo (Gilman relation) and NW stream (t = 6.0 Gyr)
    # PROBLEM=460
    # PROBLEM=461
    # PROBLEM=462
    # PROBLEM=463
    # PROBLEM=464
    # PROBLEM=465
    # PROBLEM=466
    # PROBLEM=467

    # PROBLEM=470
    # PROBLEM=471
    # PROBLEM=472
    # PROBLEM=473
    # PROBLEM=474
    # PROBLEM=475
    # PROBLEM=476
    # PROBLEM=477


    # collision between DM sub-halo (test-particle) and NW stream (t = 6.5 Gyr)
    # PROBLEM=500
    # PROBLEM=501
    # PROBLEM=502
    # PROBLEM=503
    # PROBLEM=504
    # PROBLEM=505
    # PROBLEM=506
    # PROBLEM=507

    # PROBLEM=510
    # PROBLEM=511
    # PROBLEM=512
    # PROBLEM=513
    # PROBLEM=514
    # PROBLEM=515
    # PROBLEM=516
    # PROBLEM=517


    # collision between DM sub-halo (Prada relation) and NW stream (t = 6.5 Gyr)
    # PROBLEM=520
    # PROBLEM=521
    # PROBLEM=522
    # PROBLEM=523
    # PROBLEM=524
    # PROBLEM=525
    # PROBLEM=526
    # PROBLEM=527

    # PROBLEM=530
    # PROBLEM=531
    # PROBLEM=532
    # PROBLEM=533
    # PROBLEM=534
    # PROBLEM=535
    # PROBLEM=536
    # PROBLEM=537


    # collision between DM sub-halo (Ishiyama & Ando relation) and NW stream (t = 6.5 Gyr)
    # PROBLEM=540
    # PROBLEM=541
    # PROBLEM=542
    # PROBLEM=543
    # PROBLEM=544
    # PROBLEM=545
    # PROBLEM=546
    # PROBLEM=547

    # PROBLEM=550
    # PROBLEM=551
    # PROBLEM=552
    # PROBLEM=553
    # PROBLEM=554
    # PROBLEM=555
    # PROBLEM=556
    # PROBLEM=557


    # collision between DM sub-halo (Gilman relation) and NW stream (t = 6.5 Gyr)
    # PROBLEM=560
    # PROBLEM=561
    # PROBLEM=562
    # PROBLEM=563
    # PROBLEM=564
    # PROBLEM=565
    # PROBLEM=566
    # PROBLEM=567

    # PROBLEM=570
    # PROBLEM=571
    # PROBLEM=572
    # PROBLEM=573
    # PROBLEM=574
    # PROBLEM=575
    # PROBLEM=576
    # PROBLEM=577


    # collision between DM sub-halo (test-particle) and NW stream (t = 7.0 Gyr)
    # PROBLEM=600
    # PROBLEM=601
    # PROBLEM=602
    # PROBLEM=603
    # PROBLEM=604
    # PROBLEM=605
    # PROBLEM=606
    # PROBLEM=607

    # PROBLEM=610
    # PROBLEM=611
    # PROBLEM=612
    # PROBLEM=613
    # PROBLEM=614
    # PROBLEM=615
    # PROBLEM=616
    # PROBLEM=617


    # collision between DM sub-halo (Prada relation) and NW stream (t = 7.0 Gyr)
    # PROBLEM=620
    # PROBLEM=621
    # PROBLEM=622
    # PROBLEM=623
    # PROBLEM=624
    # PROBLEM=625
    # PROBLEM=626
    # PROBLEM=627

    # PROBLEM=630
    # PROBLEM=631
    # PROBLEM=632
    # PROBLEM=633
    # PROBLEM=634
    # PROBLEM=635
    # PROBLEM=636
    # PROBLEM=637


    # collision between DM sub-halo (Ishiyama & Ando relation) and NW stream (t = 7.0 Gyr)
    # PROBLEM=640
    # PROBLEM=641
    # PROBLEM=642
    # PROBLEM=643
    # PROBLEM=644
    # PROBLEM=645
    # PROBLEM=646
    # PROBLEM=647

    # PROBLEM=650
    # PROBLEM=651
    # PROBLEM=652
    # PROBLEM=653
    # PROBLEM=654
    # PROBLEM=655
    # PROBLEM=656
    # PROBLEM=657


    # collision between DM sub-halo (Gilman relation) and NW stream (t = 7.0 Gyr)
    # PROBLEM=660
    # PROBLEM=661
    # PROBLEM=662
    # PROBLEM=663
    # PROBLEM=664
    # PROBLEM=665
    # PROBLEM=666
    # PROBLEM=667

    # PROBLEM=670
    # PROBLEM=671
    # PROBLEM=672
    # PROBLEM=673
    # PROBLEM=674
    # PROBLEM=675
    # PROBLEM=676
    # PROBLEM=677


    # collision between DM sub-halo (test-particle) and NW stream (t = 7.5 Gyr)
    # PROBLEM=700
    # PROBLEM=701
    # PROBLEM=702
    # PROBLEM=703
    # PROBLEM=704
    # PROBLEM=705
    # PROBLEM=706
    # PROBLEM=707

    # PROBLEM=710
    # PROBLEM=711
    # PROBLEM=712
    # PROBLEM=713
    # PROBLEM=714
    # PROBLEM=715
    # PROBLEM=716
    # PROBLEM=717


    # collision between DM sub-halo (Prada relation) and NW stream (t = 7.5 Gyr)
    # PROBLEM=720
    # PROBLEM=721
    # PROBLEM=722
    # PROBLEM=723
    # PROBLEM=724
    # PROBLEM=725
    # PROBLEM=726
    # PROBLEM=727

    # PROBLEM=730
    # PROBLEM=731
    # PROBLEM=732
    # PROBLEM=733
    # PROBLEM=734
    # PROBLEM=735
    # PROBLEM=736
    # PROBLEM=737


    # collision between DM sub-halo (Ishiyama & Ando relation) and NW stream (t = 7.5 Gyr)
    # PROBLEM=740
    # PROBLEM=741
    # PROBLEM=742
    # PROBLEM=743
    # PROBLEM=744
    # PROBLEM=745
    # PROBLEM=746
    # PROBLEM=747

    # PROBLEM=750
    # PROBLEM=751
    # PROBLEM=752
    # PROBLEM=753
    # PROBLEM=754
    # PROBLEM=755
    # PROBLEM=756
    # PROBLEM=757


    # collision between DM sub-halo (Gilman relation) and NW stream (t = 7.5 Gyr)
    # PROBLEM=760
    # PROBLEM=761
    # PROBLEM=762
    # PROBLEM=763
    # PROBLEM=764
    # PROBLEM=765
    # PROBLEM=766
    # PROBLEM=767

    # PROBLEM=770
    # PROBLEM=771
    # PROBLEM=772
    # PROBLEM=773
    # PROBLEM=774
    # PROBLEM=775
    # PROBLEM=776
    # PROBLEM=777


    # collision between DM sub-halo (test-particle) and NW stream (t = 8.0 Gyr)
    # PROBLEM=800
    # PROBLEM=801
    # PROBLEM=802
    # PROBLEM=803
    # PROBLEM=804
    # PROBLEM=805
    # PROBLEM=806
    # PROBLEM=807

    # PROBLEM=810
    # PROBLEM=811
    # PROBLEM=812
    # PROBLEM=813
    # PROBLEM=814
    # PROBLEM=815
    # PROBLEM=816
    # PROBLEM=817


    # collision between DM sub-halo (Prada relation) and NW stream (t = 8.0 Gyr)
    # PROBLEM=820
    # PROBLEM=821
    # PROBLEM=822
    # PROBLEM=823
    # PROBLEM=824
    # PROBLEM=825
    # PROBLEM=826
    # PROBLEM=827

    # PROBLEM=830
    # PROBLEM=831
    # PROBLEM=832
    # PROBLEM=833
    # PROBLEM=834
    # PROBLEM=835
    # PROBLEM=836
    # PROBLEM=837


    # collision between DM sub-halo (Ishiyama & Ando relation) and NW stream (t = 8.0 Gyr)
    # PROBLEM=840
    # PROBLEM=841
    # PROBLEM=842
    # PROBLEM=843
    # PROBLEM=844
    # PROBLEM=845
    # PROBLEM=846
    # PROBLEM=847

    # PROBLEM=850
    # PROBLEM=851
    # PROBLEM=852
    # PROBLEM=853
    # PROBLEM=854
    # PROBLEM=855
    # PROBLEM=856
    # PROBLEM=857


    # collision between DM sub-halo (Gilman relation) and NW stream (t = 8.0 Gyr)
    # PROBLEM=860
    # PROBLEM=861
    # PROBLEM=862
    # PROBLEM=863
    # PROBLEM=864
    # PROBLEM=865
    # PROBLEM=866
    # PROBLEM=867

    # PROBLEM=870
    # PROBLEM=871
    # PROBLEM=872
    # PROBLEM=873
    # PROBLEM=874
    # PROBLEM=875
    # PROBLEM=876
    # PROBLEM=877
fi
###############################################################


###############################################################
# problem specific configurations
###############################################################
# reproduction of Komiyama et al. (2018) in the disk coordinate system
if [ $PROBLEM -eq 13 ]; then
    FILE=k18nws
    CFG=nws/satellite.cfg
fi
###############################################################
# preparation for throwing DM sub-halo into NW stream
if [ $PROBLEM -eq 30 ]; then
    FILE=k18nws
    CFG=nws/satellite.cfg
    FINISH=4500.0
fi
###############################################################
# preparation for throwing DM sub-halo into NW stream
if [ $PROBLEM -eq 31 ]; then
    FILE=k18nws
    CFG=nws/satellite.cfg
    FINISH=5000.0
fi
###############################################################
# preparation for throwing DM sub-halo into NW stream
if [ $PROBLEM -eq 32 ]; then
    FILE=k18nws
    CFG=nws/satellite.cfg
    FINISH=5500.0
fi
###############################################################
# preparation for throwing DM sub-halo into NW stream
if [ $PROBLEM -eq 33 ]; then
    FILE=k18nws
    CFG=nws/satellite.cfg
    FINISH=6000.0
fi
###############################################################
# preparation for throwing DM sub-halo into NW stream
if [ $PROBLEM -eq 34 ]; then
    FILE=k18nws
    CFG=nws/satellite.cfg
    FINISH=6500.0
fi
###############################################################
# preparation for throwing DM sub-halo into NW stream
if [ $PROBLEM -eq 35 ]; then
    FILE=k18nws
    CFG=nws/satellite.cfg
    FINISH=7000.0
fi
###############################################################
# preparation for throwing DM sub-halo into NW stream
if [ $PROBLEM -eq 36 ]; then
    FILE=k18nws
    CFG=nws/satellite.cfg
    FINISH=7500.0
fi
###############################################################
# continue N-body simulation to reproduce NW stream
if [ $PROBLEM -eq 40 ]; then
    FILE=nws-continue
    CFG=nws/continue.cfg
fi
###############################################################
# collision between DM sub-halo (test-particle) and NW stream (t = 4.5 Gyr)
if [ $PROBLEM -eq 100 ]; then
    FILE=mass80-test-snap160_vel140_orbit0
    CFG=nws/pickup/test/snap160_vel140_orbit0_mass80-test.cfg
fi
###############################################################
# collision between DM sub-halo (test-particle) and NW stream (t = 4.5 Gyr)
if [ $PROBLEM -eq 101 ]; then
    FILE=mass80-test-snap160_vel140_orbit1
    CFG=nws/pickup/test/snap160_vel140_orbit1_mass80-test.cfg
fi
###############################################################
# collision between DM sub-halo (test-particle) and NW stream (t = 4.5 Gyr)
if [ $PROBLEM -eq 102 ]; then
    FILE=mass80-test-snap160_vel140_orbit2
    CFG=nws/pickup/test/snap160_vel140_orbit2_mass80-test.cfg
fi
###############################################################
# collision between DM sub-halo (test-particle) and NW stream (t = 4.5 Gyr)
if [ $PROBLEM -eq 103 ]; then
    FILE=mass80-test-snap160_vel140_orbit3
    CFG=nws/pickup/test/snap160_vel140_orbit3_mass80-test.cfg
fi
###############################################################
# collision between DM sub-halo (test-particle) and NW stream (t = 4.5 Gyr)
if [ $PROBLEM -eq 104 ]; then
    FILE=mass85-test-snap160_vel140_orbit0
    CFG=nws/pickup/test/snap160_vel140_orbit0_mass85-test.cfg
fi
###############################################################
# collision between DM sub-halo (test-particle) and NW stream (t = 4.5 Gyr)
if [ $PROBLEM -eq 105 ]; then
    FILE=mass85-test-snap160_vel140_orbit1
    CFG=nws/pickup/test/snap160_vel140_orbit1_mass85-test.cfg
fi
###############################################################
# collision between DM sub-halo (test-particle) and NW stream (t = 4.5 Gyr)
if [ $PROBLEM -eq 106 ]; then
    FILE=mass85-test-snap160_vel140_orbit2
    CFG=nws/pickup/test/snap160_vel140_orbit2_mass85-test.cfg
fi
###############################################################
# collision between DM sub-halo (test-particle) and NW stream (t = 4.5 Gyr)
if [ $PROBLEM -eq 107 ]; then
    FILE=mass85-test-snap160_vel140_orbit3
    CFG=nws/pickup/test/snap160_vel140_orbit3_mass85-test.cfg
fi
###############################################################
# collision between DM sub-halo (test-particle) and NW stream (t = 4.5 Gyr)
if [ $PROBLEM -eq 110 ]; then
    FILE=mass90-test-snap160_vel140_orbit0
    CFG=nws/pickup/test/snap160_vel140_orbit0_mass90-test.cfg
fi
###############################################################
# collision between DM sub-halo (test-particle) and NW stream (t = 4.5 Gyr)
if [ $PROBLEM -eq 111 ]; then
    FILE=mass90-test-snap160_vel140_orbit1
    CFG=nws/pickup/test/snap160_vel140_orbit1_mass90-test.cfg
fi
###############################################################
# collision between DM sub-halo (test-particle) and NW stream (t = 4.5 Gyr)
if [ $PROBLEM -eq 112 ]; then
    FILE=mass90-test-snap160_vel140_orbit2
    CFG=nws/pickup/test/snap160_vel140_orbit2_mass90-test.cfg
fi
###############################################################
# collision between DM sub-halo (test-particle) and NW stream (t = 4.5 Gyr)
if [ $PROBLEM -eq 113 ]; then
    FILE=mass90-test-snap160_vel140_orbit3
    CFG=nws/pickup/test/snap160_vel140_orbit3_mass90-test.cfg
fi
###############################################################
# collision between DM sub-halo (test-particle) and NW stream (t = 4.5 Gyr)
if [ $PROBLEM -eq 114 ]; then
    FILE=mass95-test-snap160_vel140_orbit0
    CFG=nws/pickup/test/snap160_vel140_orbit0_mass95-test.cfg
fi
###############################################################
# collision between DM sub-halo (test-particle) and NW stream (t = 4.5 Gyr)
if [ $PROBLEM -eq 115 ]; then
    FILE=mass95-test-snap160_vel140_orbit1
    CFG=nws/pickup/test/snap160_vel140_orbit1_mass95-test.cfg
fi
###############################################################
# collision between DM sub-halo (test-particle) and NW stream (t = 4.5 Gyr)
if [ $PROBLEM -eq 116 ]; then
    FILE=mass95-test-snap160_vel140_orbit2
    CFG=nws/pickup/test/snap160_vel140_orbit2_mass95-test.cfg
fi
###############################################################
# collision between DM sub-halo (test-particle) and NW stream (t = 4.5 Gyr)
if [ $PROBLEM -eq 117 ]; then
    FILE=mass95-test-snap160_vel140_orbit3
    CFG=nws/pickup/test/snap160_vel140_orbit3_mass95-test.cfg
fi
###############################################################
# collision between DM sub-halo (Prada relation) and NW stream (t = 4.5 Gyr)
if [ $PROBLEM -eq 120 ]; then
    FILE=mass80-prada-snap160_vel140_orbit0
    CFG=nws/pickup/live/snap160_vel140_orbit0_mass80-prada.cfg
fi
###############################################################
# collision between DM sub-halo (Prada relation) and NW stream (t = 4.5 Gyr)
if [ $PROBLEM -eq 121 ]; then
    FILE=mass80-prada-snap160_vel140_orbit1
    CFG=nws/pickup/live/snap160_vel140_orbit1_mass80-prada.cfg
fi
###############################################################
# collision between DM sub-halo (Prada relation) and NW stream (t = 4.5 Gyr)
if [ $PROBLEM -eq 122 ]; then
    FILE=mass80-prada-snap160_vel140_orbit2
    CFG=nws/pickup/live/snap160_vel140_orbit2_mass80-prada.cfg
fi
###############################################################
# collision between DM sub-halo (Prada relation) and NW stream (t = 4.5 Gyr)
if [ $PROBLEM -eq 123 ]; then
    FILE=mass80-prada-snap160_vel140_orbit3
    CFG=nws/pickup/live/snap160_vel140_orbit3_mass80-prada.cfg
fi
###############################################################
# collision between DM sub-halo (Prada relation) and NW stream (t = 4.5 Gyr)
if [ $PROBLEM -eq 124 ]; then
    FILE=mass85-prada-snap160_vel140_orbit0
    CFG=nws/pickup/live/snap160_vel140_orbit0_mass85-prada.cfg
fi
###############################################################
# collision between DM sub-halo (Prada relation) and NW stream (t = 4.5 Gyr)
if [ $PROBLEM -eq 125 ]; then
    FILE=mass85-prada-snap160_vel140_orbit1
    CFG=nws/pickup/live/snap160_vel140_orbit1_mass85-prada.cfg
fi
###############################################################
# collision between DM sub-halo (Prada relation) and NW stream (t = 4.5 Gyr)
if [ $PROBLEM -eq 126 ]; then
    FILE=mass85-prada-snap160_vel140_orbit2
    CFG=nws/pickup/live/snap160_vel140_orbit2_mass85-prada.cfg
fi
###############################################################
# collision between DM sub-halo (Prada relation) and NW stream (t = 4.5 Gyr)
if [ $PROBLEM -eq 127 ]; then
    FILE=mass85-prada-snap160_vel140_orbit3
    CFG=nws/pickup/live/snap160_vel140_orbit3_mass85-prada.cfg
fi
###############################################################
# collision between DM sub-halo (Prada relation) and NW stream (t = 4.5 Gyr)
if [ $PROBLEM -eq 130 ]; then
    FILE=mass90-prada-snap160_vel140_orbit0
    CFG=nws/pickup/live/snap160_vel140_orbit0_mass90-prada.cfg
fi
###############################################################
# collision between DM sub-halo (Prada relation) and NW stream (t = 4.5 Gyr)
if [ $PROBLEM -eq 131 ]; then
    FILE=mass90-prada-snap160_vel140_orbit1
    CFG=nws/pickup/live/snap160_vel140_orbit1_mass90-prada.cfg
fi
###############################################################
# collision between DM sub-halo (Prada relation) and NW stream (t = 4.5 Gyr)
if [ $PROBLEM -eq 132 ]; then
    FILE=mass90-prada-snap160_vel140_orbit2
    CFG=nws/pickup/live/snap160_vel140_orbit2_mass90-prada.cfg
fi
###############################################################
# collision between DM sub-halo (Prada relation) and NW stream (t = 4.5 Gyr)
if [ $PROBLEM -eq 133 ]; then
    FILE=mass90-prada-snap160_vel140_orbit3
    CFG=nws/pickup/live/snap160_vel140_orbit3_mass90-prada.cfg
fi
###############################################################
# collision between DM sub-halo (Prada relation) and NW stream (t = 4.5 Gyr)
if [ $PROBLEM -eq 134 ]; then
    FILE=mass95-prada-snap160_vel140_orbit0
    CFG=nws/pickup/live/snap160_vel140_orbit0_mass95-prada.cfg
fi
###############################################################
# collision between DM sub-halo (Prada relation) and NW stream (t = 4.5 Gyr)
if [ $PROBLEM -eq 135 ]; then
    FILE=mass95-prada-snap160_vel140_orbit1
    CFG=nws/pickup/live/snap160_vel140_orbit1_mass95-prada.cfg
fi
###############################################################
# collision between DM sub-halo (Prada relation) and NW stream (t = 4.5 Gyr)
if [ $PROBLEM -eq 136 ]; then
    FILE=mass95-prada-snap160_vel140_orbit2
    CFG=nws/pickup/live/snap160_vel140_orbit2_mass95-prada.cfg
fi
###############################################################
# collision between DM sub-halo (Prada relation) and NW stream (t = 4.5 Gyr)
if [ $PROBLEM -eq 137 ]; then
    FILE=mass95-prada-snap160_vel140_orbit3
    CFG=nws/pickup/live/snap160_vel140_orbit3_mass95-prada.cfg
fi
###############################################################
# collision between DM sub-halo (Ishiyama & Ando relation) and NW stream (t = 4.5 Gyr)
if [ $PROBLEM -eq 140 ]; then
    FILE=mass80-ishiyama-snap160_vel140_orbit0
    CFG=nws/pickup/live/snap160_vel140_orbit0_mass80-ishiyama.cfg
fi
###############################################################
# collision between DM sub-halo (Ishiyama & Ando relation) and NW stream (t = 4.5 Gyr)
if [ $PROBLEM -eq 141 ]; then
    FILE=mass80-ishiyama-snap160_vel140_orbit1
    CFG=nws/pickup/live/snap160_vel140_orbit1_mass80-ishiyama.cfg
fi
###############################################################
# collision between DM sub-halo (Ishiyama & Ando relation) and NW stream (t = 4.5 Gyr)
if [ $PROBLEM -eq 142 ]; then
    FILE=mass80-ishiyama-snap160_vel140_orbit2
    CFG=nws/pickup/live/snap160_vel140_orbit2_mass80-ishiyama.cfg
fi
###############################################################
# collision between DM sub-halo (Ishiyama & Ando relation) and NW stream (t = 4.5 Gyr)
if [ $PROBLEM -eq 143 ]; then
    FILE=mass80-ishiyama-snap160_vel140_orbit3
    CFG=nws/pickup/live/snap160_vel140_orbit3_mass80-ishiyama.cfg
fi
###############################################################
# collision between DM sub-halo (Ishiyama & Ando relation) and NW stream (t = 4.5 Gyr)
if [ $PROBLEM -eq 144 ]; then
    FILE=mass85-ishiyama-snap160_vel140_orbit0
    CFG=nws/pickup/live/snap160_vel140_orbit0_mass85-ishiyama.cfg
fi
###############################################################
# collision between DM sub-halo (Ishiyama & Ando relation) and NW stream (t = 4.5 Gyr)
if [ $PROBLEM -eq 145 ]; then
    FILE=mass85-ishiyama-snap160_vel140_orbit1
    CFG=nws/pickup/live/snap160_vel140_orbit1_mass85-ishiyama.cfg
fi
###############################################################
# collision between DM sub-halo (Ishiyama & Ando relation) and NW stream (t = 4.5 Gyr)
if [ $PROBLEM -eq 146 ]; then
    FILE=mass85-ishiyama-snap160_vel140_orbit2
    CFG=nws/pickup/live/snap160_vel140_orbit2_mass85-ishiyama.cfg
fi
###############################################################
# collision between DM sub-halo (Ishiyama & Ando relation) and NW stream (t = 4.5 Gyr)
if [ $PROBLEM -eq 147 ]; then
    FILE=mass85-ishiyama-snap160_vel140_orbit3
    CFG=nws/pickup/live/snap160_vel140_orbit3_mass85-ishiyama.cfg
fi
###############################################################
# collision between DM sub-halo (Ishiyama & Ando relation) and NW stream (t = 4.5 Gyr)
if [ $PROBLEM -eq 150 ]; then
    FILE=mass90-ishiyama-snap160_vel140_orbit0
    CFG=nws/pickup/live/snap160_vel140_orbit0_mass90-ishiyama.cfg
fi
###############################################################
# collision between DM sub-halo (Ishiyama & Ando relation) and NW stream (t = 4.5 Gyr)
if [ $PROBLEM -eq 151 ]; then
    FILE=mass90-ishiyama-snap160_vel140_orbit1
    CFG=nws/pickup/live/snap160_vel140_orbit1_mass90-ishiyama.cfg
fi
###############################################################
# collision between DM sub-halo (Ishiyama & Ando relation) and NW stream (t = 4.5 Gyr)
if [ $PROBLEM -eq 152 ]; then
    FILE=mass90-ishiyama-snap160_vel140_orbit2
    CFG=nws/pickup/live/snap160_vel140_orbit2_mass90-ishiyama.cfg
fi
###############################################################
# collision between DM sub-halo (Ishiyama & Ando relation) and NW stream (t = 4.5 Gyr)
if [ $PROBLEM -eq 153 ]; then
    FILE=mass90-ishiyama-snap160_vel140_orbit3
    CFG=nws/pickup/live/snap160_vel140_orbit3_mass90-ishiyama.cfg
fi
###############################################################
# collision between DM sub-halo (Ishiyama & Ando relation) and NW stream (t = 4.5 Gyr)
if [ $PROBLEM -eq 154 ]; then
    FILE=mass95-ishiyama-snap160_vel140_orbit0
    CFG=nws/pickup/live/snap160_vel140_orbit0_mass95-ishiyama.cfg
fi
###############################################################
# collision between DM sub-halo (Ishiyama & Ando relation) and NW stream (t = 4.5 Gyr)
if [ $PROBLEM -eq 155 ]; then
    FILE=mass95-ishiyama-snap160_vel140_orbit1
    CFG=nws/pickup/live/snap160_vel140_orbit1_mass95-ishiyama.cfg
fi
###############################################################
# collision between DM sub-halo (Ishiyama & Ando relation) and NW stream (t = 4.5 Gyr)
if [ $PROBLEM -eq 156 ]; then
    FILE=mass95-ishiyama-snap160_vel140_orbit2
    CFG=nws/pickup/live/snap160_vel140_orbit2_mass95-ishiyama.cfg
fi
###############################################################
# collision between DM sub-halo (Ishiyama & Ando relation) and NW stream (t = 4.5 Gyr)
if [ $PROBLEM -eq 157 ]; then
    FILE=mass95-ishiyama-snap160_vel140_orbit3
    CFG=nws/pickup/live/snap160_vel140_orbit3_mass95-ishiyama.cfg
fi
###############################################################
# collision between DM sub-halo (Gilman relation) and NW stream (t = 4.5 Gyr)
if [ $PROBLEM -eq 160 ]; then
    FILE=mass80-gilman-snap160_vel140_orbit0
    CFG=nws/pickup/live/snap160_vel140_orbit0_mass80-gilman.cfg
fi
###############################################################
# collision between DM sub-halo (Gilman relation) and NW stream (t = 4.5 Gyr)
if [ $PROBLEM -eq 161 ]; then
    FILE=mass80-gilman-snap160_vel140_orbit1
    CFG=nws/pickup/live/snap160_vel140_orbit1_mass80-gilman.cfg
fi
###############################################################
# collision between DM sub-halo (Gilman relation) and NW stream (t = 4.5 Gyr)
if [ $PROBLEM -eq 162 ]; then
    FILE=mass80-gilman-snap160_vel140_orbit2
    CFG=nws/pickup/live/snap160_vel140_orbit2_mass80-gilman.cfg
fi
###############################################################
# collision between DM sub-halo (Gilman relation) and NW stream (t = 4.5 Gyr)
if [ $PROBLEM -eq 163 ]; then
    FILE=mass80-gilman-snap160_vel140_orbit3
    CFG=nws/pickup/live/snap160_vel140_orbit3_mass80-gilman.cfg
fi
###############################################################
# collision between DM sub-halo (Gilman relation) and NW stream (t = 4.5 Gyr)
if [ $PROBLEM -eq 164 ]; then
    FILE=mass85-gilman-snap160_vel140_orbit0
    CFG=nws/pickup/live/snap160_vel140_orbit0_mass85-gilman.cfg
fi
###############################################################
# collision between DM sub-halo (Gilman relation) and NW stream (t = 4.5 Gyr)
if [ $PROBLEM -eq 165 ]; then
    FILE=mass85-gilman-snap160_vel140_orbit1
    CFG=nws/pickup/live/snap160_vel140_orbit1_mass85-gilman.cfg
fi
###############################################################
# collision between DM sub-halo (Gilman relation) and NW stream (t = 4.5 Gyr)
if [ $PROBLEM -eq 166 ]; then
    FILE=mass85-gilman-snap160_vel140_orbit2
    CFG=nws/pickup/live/snap160_vel140_orbit2_mass85-gilman.cfg
fi
###############################################################
# collision between DM sub-halo (Gilman relation) and NW stream (t = 4.5 Gyr)
if [ $PROBLEM -eq 167 ]; then
    FILE=mass85-gilman-snap160_vel140_orbit3
    CFG=nws/pickup/live/snap160_vel140_orbit3_mass85-gilman.cfg
fi
###############################################################
# collision between DM sub-halo (Gilman relation) and NW stream (t = 4.5 Gyr)
if [ $PROBLEM -eq 170 ]; then
    FILE=mass90-gilman-snap160_vel140_orbit0
    CFG=nws/pickup/live/snap160_vel140_orbit0_mass90-gilman.cfg
fi
###############################################################
# collision between DM sub-halo (Gilman relation) and NW stream (t = 4.5 Gyr)
if [ $PROBLEM -eq 171 ]; then
    FILE=mass90-gilman-snap160_vel140_orbit1
    CFG=nws/pickup/live/snap160_vel140_orbit1_mass90-gilman.cfg
fi
###############################################################
# collision between DM sub-halo (Gilman relation) and NW stream (t = 4.5 Gyr)
if [ $PROBLEM -eq 172 ]; then
    FILE=mass90-gilman-snap160_vel140_orbit2
    CFG=nws/pickup/live/snap160_vel140_orbit2_mass90-gilman.cfg
fi
###############################################################
# collision between DM sub-halo (Gilman relation) and NW stream (t = 4.5 Gyr)
if [ $PROBLEM -eq 173 ]; then
    FILE=mass90-gilman-snap160_vel140_orbit3
    CFG=nws/pickup/live/snap160_vel140_orbit3_mass90-gilman.cfg
fi
###############################################################
# collision between DM sub-halo (Gilman relation) and NW stream (t = 4.5 Gyr)
if [ $PROBLEM -eq 174 ]; then
    FILE=mass95-gilman-snap160_vel140_orbit0
    CFG=nws/pickup/live/snap160_vel140_orbit0_mass95-gilman.cfg
fi
###############################################################
# collision between DM sub-halo (Gilman relation) and NW stream (t = 4.5 Gyr)
if [ $PROBLEM -eq 175 ]; then
    FILE=mass95-gilman-snap160_vel140_orbit1
    CFG=nws/pickup/live/snap160_vel140_orbit1_mass95-gilman.cfg
fi
###############################################################
# collision between DM sub-halo (Gilman relation) and NW stream (t = 4.5 Gyr)
if [ $PROBLEM -eq 176 ]; then
    FILE=mass95-gilman-snap160_vel140_orbit2
    CFG=nws/pickup/live/snap160_vel140_orbit2_mass95-gilman.cfg
fi
###############################################################
# collision between DM sub-halo (Gilman relation) and NW stream (t = 4.5 Gyr)
if [ $PROBLEM -eq 177 ]; then
    FILE=mass95-gilman-snap160_vel140_orbit3
    CFG=nws/pickup/live/snap160_vel140_orbit3_mass95-gilman.cfg
fi
###############################################################
# collision between DM sub-halo (test-particle) and NW stream (t = 5.0 Gyr)
if [ $PROBLEM -eq 200 ]; then
    FILE=mass80-test-snap180_vel140_orbit0
    CFG=nws/pickup/test/snap180_vel140_orbit0_mass80-test.cfg
fi
###############################################################
# collision between DM sub-halo (test-particle) and NW stream (t = 5.0 Gyr)
if [ $PROBLEM -eq 201 ]; then
    FILE=mass80-test-snap180_vel140_orbit1
    CFG=nws/pickup/test/snap180_vel140_orbit1_mass80-test.cfg
fi
###############################################################
# collision between DM sub-halo (test-particle) and NW stream (t = 5.0 Gyr)
if [ $PROBLEM -eq 202 ]; then
    FILE=mass80-test-snap180_vel140_orbit2
    CFG=nws/pickup/test/snap180_vel140_orbit2_mass80-test.cfg
fi
###############################################################
# collision between DM sub-halo (test-particle) and NW stream (t = 5.0 Gyr)
if [ $PROBLEM -eq 203 ]; then
    FILE=mass80-test-snap180_vel140_orbit3
    CFG=nws/pickup/test/snap180_vel140_orbit3_mass80-test.cfg
fi
###############################################################
# collision between DM sub-halo (test-particle) and NW stream (t = 5.0 Gyr)
if [ $PROBLEM -eq 204 ]; then
    FILE=mass85-test-snap180_vel140_orbit0
    CFG=nws/pickup/test/snap180_vel140_orbit0_mass85-test.cfg
fi
###############################################################
# collision between DM sub-halo (test-particle) and NW stream (t = 5.0 Gyr)
if [ $PROBLEM -eq 205 ]; then
    FILE=mass85-test-snap180_vel140_orbit1
    CFG=nws/pickup/test/snap180_vel140_orbit1_mass85-test.cfg
fi
###############################################################
# collision between DM sub-halo (test-particle) and NW stream (t = 5.0 Gyr)
if [ $PROBLEM -eq 206 ]; then
    FILE=mass85-test-snap180_vel140_orbit2
    CFG=nws/pickup/test/snap180_vel140_orbit2_mass85-test.cfg
fi
###############################################################
# collision between DM sub-halo (test-particle) and NW stream (t = 5.0 Gyr)
if [ $PROBLEM -eq 207 ]; then
    FILE=mass85-test-snap180_vel140_orbit3
    CFG=nws/pickup/test/snap180_vel140_orbit3_mass85-test.cfg
fi
###############################################################
# collision between DM sub-halo (test-particle) and NW stream (t = 5.0 Gyr)
if [ $PROBLEM -eq 210 ]; then
    FILE=mass90-test-snap180_vel140_orbit0
    CFG=nws/pickup/test/snap180_vel140_orbit0_mass90-test.cfg
fi
###############################################################
# collision between DM sub-halo (test-particle) and NW stream (t = 5.0 Gyr)
if [ $PROBLEM -eq 211 ]; then
    FILE=mass90-test-snap180_vel140_orbit1
    CFG=nws/pickup/test/snap180_vel140_orbit1_mass90-test.cfg
fi
###############################################################
# collision between DM sub-halo (test-particle) and NW stream (t = 5.0 Gyr)
if [ $PROBLEM -eq 212 ]; then
    FILE=mass90-test-snap180_vel140_orbit2
    CFG=nws/pickup/test/snap180_vel140_orbit2_mass90-test.cfg
fi
###############################################################
# collision between DM sub-halo (test-particle) and NW stream (t = 5.0 Gyr)
if [ $PROBLEM -eq 213 ]; then
    FILE=mass90-test-snap180_vel140_orbit3
    CFG=nws/pickup/test/snap180_vel140_orbit3_mass90-test.cfg
fi
###############################################################
# collision between DM sub-halo (test-particle) and NW stream (t = 5.0 Gyr)
if [ $PROBLEM -eq 214 ]; then
    FILE=mass95-test-snap180_vel140_orbit0
    CFG=nws/pickup/test/snap180_vel140_orbit0_mass95-test.cfg
fi
###############################################################
# collision between DM sub-halo (test-particle) and NW stream (t = 5.0 Gyr)
if [ $PROBLEM -eq 215 ]; then
    FILE=mass95-test-snap180_vel140_orbit1
    CFG=nws/pickup/test/snap180_vel140_orbit1_mass95-test.cfg
fi
###############################################################
# collision between DM sub-halo (test-particle) and NW stream (t = 5.0 Gyr)
if [ $PROBLEM -eq 216 ]; then
    FILE=mass95-test-snap180_vel140_orbit2
    CFG=nws/pickup/test/snap180_vel140_orbit2_mass95-test.cfg
fi
###############################################################
# collision between DM sub-halo (test-particle) and NW stream (t = 5.0 Gyr)
if [ $PROBLEM -eq 217 ]; then
    FILE=mass95-test-snap180_vel140_orbit3
    CFG=nws/pickup/test/snap180_vel140_orbit3_mass95-test.cfg
fi
###############################################################
# collision between DM sub-halo (Prada relation) and NW stream (t = 5.0 Gyr)
if [ $PROBLEM -eq 220 ]; then
    FILE=mass80-prada-snap180_vel140_orbit0
    CFG=nws/pickup/live/snap180_vel140_orbit0_mass80-prada.cfg
fi
###############################################################
# collision between DM sub-halo (Prada relation) and NW stream (t = 5.0 Gyr)
if [ $PROBLEM -eq 221 ]; then
    FILE=mass80-prada-snap180_vel140_orbit1
    CFG=nws/pickup/live/snap180_vel140_orbit1_mass80-prada.cfg
fi
###############################################################
# collision between DM sub-halo (Prada relation) and NW stream (t = 5.0 Gyr)
if [ $PROBLEM -eq 222 ]; then
    FILE=mass80-prada-snap180_vel140_orbit2
    CFG=nws/pickup/live/snap180_vel140_orbit2_mass80-prada.cfg
fi
###############################################################
# collision between DM sub-halo (Prada relation) and NW stream (t = 5.0 Gyr)
if [ $PROBLEM -eq 223 ]; then
    FILE=mass80-prada-snap180_vel140_orbit3
    CFG=nws/pickup/live/snap180_vel140_orbit3_mass80-prada.cfg
fi
###############################################################
# collision between DM sub-halo (Prada relation) and NW stream (t = 5.0 Gyr)
if [ $PROBLEM -eq 224 ]; then
    FILE=mass85-prada-snap180_vel140_orbit0
    CFG=nws/pickup/live/snap180_vel140_orbit0_mass85-prada.cfg
fi
###############################################################
# collision between DM sub-halo (Prada relation) and NW stream (t = 5.0 Gyr)
if [ $PROBLEM -eq 225 ]; then
    FILE=mass85-prada-snap180_vel140_orbit1
    CFG=nws/pickup/live/snap180_vel140_orbit1_mass85-prada.cfg
fi
###############################################################
# collision between DM sub-halo (Prada relation) and NW stream (t = 5.0 Gyr)
if [ $PROBLEM -eq 226 ]; then
    FILE=mass85-prada-snap180_vel140_orbit2
    CFG=nws/pickup/live/snap180_vel140_orbit2_mass85-prada.cfg
fi
###############################################################
# collision between DM sub-halo (Prada relation) and NW stream (t = 5.0 Gyr)
if [ $PROBLEM -eq 227 ]; then
    FILE=mass85-prada-snap180_vel140_orbit3
    CFG=nws/pickup/live/snap180_vel140_orbit3_mass85-prada.cfg
fi
###############################################################
# collision between DM sub-halo (Prada relation) and NW stream (t = 5.0 Gyr)
if [ $PROBLEM -eq 230 ]; then
    FILE=mass90-prada-snap180_vel140_orbit0
    CFG=nws/pickup/live/snap180_vel140_orbit0_mass90-prada.cfg
fi
###############################################################
# collision between DM sub-halo (Prada relation) and NW stream (t = 5.0 Gyr)
if [ $PROBLEM -eq 231 ]; then
    FILE=mass90-prada-snap180_vel140_orbit1
    CFG=nws/pickup/live/snap180_vel140_orbit1_mass90-prada.cfg
fi
###############################################################
# collision between DM sub-halo (Prada relation) and NW stream (t = 5.0 Gyr)
if [ $PROBLEM -eq 232 ]; then
    FILE=mass90-prada-snap180_vel140_orbit2
    CFG=nws/pickup/live/snap180_vel140_orbit2_mass90-prada.cfg
fi
###############################################################
# collision between DM sub-halo (Prada relation) and NW stream (t = 5.0 Gyr)
if [ $PROBLEM -eq 233 ]; then
    FILE=mass90-prada-snap180_vel140_orbit3
    CFG=nws/pickup/live/snap180_vel140_orbit3_mass90-prada.cfg
fi
###############################################################
# collision between DM sub-halo (Prada relation) and NW stream (t = 5.0 Gyr)
if [ $PROBLEM -eq 234 ]; then
    FILE=mass95-prada-snap180_vel140_orbit0
    CFG=nws/pickup/live/snap180_vel140_orbit0_mass95-prada.cfg
fi
###############################################################
# collision between DM sub-halo (Prada relation) and NW stream (t = 5.0 Gyr)
if [ $PROBLEM -eq 235 ]; then
    FILE=mass95-prada-snap180_vel140_orbit1
    CFG=nws/pickup/live/snap180_vel140_orbit1_mass95-prada.cfg
fi
###############################################################
# collision between DM sub-halo (Prada relation) and NW stream (t = 5.0 Gyr)
if [ $PROBLEM -eq 236 ]; then
    FILE=mass95-prada-snap180_vel140_orbit2
    CFG=nws/pickup/live/snap180_vel140_orbit2_mass95-prada.cfg
fi
###############################################################
# collision between DM sub-halo (Prada relation) and NW stream (t = 5.0 Gyr)
if [ $PROBLEM -eq 237 ]; then
    FILE=mass95-prada-snap180_vel140_orbit3
    CFG=nws/pickup/live/snap180_vel140_orbit3_mass95-prada.cfg
fi
###############################################################
# collision between DM sub-halo (Ishiyama & Ando relation) and NW stream (t = 5.0 Gyr)
if [ $PROBLEM -eq 240 ]; then
    FILE=mass80-ishiyama-snap180_vel140_orbit0
    CFG=nws/pickup/live/snap180_vel140_orbit0_mass80-ishiyama.cfg
fi
###############################################################
# collision between DM sub-halo (Ishiyama & Ando relation) and NW stream (t = 5.0 Gyr)
if [ $PROBLEM -eq 241 ]; then
    FILE=mass80-ishiyama-snap180_vel140_orbit1
    CFG=nws/pickup/live/snap180_vel140_orbit1_mass80-ishiyama.cfg
fi
###############################################################
# collision between DM sub-halo (Ishiyama & Ando relation) and NW stream (t = 5.0 Gyr)
if [ $PROBLEM -eq 242 ]; then
    FILE=mass80-ishiyama-snap180_vel140_orbit2
    CFG=nws/pickup/live/snap180_vel140_orbit2_mass80-ishiyama.cfg
fi
###############################################################
# collision between DM sub-halo (Ishiyama & Ando relation) and NW stream (t = 5.0 Gyr)
if [ $PROBLEM -eq 243 ]; then
    FILE=mass80-ishiyama-snap180_vel140_orbit3
    CFG=nws/pickup/live/snap180_vel140_orbit3_mass80-ishiyama.cfg
fi
###############################################################
# collision between DM sub-halo (Ishiyama & Ando relation) and NW stream (t = 5.0 Gyr)
if [ $PROBLEM -eq 244 ]; then
    FILE=mass85-ishiyama-snap180_vel140_orbit0
    CFG=nws/pickup/live/snap180_vel140_orbit0_mass85-ishiyama.cfg
fi
###############################################################
# collision between DM sub-halo (Ishiyama & Ando relation) and NW stream (t = 5.0 Gyr)
if [ $PROBLEM -eq 245 ]; then
    FILE=mass85-ishiyama-snap180_vel140_orbit1
    CFG=nws/pickup/live/snap180_vel140_orbit1_mass85-ishiyama.cfg
fi
###############################################################
# collision between DM sub-halo (Ishiyama & Ando relation) and NW stream (t = 5.0 Gyr)
if [ $PROBLEM -eq 246 ]; then
    FILE=mass85-ishiyama-snap180_vel140_orbit2
    CFG=nws/pickup/live/snap180_vel140_orbit2_mass85-ishiyama.cfg
fi
###############################################################
# collision between DM sub-halo (Ishiyama & Ando relation) and NW stream (t = 5.0 Gyr)
if [ $PROBLEM -eq 247 ]; then
    FILE=mass85-ishiyama-snap180_vel140_orbit3
    CFG=nws/pickup/live/snap180_vel140_orbit3_mass85-ishiyama.cfg
fi
###############################################################
# collision between DM sub-halo (Ishiyama & Ando relation) and NW stream (t = 5.0 Gyr)
if [ $PROBLEM -eq 250 ]; then
    FILE=mass90-ishiyama-snap180_vel140_orbit0
    CFG=nws/pickup/live/snap180_vel140_orbit0_mass90-ishiyama.cfg
fi
###############################################################
# collision between DM sub-halo (Ishiyama & Ando relation) and NW stream (t = 5.0 Gyr)
if [ $PROBLEM -eq 251 ]; then
    FILE=mass90-ishiyama-snap180_vel140_orbit1
    CFG=nws/pickup/live/snap180_vel140_orbit1_mass90-ishiyama.cfg
fi
###############################################################
# collision between DM sub-halo (Ishiyama & Ando relation) and NW stream (t = 5.0 Gyr)
if [ $PROBLEM -eq 252 ]; then
    FILE=mass90-ishiyama-snap180_vel140_orbit2
    CFG=nws/pickup/live/snap180_vel140_orbit2_mass90-ishiyama.cfg
fi
###############################################################
# collision between DM sub-halo (Ishiyama & Ando relation) and NW stream (t = 5.0 Gyr)
if [ $PROBLEM -eq 253 ]; then
    FILE=mass90-ishiyama-snap180_vel140_orbit3
    CFG=nws/pickup/live/snap180_vel140_orbit3_mass90-ishiyama.cfg
fi
###############################################################
# collision between DM sub-halo (Ishiyama & Ando relation) and NW stream (t = 5.0 Gyr)
if [ $PROBLEM -eq 254 ]; then
    FILE=mass95-ishiyama-snap180_vel140_orbit0
    CFG=nws/pickup/live/snap180_vel140_orbit0_mass95-ishiyama.cfg
fi
###############################################################
# collision between DM sub-halo (Ishiyama & Ando relation) and NW stream (t = 5.0 Gyr)
if [ $PROBLEM -eq 255 ]; then
    FILE=mass95-ishiyama-snap180_vel140_orbit1
    CFG=nws/pickup/live/snap180_vel140_orbit1_mass95-ishiyama.cfg
fi
###############################################################
# collision between DM sub-halo (Ishiyama & Ando relation) and NW stream (t = 5.0 Gyr)
if [ $PROBLEM -eq 256 ]; then
    FILE=mass95-ishiyama-snap180_vel140_orbit2
    CFG=nws/pickup/live/snap180_vel140_orbit2_mass95-ishiyama.cfg
fi
###############################################################
# collision between DM sub-halo (Ishiyama & Ando relation) and NW stream (t = 5.0 Gyr)
if [ $PROBLEM -eq 257 ]; then
    FILE=mass95-ishiyama-snap180_vel140_orbit3
    CFG=nws/pickup/live/snap180_vel140_orbit3_mass95-ishiyama.cfg
fi
###############################################################
# collision between DM sub-halo (Gilman relation) and NW stream (t = 5.0 Gyr)
if [ $PROBLEM -eq 260 ]; then
    FILE=mass80-gilman-snap180_vel140_orbit0
    CFG=nws/pickup/live/snap180_vel140_orbit0_mass80-gilman.cfg
fi
###############################################################
# collision between DM sub-halo (Gilman relation) and NW stream (t = 5.0 Gyr)
if [ $PROBLEM -eq 261 ]; then
    FILE=mass80-gilman-snap180_vel140_orbit1
    CFG=nws/pickup/live/snap180_vel140_orbit1_mass80-gilman.cfg
fi
###############################################################
# collision between DM sub-halo (Gilman relation) and NW stream (t = 5.0 Gyr)
if [ $PROBLEM -eq 262 ]; then
    FILE=mass80-gilman-snap180_vel140_orbit2
    CFG=nws/pickup/live/snap180_vel140_orbit2_mass80-gilman.cfg
fi
###############################################################
# collision between DM sub-halo (Gilman relation) and NW stream (t = 5.0 Gyr)
if [ $PROBLEM -eq 263 ]; then
    FILE=mass80-gilman-snap180_vel140_orbit3
    CFG=nws/pickup/live/snap180_vel140_orbit3_mass80-gilman.cfg
fi
###############################################################
# collision between DM sub-halo (Gilman relation) and NW stream (t = 5.0 Gyr)
if [ $PROBLEM -eq 264 ]; then
    FILE=mass85-gilman-snap180_vel140_orbit0
    CFG=nws/pickup/live/snap180_vel140_orbit0_mass85-gilman.cfg
fi
###############################################################
# collision between DM sub-halo (Gilman relation) and NW stream (t = 5.0 Gyr)
if [ $PROBLEM -eq 265 ]; then
    FILE=mass85-gilman-snap180_vel140_orbit1
    CFG=nws/pickup/live/snap180_vel140_orbit1_mass85-gilman.cfg
fi
###############################################################
# collision between DM sub-halo (Gilman relation) and NW stream (t = 5.0 Gyr)
if [ $PROBLEM -eq 266 ]; then
    FILE=mass85-gilman-snap180_vel140_orbit2
    CFG=nws/pickup/live/snap180_vel140_orbit2_mass85-gilman.cfg
fi
###############################################################
# collision between DM sub-halo (Gilman relation) and NW stream (t = 5.0 Gyr)
if [ $PROBLEM -eq 267 ]; then
    FILE=mass85-gilman-snap180_vel140_orbit3
    CFG=nws/pickup/live/snap180_vel140_orbit3_mass85-gilman.cfg
fi
###############################################################
# collision between DM sub-halo (Gilman relation) and NW stream (t = 5.0 Gyr)
if [ $PROBLEM -eq 270 ]; then
    FILE=mass90-gilman-snap180_vel140_orbit0
    CFG=nws/pickup/live/snap180_vel140_orbit0_mass90-gilman.cfg
fi
###############################################################
# collision between DM sub-halo (Gilman relation) and NW stream (t = 5.0 Gyr)
if [ $PROBLEM -eq 271 ]; then
    FILE=mass90-gilman-snap180_vel140_orbit1
    CFG=nws/pickup/live/snap180_vel140_orbit1_mass90-gilman.cfg
fi
###############################################################
# collision between DM sub-halo (Gilman relation) and NW stream (t = 5.0 Gyr)
if [ $PROBLEM -eq 272 ]; then
    FILE=mass90-gilman-snap180_vel140_orbit2
    CFG=nws/pickup/live/snap180_vel140_orbit2_mass90-gilman.cfg
fi
###############################################################
# collision between DM sub-halo (Gilman relation) and NW stream (t = 5.0 Gyr)
if [ $PROBLEM -eq 273 ]; then
    FILE=mass90-gilman-snap180_vel140_orbit3
    CFG=nws/pickup/live/snap180_vel140_orbit3_mass90-gilman.cfg
fi
###############################################################
# collision between DM sub-halo (Gilman relation) and NW stream (t = 5.0 Gyr)
if [ $PROBLEM -eq 274 ]; then
    FILE=mass95-gilman-snap180_vel140_orbit0
    CFG=nws/pickup/live/snap180_vel140_orbit0_mass95-gilman.cfg
fi
###############################################################
# collision between DM sub-halo (Gilman relation) and NW stream (t = 5.0 Gyr)
if [ $PROBLEM -eq 275 ]; then
    FILE=mass95-gilman-snap180_vel140_orbit1
    CFG=nws/pickup/live/snap180_vel140_orbit1_mass95-gilman.cfg
fi
###############################################################
# collision between DM sub-halo (Gilman relation) and NW stream (t = 5.0 Gyr)
if [ $PROBLEM -eq 276 ]; then
    FILE=mass95-gilman-snap180_vel140_orbit2
    CFG=nws/pickup/live/snap180_vel140_orbit2_mass95-gilman.cfg
fi
###############################################################
# collision between DM sub-halo (Gilman relation) and NW stream (t = 5.0 Gyr)
if [ $PROBLEM -eq 277 ]; then
    FILE=mass95-gilman-snap180_vel140_orbit3
    CFG=nws/pickup/live/snap180_vel140_orbit3_mass95-gilman.cfg
fi
###############################################################
# collision between DM sub-halo (test-particle) and NW stream (t = 5.5 Gyr)
if [ $PROBLEM -eq 300 ]; then
    FILE=mass80-test-snap200_vel140_orbit0
    CFG=nws/pickup/test/snap200_vel140_orbit0_mass80-test.cfg
fi
###############################################################
# collision between DM sub-halo (test-particle) and NW stream (t = 5.5 Gyr)
if [ $PROBLEM -eq 301 ]; then
    FILE=mass80-test-snap200_vel140_orbit1
    CFG=nws/pickup/test/snap200_vel140_orbit1_mass80-test.cfg
fi
###############################################################
# collision between DM sub-halo (test-particle) and NW stream (t = 5.5 Gyr)
if [ $PROBLEM -eq 302 ]; then
    FILE=mass80-test-snap200_vel140_orbit2
    CFG=nws/pickup/test/snap200_vel140_orbit2_mass80-test.cfg
fi
###############################################################
# collision between DM sub-halo (test-particle) and NW stream (t = 5.5 Gyr)
if [ $PROBLEM -eq 303 ]; then
    FILE=mass80-test-snap200_vel140_orbit3
    CFG=nws/pickup/test/snap200_vel140_orbit3_mass80-test.cfg
fi
###############################################################
# collision between DM sub-halo (test-particle) and NW stream (t = 5.5 Gyr)
if [ $PROBLEM -eq 304 ]; then
    FILE=mass85-test-snap200_vel140_orbit0
    CFG=nws/pickup/test/snap200_vel140_orbit0_mass85-test.cfg
fi
###############################################################
# collision between DM sub-halo (test-particle) and NW stream (t = 5.5 Gyr)
if [ $PROBLEM -eq 305 ]; then
    FILE=mass85-test-snap200_vel140_orbit1
    CFG=nws/pickup/test/snap200_vel140_orbit1_mass85-test.cfg
fi
###############################################################
# collision between DM sub-halo (test-particle) and NW stream (t = 5.5 Gyr)
if [ $PROBLEM -eq 306 ]; then
    FILE=mass85-test-snap200_vel140_orbit2
    CFG=nws/pickup/test/snap200_vel140_orbit2_mass85-test.cfg
fi
###############################################################
# collision between DM sub-halo (test-particle) and NW stream (t = 5.5 Gyr)
if [ $PROBLEM -eq 307 ]; then
    FILE=mass85-test-snap200_vel140_orbit3
    CFG=nws/pickup/test/snap200_vel140_orbit3_mass85-test.cfg
fi
###############################################################
# collision between DM sub-halo (test-particle) and NW stream (t = 5.5 Gyr)
if [ $PROBLEM -eq 310 ]; then
    FILE=mass90-test-snap200_vel140_orbit0
    CFG=nws/pickup/test/snap200_vel140_orbit0_mass90-test.cfg
fi
###############################################################
# collision between DM sub-halo (test-particle) and NW stream (t = 5.5 Gyr)
if [ $PROBLEM -eq 311 ]; then
    FILE=mass90-test-snap200_vel140_orbit1
    CFG=nws/pickup/test/snap200_vel140_orbit1_mass90-test.cfg
fi
###############################################################
# collision between DM sub-halo (test-particle) and NW stream (t = 5.5 Gyr)
if [ $PROBLEM -eq 312 ]; then
    FILE=mass90-test-snap200_vel140_orbit2
    CFG=nws/pickup/test/snap200_vel140_orbit2_mass90-test.cfg
fi
###############################################################
# collision between DM sub-halo (test-particle) and NW stream (t = 5.5 Gyr)
if [ $PROBLEM -eq 313 ]; then
    FILE=mass90-test-snap200_vel140_orbit3
    CFG=nws/pickup/test/snap200_vel140_orbit3_mass90-test.cfg
fi
###############################################################
# collision between DM sub-halo (test-particle) and NW stream (t = 5.5 Gyr)
if [ $PROBLEM -eq 314 ]; then
    FILE=mass95-test-snap200_vel140_orbit0
    CFG=nws/pickup/test/snap200_vel140_orbit0_mass95-test.cfg
fi
###############################################################
# collision between DM sub-halo (test-particle) and NW stream (t = 5.5 Gyr)
if [ $PROBLEM -eq 315 ]; then
    FILE=mass95-test-snap200_vel140_orbit1
    CFG=nws/pickup/test/snap200_vel140_orbit1_mass95-test.cfg
fi
###############################################################
# collision between DM sub-halo (test-particle) and NW stream (t = 5.5 Gyr)
if [ $PROBLEM -eq 316 ]; then
    FILE=mass95-test-snap200_vel140_orbit2
    CFG=nws/pickup/test/snap200_vel140_orbit2_mass95-test.cfg
fi
###############################################################
# collision between DM sub-halo (test-particle) and NW stream (t = 5.5 Gyr)
if [ $PROBLEM -eq 317 ]; then
    FILE=mass95-test-snap200_vel140_orbit3
    CFG=nws/pickup/test/snap200_vel140_orbit3_mass95-test.cfg
fi
###############################################################
# collision between DM sub-halo (Prada relation) and NW stream (t = 5.5 Gyr)
if [ $PROBLEM -eq 320 ]; then
    FILE=mass80-prada-snap200_vel140_orbit0
    CFG=nws/pickup/live/snap200_vel140_orbit0_mass80-prada.cfg
fi
###############################################################
# collision between DM sub-halo (Prada relation) and NW stream (t = 5.5 Gyr)
if [ $PROBLEM -eq 321 ]; then
    FILE=mass80-prada-snap200_vel140_orbit1
    CFG=nws/pickup/live/snap200_vel140_orbit1_mass80-prada.cfg
fi
###############################################################
# collision between DM sub-halo (Prada relation) and NW stream (t = 5.5 Gyr)
if [ $PROBLEM -eq 322 ]; then
    FILE=mass80-prada-snap200_vel140_orbit2
    CFG=nws/pickup/live/snap200_vel140_orbit2_mass80-prada.cfg
fi
###############################################################
# collision between DM sub-halo (Prada relation) and NW stream (t = 5.5 Gyr)
if [ $PROBLEM -eq 323 ]; then
    FILE=mass80-prada-snap200_vel140_orbit3
    CFG=nws/pickup/live/snap200_vel140_orbit3_mass80-prada.cfg
fi
###############################################################
# collision between DM sub-halo (Prada relation) and NW stream (t = 5.5 Gyr)
if [ $PROBLEM -eq 324 ]; then
    FILE=mass85-prada-snap200_vel140_orbit0
    CFG=nws/pickup/live/snap200_vel140_orbit0_mass85-prada.cfg
fi
###############################################################
# collision between DM sub-halo (Prada relation) and NW stream (t = 5.5 Gyr)
if [ $PROBLEM -eq 325 ]; then
    FILE=mass85-prada-snap200_vel140_orbit1
    CFG=nws/pickup/live/snap200_vel140_orbit1_mass85-prada.cfg
fi
###############################################################
# collision between DM sub-halo (Prada relation) and NW stream (t = 5.5 Gyr)
if [ $PROBLEM -eq 326 ]; then
    FILE=mass85-prada-snap200_vel140_orbit2
    CFG=nws/pickup/live/snap200_vel140_orbit2_mass85-prada.cfg
fi
###############################################################
# collision between DM sub-halo (Prada relation) and NW stream (t = 5.5 Gyr)
if [ $PROBLEM -eq 327 ]; then
    FILE=mass85-prada-snap200_vel140_orbit3
    CFG=nws/pickup/live/snap200_vel140_orbit3_mass85-prada.cfg
fi
###############################################################
# collision between DM sub-halo (Prada relation) and NW stream (t = 5.5 Gyr)
if [ $PROBLEM -eq 330 ]; then
    FILE=mass90-prada-snap200_vel140_orbit0
    CFG=nws/pickup/live/snap200_vel140_orbit0_mass90-prada.cfg
fi
###############################################################
# collision between DM sub-halo (Prada relation) and NW stream (t = 5.5 Gyr)
if [ $PROBLEM -eq 331 ]; then
    FILE=mass90-prada-snap200_vel140_orbit1
    CFG=nws/pickup/live/snap200_vel140_orbit1_mass90-prada.cfg
fi
###############################################################
# collision between DM sub-halo (Prada relation) and NW stream (t = 5.5 Gyr)
if [ $PROBLEM -eq 332 ]; then
    FILE=mass90-prada-snap200_vel140_orbit2
    CFG=nws/pickup/live/snap200_vel140_orbit2_mass90-prada.cfg
fi
###############################################################
# collision between DM sub-halo (Prada relation) and NW stream (t = 5.5 Gyr)
if [ $PROBLEM -eq 333 ]; then
    FILE=mass90-prada-snap200_vel140_orbit3
    CFG=nws/pickup/live/snap200_vel140_orbit3_mass90-prada.cfg
fi
###############################################################
# collision between DM sub-halo (Prada relation) and NW stream (t = 5.5 Gyr)
if [ $PROBLEM -eq 334 ]; then
    FILE=mass95-prada-snap200_vel140_orbit0
    CFG=nws/pickup/live/snap200_vel140_orbit0_mass95-prada.cfg
fi
###############################################################
# collision between DM sub-halo (Prada relation) and NW stream (t = 5.5 Gyr)
if [ $PROBLEM -eq 335 ]; then
    FILE=mass95-prada-snap200_vel140_orbit1
    CFG=nws/pickup/live/snap200_vel140_orbit1_mass95-prada.cfg
fi
###############################################################
# collision between DM sub-halo (Prada relation) and NW stream (t = 5.5 Gyr)
if [ $PROBLEM -eq 336 ]; then
    FILE=mass95-prada-snap200_vel140_orbit2
    CFG=nws/pickup/live/snap200_vel140_orbit2_mass95-prada.cfg
fi
###############################################################
# collision between DM sub-halo (Prada relation) and NW stream (t = 5.5 Gyr)
if [ $PROBLEM -eq 337 ]; then
    FILE=mass95-prada-snap200_vel140_orbit3
    CFG=nws/pickup/live/snap200_vel140_orbit3_mass95-prada.cfg
fi
###############################################################
# collision between DM sub-halo (Ishiyama & Ando relation) and NW stream (t = 5.5 Gyr)
if [ $PROBLEM -eq 340 ]; then
    FILE=mass80-ishiyama-snap200_vel140_orbit0
    CFG=nws/pickup/live/snap200_vel140_orbit0_mass80-ishiyama.cfg
fi
###############################################################
# collision between DM sub-halo (Ishiyama & Ando relation) and NW stream (t = 5.5 Gyr)
if [ $PROBLEM -eq 341 ]; then
    FILE=mass80-ishiyama-snap200_vel140_orbit1
    CFG=nws/pickup/live/snap200_vel140_orbit1_mass80-ishiyama.cfg
fi
###############################################################
# collision between DM sub-halo (Ishiyama & Ando relation) and NW stream (t = 5.5 Gyr)
if [ $PROBLEM -eq 342 ]; then
    FILE=mass80-ishiyama-snap200_vel140_orbit2
    CFG=nws/pickup/live/snap200_vel140_orbit2_mass80-ishiyama.cfg
fi
###############################################################
# collision between DM sub-halo (Ishiyama & Ando relation) and NW stream (t = 5.5 Gyr)
if [ $PROBLEM -eq 343 ]; then
    FILE=mass80-ishiyama-snap200_vel140_orbit3
    CFG=nws/pickup/live/snap200_vel140_orbit3_mass80-ishiyama.cfg
fi
###############################################################
# collision between DM sub-halo (Ishiyama & Ando relation) and NW stream (t = 5.5 Gyr)
if [ $PROBLEM -eq 344 ]; then
    FILE=mass85-ishiyama-snap200_vel140_orbit0
    CFG=nws/pickup/live/snap200_vel140_orbit0_mass85-ishiyama.cfg
fi
###############################################################
# collision between DM sub-halo (Ishiyama & Ando relation) and NW stream (t = 5.5 Gyr)
if [ $PROBLEM -eq 345 ]; then
    FILE=mass85-ishiyama-snap200_vel140_orbit1
    CFG=nws/pickup/live/snap200_vel140_orbit1_mass85-ishiyama.cfg
fi
###############################################################
# collision between DM sub-halo (Ishiyama & Ando relation) and NW stream (t = 5.5 Gyr)
if [ $PROBLEM -eq 346 ]; then
    FILE=mass85-ishiyama-snap200_vel140_orbit2
    CFG=nws/pickup/live/snap200_vel140_orbit2_mass85-ishiyama.cfg
fi
###############################################################
# collision between DM sub-halo (Ishiyama & Ando relation) and NW stream (t = 5.5 Gyr)
if [ $PROBLEM -eq 347 ]; then
    FILE=mass85-ishiyama-snap200_vel140_orbit3
    CFG=nws/pickup/live/snap200_vel140_orbit3_mass85-ishiyama.cfg
fi
###############################################################
# collision between DM sub-halo (Ishiyama & Ando relation) and NW stream (t = 5.5 Gyr)
if [ $PROBLEM -eq 350 ]; then
    FILE=mass90-ishiyama-snap200_vel140_orbit0
    CFG=nws/pickup/live/snap200_vel140_orbit0_mass90-ishiyama.cfg
fi
###############################################################
# collision between DM sub-halo (Ishiyama & Ando relation) and NW stream (t = 5.5 Gyr)
if [ $PROBLEM -eq 351 ]; then
    FILE=mass90-ishiyama-snap200_vel140_orbit1
    CFG=nws/pickup/live/snap200_vel140_orbit1_mass90-ishiyama.cfg
fi
###############################################################
# collision between DM sub-halo (Ishiyama & Ando relation) and NW stream (t = 5.5 Gyr)
if [ $PROBLEM -eq 352 ]; then
    FILE=mass90-ishiyama-snap200_vel140_orbit2
    CFG=nws/pickup/live/snap200_vel140_orbit2_mass90-ishiyama.cfg
fi
###############################################################
# collision between DM sub-halo (Ishiyama & Ando relation) and NW stream (t = 5.5 Gyr)
if [ $PROBLEM -eq 353 ]; then
    FILE=mass90-ishiyama-snap200_vel140_orbit3
    CFG=nws/pickup/live/snap200_vel140_orbit3_mass90-ishiyama.cfg
fi
###############################################################
# collision between DM sub-halo (Ishiyama & Ando relation) and NW stream (t = 5.5 Gyr)
if [ $PROBLEM -eq 354 ]; then
    FILE=mass95-ishiyama-snap200_vel140_orbit0
    CFG=nws/pickup/live/snap200_vel140_orbit0_mass95-ishiyama.cfg
fi
###############################################################
# collision between DM sub-halo (Ishiyama & Ando relation) and NW stream (t = 5.5 Gyr)
if [ $PROBLEM -eq 355 ]; then
    FILE=mass95-ishiyama-snap200_vel140_orbit1
    CFG=nws/pickup/live/snap200_vel140_orbit1_mass95-ishiyama.cfg
fi
###############################################################
# collision between DM sub-halo (Ishiyama & Ando relation) and NW stream (t = 5.5 Gyr)
if [ $PROBLEM -eq 356 ]; then
    FILE=mass95-ishiyama-snap200_vel140_orbit2
    CFG=nws/pickup/live/snap200_vel140_orbit2_mass95-ishiyama.cfg
fi
###############################################################
# collision between DM sub-halo (Ishiyama & Ando relation) and NW stream (t = 5.5 Gyr)
if [ $PROBLEM -eq 357 ]; then
    FILE=mass95-ishiyama-snap200_vel140_orbit3
    CFG=nws/pickup/live/snap200_vel140_orbit3_mass95-ishiyama.cfg
fi
###############################################################
# collision between DM sub-halo (Gilman relation) and NW stream (t = 5.5 Gyr)
if [ $PROBLEM -eq 360 ]; then
    FILE=mass80-gilman-snap200_vel140_orbit0
    CFG=nws/pickup/live/snap200_vel140_orbit0_mass80-gilman.cfg
fi
###############################################################
# collision between DM sub-halo (Gilman relation) and NW stream (t = 5.5 Gyr)
if [ $PROBLEM -eq 361 ]; then
    FILE=mass80-gilman-snap200_vel140_orbit1
    CFG=nws/pickup/live/snap200_vel140_orbit1_mass80-gilman.cfg
fi
###############################################################
# collision between DM sub-halo (Gilman relation) and NW stream (t = 5.5 Gyr)
if [ $PROBLEM -eq 362 ]; then
    FILE=mass80-gilman-snap200_vel140_orbit2
    CFG=nws/pickup/live/snap200_vel140_orbit2_mass80-gilman.cfg
fi
###############################################################
# collision between DM sub-halo (Gilman relation) and NW stream (t = 5.5 Gyr)
if [ $PROBLEM -eq 363 ]; then
    FILE=mass80-gilman-snap200_vel140_orbit3
    CFG=nws/pickup/live/snap200_vel140_orbit3_mass80-gilman.cfg
fi
###############################################################
# collision between DM sub-halo (Gilman relation) and NW stream (t = 5.5 Gyr)
if [ $PROBLEM -eq 364 ]; then
    FILE=mass85-gilman-snap200_vel140_orbit0
    CFG=nws/pickup/live/snap200_vel140_orbit0_mass85-gilman.cfg
fi
###############################################################
# collision between DM sub-halo (Gilman relation) and NW stream (t = 5.5 Gyr)
if [ $PROBLEM -eq 365 ]; then
    FILE=mass85-gilman-snap200_vel140_orbit1
    CFG=nws/pickup/live/snap200_vel140_orbit1_mass85-gilman.cfg
fi
###############################################################
# collision between DM sub-halo (Gilman relation) and NW stream (t = 5.5 Gyr)
if [ $PROBLEM -eq 366 ]; then
    FILE=mass85-gilman-snap200_vel140_orbit2
    CFG=nws/pickup/live/snap200_vel140_orbit2_mass85-gilman.cfg
fi
###############################################################
# collision between DM sub-halo (Gilman relation) and NW stream (t = 5.5 Gyr)
if [ $PROBLEM -eq 367 ]; then
    FILE=mass85-gilman-snap200_vel140_orbit3
    CFG=nws/pickup/live/snap200_vel140_orbit3_mass85-gilman.cfg
fi
###############################################################
# collision between DM sub-halo (Gilman relation) and NW stream (t = 5.5 Gyr)
if [ $PROBLEM -eq 370 ]; then
    FILE=mass90-gilman-snap200_vel140_orbit0
    CFG=nws/pickup/live/snap200_vel140_orbit0_mass90-gilman.cfg
fi
###############################################################
# collision between DM sub-halo (Gilman relation) and NW stream (t = 5.5 Gyr)
if [ $PROBLEM -eq 371 ]; then
    FILE=mass90-gilman-snap200_vel140_orbit1
    CFG=nws/pickup/live/snap200_vel140_orbit1_mass90-gilman.cfg
fi
###############################################################
# collision between DM sub-halo (Gilman relation) and NW stream (t = 5.5 Gyr)
if [ $PROBLEM -eq 372 ]; then
    FILE=mass90-gilman-snap200_vel140_orbit2
    CFG=nws/pickup/live/snap200_vel140_orbit2_mass90-gilman.cfg
fi
###############################################################
# collision between DM sub-halo (Gilman relation) and NW stream (t = 5.5 Gyr)
if [ $PROBLEM -eq 373 ]; then
    FILE=mass90-gilman-snap200_vel140_orbit3
    CFG=nws/pickup/live/snap200_vel140_orbit3_mass90-gilman.cfg
fi
###############################################################
# collision between DM sub-halo (Gilman relation) and NW stream (t = 5.5 Gyr)
if [ $PROBLEM -eq 374 ]; then
    FILE=mass95-gilman-snap200_vel140_orbit0
    CFG=nws/pickup/live/snap200_vel140_orbit0_mass95-gilman.cfg
fi
###############################################################
# collision between DM sub-halo (Gilman relation) and NW stream (t = 5.5 Gyr)
if [ $PROBLEM -eq 375 ]; then
    FILE=mass95-gilman-snap200_vel140_orbit1
    CFG=nws/pickup/live/snap200_vel140_orbit1_mass95-gilman.cfg
fi
###############################################################
# collision between DM sub-halo (Gilman relation) and NW stream (t = 5.5 Gyr)
if [ $PROBLEM -eq 376 ]; then
    FILE=mass95-gilman-snap200_vel140_orbit2
    CFG=nws/pickup/live/snap200_vel140_orbit2_mass95-gilman.cfg
fi
###############################################################
# collision between DM sub-halo (Gilman relation) and NW stream (t = 5.5 Gyr)
if [ $PROBLEM -eq 377 ]; then
    FILE=mass95-gilman-snap200_vel140_orbit3
    CFG=nws/pickup/live/snap200_vel140_orbit3_mass95-gilman.cfg
fi
###############################################################
# collision between DM sub-halo (test-particle) and NW stream (t = 6.0 Gyr)
if [ $PROBLEM -eq 400 ]; then
    FILE=mass80-test-snap220_vel140_orbit0
    CFG=nws/pickup/test/snap220_vel140_orbit0_mass80-test.cfg
fi
###############################################################
# collision between DM sub-halo (test-particle) and NW stream (t = 6.0 Gyr)
if [ $PROBLEM -eq 401 ]; then
    FILE=mass80-test-snap220_vel140_orbit1
    CFG=nws/pickup/test/snap220_vel140_orbit1_mass80-test.cfg
fi
###############################################################
# collision between DM sub-halo (test-particle) and NW stream (t = 6.0 Gyr)
if [ $PROBLEM -eq 402 ]; then
    FILE=mass80-test-snap220_vel140_orbit2
    CFG=nws/pickup/test/snap220_vel140_orbit2_mass80-test.cfg
fi
###############################################################
# collision between DM sub-halo (test-particle) and NW stream (t = 6.0 Gyr)
if [ $PROBLEM -eq 403 ]; then
    FILE=mass80-test-snap220_vel140_orbit3
    CFG=nws/pickup/test/snap220_vel140_orbit3_mass80-test.cfg
fi
###############################################################
# collision between DM sub-halo (test-particle) and NW stream (t = 6.0 Gyr)
if [ $PROBLEM -eq 404 ]; then
    FILE=mass85-test-snap220_vel140_orbit0
    CFG=nws/pickup/test/snap220_vel140_orbit0_mass85-test.cfg
fi
###############################################################
# collision between DM sub-halo (test-particle) and NW stream (t = 6.0 Gyr)
if [ $PROBLEM -eq 405 ]; then
    FILE=mass85-test-snap220_vel140_orbit1
    CFG=nws/pickup/test/snap220_vel140_orbit1_mass85-test.cfg
fi
###############################################################
# collision between DM sub-halo (test-particle) and NW stream (t = 6.0 Gyr)
if [ $PROBLEM -eq 406 ]; then
    FILE=mass85-test-snap220_vel140_orbit2
    CFG=nws/pickup/test/snap220_vel140_orbit2_mass85-test.cfg
fi
###############################################################
# collision between DM sub-halo (test-particle) and NW stream (t = 6.0 Gyr)
if [ $PROBLEM -eq 407 ]; then
    FILE=mass85-test-snap220_vel140_orbit3
    CFG=nws/pickup/test/snap220_vel140_orbit3_mass85-test.cfg
fi
###############################################################
# collision between DM sub-halo (test-particle) and NW stream (t = 6.0 Gyr)
if [ $PROBLEM -eq 410 ]; then
    FILE=mass90-test-snap220_vel140_orbit0
    CFG=nws/pickup/test/snap220_vel140_orbit0_mass90-test.cfg
fi
###############################################################
# collision between DM sub-halo (test-particle) and NW stream (t = 6.0 Gyr)
if [ $PROBLEM -eq 411 ]; then
    FILE=mass90-test-snap220_vel140_orbit1
    CFG=nws/pickup/test/snap220_vel140_orbit1_mass90-test.cfg
fi
###############################################################
# collision between DM sub-halo (test-particle) and NW stream (t = 6.0 Gyr)
if [ $PROBLEM -eq 412 ]; then
    FILE=mass90-test-snap220_vel140_orbit2
    CFG=nws/pickup/test/snap220_vel140_orbit2_mass90-test.cfg
fi
###############################################################
# collision between DM sub-halo (test-particle) and NW stream (t = 6.0 Gyr)
if [ $PROBLEM -eq 413 ]; then
    FILE=mass90-test-snap220_vel140_orbit3
    CFG=nws/pickup/test/snap220_vel140_orbit3_mass90-test.cfg
fi
###############################################################
# collision between DM sub-halo (test-particle) and NW stream (t = 6.0 Gyr)
if [ $PROBLEM -eq 414 ]; then
    FILE=mass95-test-snap220_vel140_orbit0
    CFG=nws/pickup/test/snap220_vel140_orbit0_mass95-test.cfg
fi
###############################################################
# collision between DM sub-halo (test-particle) and NW stream (t = 6.0 Gyr)
if [ $PROBLEM -eq 415 ]; then
    FILE=mass95-test-snap220_vel140_orbit1
    CFG=nws/pickup/test/snap220_vel140_orbit1_mass95-test.cfg
fi
###############################################################
# collision between DM sub-halo (test-particle) and NW stream (t = 6.0 Gyr)
if [ $PROBLEM -eq 416 ]; then
    FILE=mass95-test-snap220_vel140_orbit2
    CFG=nws/pickup/test/snap220_vel140_orbit2_mass95-test.cfg
fi
###############################################################
# collision between DM sub-halo (test-particle) and NW stream (t = 6.0 Gyr)
if [ $PROBLEM -eq 417 ]; then
    FILE=mass95-test-snap220_vel140_orbit3
    CFG=nws/pickup/test/snap220_vel140_orbit3_mass95-test.cfg
fi
###############################################################
# collision between DM sub-halo (Prada relation) and NW stream (t = 6.0 Gyr)
if [ $PROBLEM -eq 420 ]; then
    FILE=mass80-prada-snap220_vel140_orbit0
    CFG=nws/pickup/live/snap220_vel140_orbit0_mass80-prada.cfg
fi
###############################################################
# collision between DM sub-halo (Prada relation) and NW stream (t = 6.0 Gyr)
if [ $PROBLEM -eq 421 ]; then
    FILE=mass80-prada-snap220_vel140_orbit1
    CFG=nws/pickup/live/snap220_vel140_orbit1_mass80-prada.cfg
fi
###############################################################
# collision between DM sub-halo (Prada relation) and NW stream (t = 6.0 Gyr)
if [ $PROBLEM -eq 422 ]; then
    FILE=mass80-prada-snap220_vel140_orbit2
    CFG=nws/pickup/live/snap220_vel140_orbit2_mass80-prada.cfg
fi
###############################################################
# collision between DM sub-halo (Prada relation) and NW stream (t = 6.0 Gyr)
if [ $PROBLEM -eq 423 ]; then
    FILE=mass80-prada-snap220_vel140_orbit3
    CFG=nws/pickup/live/snap220_vel140_orbit3_mass80-prada.cfg
fi
###############################################################
# collision between DM sub-halo (Prada relation) and NW stream (t = 6.0 Gyr)
if [ $PROBLEM -eq 424 ]; then
    FILE=mass85-prada-snap220_vel140_orbit0
    CFG=nws/pickup/live/snap220_vel140_orbit0_mass85-prada.cfg
fi
###############################################################
# collision between DM sub-halo (Prada relation) and NW stream (t = 6.0 Gyr)
if [ $PROBLEM -eq 425 ]; then
    FILE=mass85-prada-snap220_vel140_orbit1
    CFG=nws/pickup/live/snap220_vel140_orbit1_mass85-prada.cfg
fi
###############################################################
# collision between DM sub-halo (Prada relation) and NW stream (t = 6.0 Gyr)
if [ $PROBLEM -eq 426 ]; then
    FILE=mass85-prada-snap220_vel140_orbit2
    CFG=nws/pickup/live/snap220_vel140_orbit2_mass85-prada.cfg
fi
###############################################################
# collision between DM sub-halo (Prada relation) and NW stream (t = 6.0 Gyr)
if [ $PROBLEM -eq 427 ]; then
    FILE=mass85-prada-snap220_vel140_orbit3
    CFG=nws/pickup/live/snap220_vel140_orbit3_mass85-prada.cfg
fi
###############################################################
# collision between DM sub-halo (Prada relation) and NW stream (t = 6.0 Gyr)
if [ $PROBLEM -eq 430 ]; then
    FILE=mass90-prada-snap220_vel140_orbit0
    CFG=nws/pickup/live/snap220_vel140_orbit0_mass90-prada.cfg
fi
###############################################################
# collision between DM sub-halo (Prada relation) and NW stream (t = 6.0 Gyr)
if [ $PROBLEM -eq 431 ]; then
    FILE=mass90-prada-snap220_vel140_orbit1
    CFG=nws/pickup/live/snap220_vel140_orbit1_mass90-prada.cfg
fi
###############################################################
# collision between DM sub-halo (Prada relation) and NW stream (t = 6.0 Gyr)
if [ $PROBLEM -eq 432 ]; then
    FILE=mass90-prada-snap220_vel140_orbit2
    CFG=nws/pickup/live/snap220_vel140_orbit2_mass90-prada.cfg
fi
###############################################################
# collision between DM sub-halo (Prada relation) and NW stream (t = 6.0 Gyr)
if [ $PROBLEM -eq 433 ]; then
    FILE=mass90-prada-snap220_vel140_orbit3
    CFG=nws/pickup/live/snap220_vel140_orbit3_mass90-prada.cfg
fi
###############################################################
# collision between DM sub-halo (Prada relation) and NW stream (t = 6.0 Gyr)
if [ $PROBLEM -eq 434 ]; then
    FILE=mass95-prada-snap220_vel140_orbit0
    CFG=nws/pickup/live/snap220_vel140_orbit0_mass95-prada.cfg
fi
###############################################################
# collision between DM sub-halo (Prada relation) and NW stream (t = 6.0 Gyr)
if [ $PROBLEM -eq 435 ]; then
    FILE=mass95-prada-snap220_vel140_orbit1
    CFG=nws/pickup/live/snap220_vel140_orbit1_mass95-prada.cfg
fi
###############################################################
# collision between DM sub-halo (Prada relation) and NW stream (t = 6.0 Gyr)
if [ $PROBLEM -eq 436 ]; then
    FILE=mass95-prada-snap220_vel140_orbit2
    CFG=nws/pickup/live/snap220_vel140_orbit2_mass95-prada.cfg
fi
###############################################################
# collision between DM sub-halo (Prada relation) and NW stream (t = 6.0 Gyr)
if [ $PROBLEM -eq 437 ]; then
    FILE=mass95-prada-snap220_vel140_orbit3
    CFG=nws/pickup/live/snap220_vel140_orbit3_mass95-prada.cfg
fi
###############################################################
# collision between DM sub-halo (Ishiyama & Ando relation) and NW stream (t = 6.0 Gyr)
if [ $PROBLEM -eq 440 ]; then
    FILE=mass80-ishiyama-snap220_vel140_orbit0
    CFG=nws/pickup/live/snap220_vel140_orbit0_mass80-ishiyama.cfg
fi
###############################################################
# collision between DM sub-halo (Ishiyama & Ando relation) and NW stream (t = 6.0 Gyr)
if [ $PROBLEM -eq 441 ]; then
    FILE=mass80-ishiyama-snap220_vel140_orbit1
    CFG=nws/pickup/live/snap220_vel140_orbit1_mass80-ishiyama.cfg
fi
###############################################################
# collision between DM sub-halo (Ishiyama & Ando relation) and NW stream (t = 6.0 Gyr)
if [ $PROBLEM -eq 442 ]; then
    FILE=mass80-ishiyama-snap220_vel140_orbit2
    CFG=nws/pickup/live/snap220_vel140_orbit2_mass80-ishiyama.cfg
fi
###############################################################
# collision between DM sub-halo (Ishiyama & Ando relation) and NW stream (t = 6.0 Gyr)
if [ $PROBLEM -eq 443 ]; then
    FILE=mass80-ishiyama-snap220_vel140_orbit3
    CFG=nws/pickup/live/snap220_vel140_orbit3_mass80-ishiyama.cfg
fi
###############################################################
# collision between DM sub-halo (Ishiyama & Ando relation) and NW stream (t = 6.0 Gyr)
if [ $PROBLEM -eq 444 ]; then
    FILE=mass85-ishiyama-snap220_vel140_orbit0
    CFG=nws/pickup/live/snap220_vel140_orbit0_mass85-ishiyama.cfg
fi
###############################################################
# collision between DM sub-halo (Ishiyama & Ando relation) and NW stream (t = 6.0 Gyr)
if [ $PROBLEM -eq 445 ]; then
    FILE=mass85-ishiyama-snap220_vel140_orbit1
    CFG=nws/pickup/live/snap220_vel140_orbit1_mass85-ishiyama.cfg
fi
###############################################################
# collision between DM sub-halo (Ishiyama & Ando relation) and NW stream (t = 6.0 Gyr)
if [ $PROBLEM -eq 446 ]; then
    FILE=mass85-ishiyama-snap220_vel140_orbit2
    CFG=nws/pickup/live/snap220_vel140_orbit2_mass85-ishiyama.cfg
fi
###############################################################
# collision between DM sub-halo (Ishiyama & Ando relation) and NW stream (t = 6.0 Gyr)
if [ $PROBLEM -eq 447 ]; then
    FILE=mass85-ishiyama-snap220_vel140_orbit3
    CFG=nws/pickup/live/snap220_vel140_orbit3_mass85-ishiyama.cfg
fi
###############################################################
# collision between DM sub-halo (Ishiyama & Ando relation) and NW stream (t = 6.0 Gyr)
if [ $PROBLEM -eq 450 ]; then
    FILE=mass90-ishiyama-snap220_vel140_orbit0
    CFG=nws/pickup/live/snap220_vel140_orbit0_mass90-ishiyama.cfg
fi
###############################################################
# collision between DM sub-halo (Ishiyama & Ando relation) and NW stream (t = 6.0 Gyr)
if [ $PROBLEM -eq 451 ]; then
    FILE=mass90-ishiyama-snap220_vel140_orbit1
    CFG=nws/pickup/live/snap220_vel140_orbit1_mass90-ishiyama.cfg
fi
###############################################################
# collision between DM sub-halo (Ishiyama & Ando relation) and NW stream (t = 6.0 Gyr)
if [ $PROBLEM -eq 452 ]; then
    FILE=mass90-ishiyama-snap220_vel140_orbit2
    CFG=nws/pickup/live/snap220_vel140_orbit2_mass90-ishiyama.cfg
fi
###############################################################
# collision between DM sub-halo (Ishiyama & Ando relation) and NW stream (t = 6.0 Gyr)
if [ $PROBLEM -eq 453 ]; then
    FILE=mass90-ishiyama-snap220_vel140_orbit3
    CFG=nws/pickup/live/snap220_vel140_orbit3_mass90-ishiyama.cfg
fi
###############################################################
# collision between DM sub-halo (Ishiyama & Ando relation) and NW stream (t = 6.0 Gyr)
if [ $PROBLEM -eq 454 ]; then
    FILE=mass95-ishiyama-snap220_vel140_orbit0
    CFG=nws/pickup/live/snap220_vel140_orbit0_mass95-ishiyama.cfg
fi
###############################################################
# collision between DM sub-halo (Ishiyama & Ando relation) and NW stream (t = 6.0 Gyr)
if [ $PROBLEM -eq 455 ]; then
    FILE=mass95-ishiyama-snap220_vel140_orbit1
    CFG=nws/pickup/live/snap220_vel140_orbit1_mass95-ishiyama.cfg
fi
###############################################################
# collision between DM sub-halo (Ishiyama & Ando relation) and NW stream (t = 6.0 Gyr)
if [ $PROBLEM -eq 456 ]; then
    FILE=mass95-ishiyama-snap220_vel140_orbit2
    CFG=nws/pickup/live/snap220_vel140_orbit2_mass95-ishiyama.cfg
fi
###############################################################
# collision between DM sub-halo (Ishiyama & Ando relation) and NW stream (t = 6.0 Gyr)
if [ $PROBLEM -eq 457 ]; then
    FILE=mass95-ishiyama-snap220_vel140_orbit3
    CFG=nws/pickup/live/snap220_vel140_orbit3_mass95-ishiyama.cfg
fi
###############################################################
# collision between DM sub-halo (Gilman relation) and NW stream (t = 6.0 Gyr)
if [ $PROBLEM -eq 460 ]; then
    FILE=mass80-gilman-snap220_vel140_orbit0
    CFG=nws/pickup/live/snap220_vel140_orbit0_mass80-gilman.cfg
fi
###############################################################
# collision between DM sub-halo (Gilman relation) and NW stream (t = 6.0 Gyr)
if [ $PROBLEM -eq 461 ]; then
    FILE=mass80-gilman-snap220_vel140_orbit1
    CFG=nws/pickup/live/snap220_vel140_orbit1_mass80-gilman.cfg
fi
###############################################################
# collision between DM sub-halo (Gilman relation) and NW stream (t = 6.0 Gyr)
if [ $PROBLEM -eq 462 ]; then
    FILE=mass80-gilman-snap220_vel140_orbit2
    CFG=nws/pickup/live/snap220_vel140_orbit2_mass80-gilman.cfg
fi
###############################################################
# collision between DM sub-halo (Gilman relation) and NW stream (t = 6.0 Gyr)
if [ $PROBLEM -eq 463 ]; then
    FILE=mass80-gilman-snap220_vel140_orbit3
    CFG=nws/pickup/live/snap220_vel140_orbit3_mass80-gilman.cfg
fi
###############################################################
# collision between DM sub-halo (Gilman relation) and NW stream (t = 6.0 Gyr)
if [ $PROBLEM -eq 464 ]; then
    FILE=mass85-gilman-snap220_vel140_orbit0
    CFG=nws/pickup/live/snap220_vel140_orbit0_mass85-gilman.cfg
fi
###############################################################
# collision between DM sub-halo (Gilman relation) and NW stream (t = 6.0 Gyr)
if [ $PROBLEM -eq 465 ]; then
    FILE=mass85-gilman-snap220_vel140_orbit1
    CFG=nws/pickup/live/snap220_vel140_orbit1_mass85-gilman.cfg
fi
###############################################################
# collision between DM sub-halo (Gilman relation) and NW stream (t = 6.0 Gyr)
if [ $PROBLEM -eq 466 ]; then
    FILE=mass85-gilman-snap220_vel140_orbit2
    CFG=nws/pickup/live/snap220_vel140_orbit2_mass85-gilman.cfg
fi
###############################################################
# collision between DM sub-halo (Gilman relation) and NW stream (t = 6.0 Gyr)
if [ $PROBLEM -eq 467 ]; then
    FILE=mass85-gilman-snap220_vel140_orbit3
    CFG=nws/pickup/live/snap220_vel140_orbit3_mass85-gilman.cfg
fi
###############################################################
# collision between DM sub-halo (Gilman relation) and NW stream (t = 6.0 Gyr)
if [ $PROBLEM -eq 470 ]; then
    FILE=mass90-gilman-snap220_vel140_orbit0
    CFG=nws/pickup/live/snap220_vel140_orbit0_mass90-gilman.cfg
fi
###############################################################
# collision between DM sub-halo (Gilman relation) and NW stream (t = 6.0 Gyr)
if [ $PROBLEM -eq 471 ]; then
    FILE=mass90-gilman-snap220_vel140_orbit1
    CFG=nws/pickup/live/snap220_vel140_orbit1_mass90-gilman.cfg
fi
###############################################################
# collision between DM sub-halo (Gilman relation) and NW stream (t = 6.0 Gyr)
if [ $PROBLEM -eq 472 ]; then
    FILE=mass90-gilman-snap220_vel140_orbit2
    CFG=nws/pickup/live/snap220_vel140_orbit2_mass90-gilman.cfg
fi
###############################################################
# collision between DM sub-halo (Gilman relation) and NW stream (t = 6.0 Gyr)
if [ $PROBLEM -eq 473 ]; then
    FILE=mass90-gilman-snap220_vel140_orbit3
    CFG=nws/pickup/live/snap220_vel140_orbit3_mass90-gilman.cfg
fi
###############################################################
# collision between DM sub-halo (Gilman relation) and NW stream (t = 6.0 Gyr)
if [ $PROBLEM -eq 474 ]; then
    FILE=mass95-gilman-snap220_vel140_orbit0
    CFG=nws/pickup/live/snap220_vel140_orbit0_mass95-gilman.cfg
fi
###############################################################
# collision between DM sub-halo (Gilman relation) and NW stream (t = 6.0 Gyr)
if [ $PROBLEM -eq 475 ]; then
    FILE=mass95-gilman-snap220_vel140_orbit1
    CFG=nws/pickup/live/snap220_vel140_orbit1_mass95-gilman.cfg
fi
###############################################################
# collision between DM sub-halo (Gilman relation) and NW stream (t = 6.0 Gyr)
if [ $PROBLEM -eq 476 ]; then
    FILE=mass95-gilman-snap220_vel140_orbit2
    CFG=nws/pickup/live/snap220_vel140_orbit2_mass95-gilman.cfg
fi
###############################################################
# collision between DM sub-halo (Gilman relation) and NW stream (t = 6.0 Gyr)
if [ $PROBLEM -eq 477 ]; then
    FILE=mass95-gilman-snap220_vel140_orbit3
    CFG=nws/pickup/live/snap220_vel140_orbit3_mass95-gilman.cfg
fi
###############################################################
# collision between DM sub-halo (test-particle) and NW stream (t = 6.5 Gyr)
if [ $PROBLEM -eq 500 ]; then
    FILE=mass80-test-snap240_vel140_orbit0
    CFG=nws/pickup/test/snap240_vel140_orbit0_mass80-test.cfg
fi
###############################################################
# collision between DM sub-halo (test-particle) and NW stream (t = 6.5 Gyr)
if [ $PROBLEM -eq 501 ]; then
    FILE=mass80-test-snap240_vel140_orbit1
    CFG=nws/pickup/test/snap240_vel140_orbit1_mass80-test.cfg
fi
###############################################################
# collision between DM sub-halo (test-particle) and NW stream (t = 6.5 Gyr)
if [ $PROBLEM -eq 502 ]; then
    FILE=mass80-test-snap240_vel140_orbit2
    CFG=nws/pickup/test/snap240_vel140_orbit2_mass80-test.cfg
fi
###############################################################
# collision between DM sub-halo (test-particle) and NW stream (t = 6.5 Gyr)
if [ $PROBLEM -eq 503 ]; then
    FILE=mass80-test-snap240_vel140_orbit3
    CFG=nws/pickup/test/snap240_vel140_orbit3_mass80-test.cfg
fi
###############################################################
# collision between DM sub-halo (test-particle) and NW stream (t = 6.5 Gyr)
if [ $PROBLEM -eq 504 ]; then
    FILE=mass85-test-snap240_vel140_orbit0
    CFG=nws/pickup/test/snap240_vel140_orbit0_mass85-test.cfg
fi
###############################################################
# collision between DM sub-halo (test-particle) and NW stream (t = 6.5 Gyr)
if [ $PROBLEM -eq 505 ]; then
    FILE=mass85-test-snap240_vel140_orbit1
    CFG=nws/pickup/test/snap240_vel140_orbit1_mass85-test.cfg
fi
###############################################################
# collision between DM sub-halo (test-particle) and NW stream (t = 6.5 Gyr)
if [ $PROBLEM -eq 506 ]; then
    FILE=mass85-test-snap240_vel140_orbit2
    CFG=nws/pickup/test/snap240_vel140_orbit2_mass85-test.cfg
fi
###############################################################
# collision between DM sub-halo (test-particle) and NW stream (t = 6.5 Gyr)
if [ $PROBLEM -eq 507 ]; then
    FILE=mass85-test-snap240_vel140_orbit3
    CFG=nws/pickup/test/snap240_vel140_orbit3_mass85-test.cfg
fi
###############################################################
# collision between DM sub-halo (test-particle) and NW stream (t = 6.5 Gyr)
if [ $PROBLEM -eq 510 ]; then
    FILE=mass90-test-snap240_vel140_orbit0
    CFG=nws/pickup/test/snap240_vel140_orbit0_mass90-test.cfg
fi
###############################################################
# collision between DM sub-halo (test-particle) and NW stream (t = 6.5 Gyr)
if [ $PROBLEM -eq 511 ]; then
    FILE=mass90-test-snap240_vel140_orbit1
    CFG=nws/pickup/test/snap240_vel140_orbit1_mass90-test.cfg
fi
###############################################################
# collision between DM sub-halo (test-particle) and NW stream (t = 6.5 Gyr)
if [ $PROBLEM -eq 512 ]; then
    FILE=mass90-test-snap240_vel140_orbit2
    CFG=nws/pickup/test/snap240_vel140_orbit2_mass90-test.cfg
fi
###############################################################
# collision between DM sub-halo (test-particle) and NW stream (t = 6.5 Gyr)
if [ $PROBLEM -eq 513 ]; then
    FILE=mass90-test-snap240_vel140_orbit3
    CFG=nws/pickup/test/snap240_vel140_orbit3_mass90-test.cfg
fi
###############################################################
# collision between DM sub-halo (test-particle) and NW stream (t = 6.5 Gyr)
if [ $PROBLEM -eq 514 ]; then
    FILE=mass95-test-snap240_vel140_orbit0
    CFG=nws/pickup/test/snap240_vel140_orbit0_mass95-test.cfg
fi
###############################################################
# collision between DM sub-halo (test-particle) and NW stream (t = 6.5 Gyr)
if [ $PROBLEM -eq 515 ]; then
    FILE=mass95-test-snap240_vel140_orbit1
    CFG=nws/pickup/test/snap240_vel140_orbit1_mass95-test.cfg
fi
###############################################################
# collision between DM sub-halo (test-particle) and NW stream (t = 6.5 Gyr)
if [ $PROBLEM -eq 516 ]; then
    FILE=mass95-test-snap240_vel140_orbit2
    CFG=nws/pickup/test/snap240_vel140_orbit2_mass95-test.cfg
fi
###############################################################
# collision between DM sub-halo (test-particle) and NW stream (t = 6.5 Gyr)
if [ $PROBLEM -eq 517 ]; then
    FILE=mass95-test-snap240_vel140_orbit3
    CFG=nws/pickup/test/snap240_vel140_orbit3_mass95-test.cfg
fi
###############################################################
# collision between DM sub-halo (Prada relation) and NW stream (t = 6.5 Gyr)
if [ $PROBLEM -eq 520 ]; then
    FILE=mass80-prada-snap240_vel140_orbit0
    CFG=nws/pickup/live/snap240_vel140_orbit0_mass80-prada.cfg
fi
###############################################################
# collision between DM sub-halo (Prada relation) and NW stream (t = 6.5 Gyr)
if [ $PROBLEM -eq 521 ]; then
    FILE=mass80-prada-snap240_vel140_orbit1
    CFG=nws/pickup/live/snap240_vel140_orbit1_mass80-prada.cfg
fi
###############################################################
# collision between DM sub-halo (Prada relation) and NW stream (t = 6.5 Gyr)
if [ $PROBLEM -eq 522 ]; then
    FILE=mass80-prada-snap240_vel140_orbit2
    CFG=nws/pickup/live/snap240_vel140_orbit2_mass80-prada.cfg
fi
###############################################################
# collision between DM sub-halo (Prada relation) and NW stream (t = 6.5 Gyr)
if [ $PROBLEM -eq 523 ]; then
    FILE=mass80-prada-snap240_vel140_orbit3
    CFG=nws/pickup/live/snap240_vel140_orbit3_mass80-prada.cfg
fi
###############################################################
# collision between DM sub-halo (Prada relation) and NW stream (t = 6.5 Gyr)
if [ $PROBLEM -eq 524 ]; then
    FILE=mass85-prada-snap240_vel140_orbit0
    CFG=nws/pickup/live/snap240_vel140_orbit0_mass85-prada.cfg
fi
###############################################################
# collision between DM sub-halo (Prada relation) and NW stream (t = 6.5 Gyr)
if [ $PROBLEM -eq 525 ]; then
    FILE=mass85-prada-snap240_vel140_orbit1
    CFG=nws/pickup/live/snap240_vel140_orbit1_mass85-prada.cfg
fi
###############################################################
# collision between DM sub-halo (Prada relation) and NW stream (t = 6.5 Gyr)
if [ $PROBLEM -eq 526 ]; then
    FILE=mass85-prada-snap240_vel140_orbit2
    CFG=nws/pickup/live/snap240_vel140_orbit2_mass85-prada.cfg
fi
###############################################################
# collision between DM sub-halo (Prada relation) and NW stream (t = 6.5 Gyr)
if [ $PROBLEM -eq 527 ]; then
    FILE=mass85-prada-snap240_vel140_orbit3
    CFG=nws/pickup/live/snap240_vel140_orbit3_mass85-prada.cfg
fi
###############################################################
# collision between DM sub-halo (Prada relation) and NW stream (t = 6.5 Gyr)
if [ $PROBLEM -eq 530 ]; then
    FILE=mass90-prada-snap240_vel140_orbit0
    CFG=nws/pickup/live/snap240_vel140_orbit0_mass90-prada.cfg
fi
###############################################################
# collision between DM sub-halo (Prada relation) and NW stream (t = 6.5 Gyr)
if [ $PROBLEM -eq 531 ]; then
    FILE=mass90-prada-snap240_vel140_orbit1
    CFG=nws/pickup/live/snap240_vel140_orbit1_mass90-prada.cfg
fi
###############################################################
# collision between DM sub-halo (Prada relation) and NW stream (t = 6.5 Gyr)
if [ $PROBLEM -eq 532 ]; then
    FILE=mass90-prada-snap240_vel140_orbit2
    CFG=nws/pickup/live/snap240_vel140_orbit2_mass90-prada.cfg
fi
###############################################################
# collision between DM sub-halo (Prada relation) and NW stream (t = 6.5 Gyr)
if [ $PROBLEM -eq 533 ]; then
    FILE=mass90-prada-snap240_vel140_orbit3
    CFG=nws/pickup/live/snap240_vel140_orbit3_mass90-prada.cfg
fi
###############################################################
# collision between DM sub-halo (Prada relation) and NW stream (t = 6.5 Gyr)
if [ $PROBLEM -eq 534 ]; then
    FILE=mass95-prada-snap240_vel140_orbit0
    CFG=nws/pickup/live/snap240_vel140_orbit0_mass95-prada.cfg
fi
###############################################################
# collision between DM sub-halo (Prada relation) and NW stream (t = 6.5 Gyr)
if [ $PROBLEM -eq 535 ]; then
    FILE=mass95-prada-snap240_vel140_orbit1
    CFG=nws/pickup/live/snap240_vel140_orbit1_mass95-prada.cfg
fi
###############################################################
# collision between DM sub-halo (Prada relation) and NW stream (t = 6.5 Gyr)
if [ $PROBLEM -eq 536 ]; then
    FILE=mass95-prada-snap240_vel140_orbit2
    CFG=nws/pickup/live/snap240_vel140_orbit2_mass95-prada.cfg
fi
###############################################################
# collision between DM sub-halo (Prada relation) and NW stream (t = 6.5 Gyr)
if [ $PROBLEM -eq 537 ]; then
    FILE=mass95-prada-snap240_vel140_orbit3
    CFG=nws/pickup/live/snap240_vel140_orbit3_mass95-prada.cfg
fi
###############################################################
# collision between DM sub-halo (Ishiyama & Ando relation) and NW stream (t = 6.5 Gyr)
if [ $PROBLEM -eq 540 ]; then
    FILE=mass80-ishiyama-snap240_vel140_orbit0
    CFG=nws/pickup/live/snap240_vel140_orbit0_mass80-ishiyama.cfg
fi
###############################################################
# collision between DM sub-halo (Ishiyama & Ando relation) and NW stream (t = 6.5 Gyr)
if [ $PROBLEM -eq 541 ]; then
    FILE=mass80-ishiyama-snap240_vel140_orbit1
    CFG=nws/pickup/live/snap240_vel140_orbit1_mass80-ishiyama.cfg
fi
###############################################################
# collision between DM sub-halo (Ishiyama & Ando relation) and NW stream (t = 6.5 Gyr)
if [ $PROBLEM -eq 542 ]; then
    FILE=mass80-ishiyama-snap240_vel140_orbit2
    CFG=nws/pickup/live/snap240_vel140_orbit2_mass80-ishiyama.cfg
fi
###############################################################
# collision between DM sub-halo (Ishiyama & Ando relation) and NW stream (t = 6.5 Gyr)
if [ $PROBLEM -eq 543 ]; then
    FILE=mass80-ishiyama-snap240_vel140_orbit3
    CFG=nws/pickup/live/snap240_vel140_orbit3_mass80-ishiyama.cfg
fi
###############################################################
# collision between DM sub-halo (Ishiyama & Ando relation) and NW stream (t = 6.5 Gyr)
if [ $PROBLEM -eq 544 ]; then
    FILE=mass85-ishiyama-snap240_vel140_orbit0
    CFG=nws/pickup/live/snap240_vel140_orbit0_mass85-ishiyama.cfg
fi
###############################################################
# collision between DM sub-halo (Ishiyama & Ando relation) and NW stream (t = 6.5 Gyr)
if [ $PROBLEM -eq 545 ]; then
    FILE=mass85-ishiyama-snap240_vel140_orbit1
    CFG=nws/pickup/live/snap240_vel140_orbit1_mass85-ishiyama.cfg
fi
###############################################################
# collision between DM sub-halo (Ishiyama & Ando relation) and NW stream (t = 6.5 Gyr)
if [ $PROBLEM -eq 546 ]; then
    FILE=mass85-ishiyama-snap240_vel140_orbit2
    CFG=nws/pickup/live/snap240_vel140_orbit2_mass85-ishiyama.cfg
fi
###############################################################
# collision between DM sub-halo (Ishiyama & Ando relation) and NW stream (t = 6.5 Gyr)
if [ $PROBLEM -eq 547 ]; then
    FILE=mass85-ishiyama-snap240_vel140_orbit3
    CFG=nws/pickup/live/snap240_vel140_orbit3_mass85-ishiyama.cfg
fi
###############################################################
# collision between DM sub-halo (Ishiyama & Ando relation) and NW stream (t = 6.5 Gyr)
if [ $PROBLEM -eq 550 ]; then
    FILE=mass90-ishiyama-snap240_vel140_orbit0
    CFG=nws/pickup/live/snap240_vel140_orbit0_mass90-ishiyama.cfg
fi
###############################################################
# collision between DM sub-halo (Ishiyama & Ando relation) and NW stream (t = 6.5 Gyr)
if [ $PROBLEM -eq 551 ]; then
    FILE=mass90-ishiyama-snap240_vel140_orbit1
    CFG=nws/pickup/live/snap240_vel140_orbit1_mass90-ishiyama.cfg
fi
###############################################################
# collision between DM sub-halo (Ishiyama & Ando relation) and NW stream (t = 6.5 Gyr)
if [ $PROBLEM -eq 552 ]; then
    FILE=mass90-ishiyama-snap240_vel140_orbit2
    CFG=nws/pickup/live/snap240_vel140_orbit2_mass90-ishiyama.cfg
fi
###############################################################
# collision between DM sub-halo (Ishiyama & Ando relation) and NW stream (t = 6.5 Gyr)
if [ $PROBLEM -eq 553 ]; then
    FILE=mass90-ishiyama-snap240_vel140_orbit3
    CFG=nws/pickup/live/snap240_vel140_orbit3_mass90-ishiyama.cfg
fi
###############################################################
# collision between DM sub-halo (Ishiyama & Ando relation) and NW stream (t = 6.5 Gyr)
if [ $PROBLEM -eq 554 ]; then
    FILE=mass95-ishiyama-snap240_vel140_orbit0
    CFG=nws/pickup/live/snap240_vel140_orbit0_mass95-ishiyama.cfg
fi
###############################################################
# collision between DM sub-halo (Ishiyama & Ando relation) and NW stream (t = 6.5 Gyr)
if [ $PROBLEM -eq 555 ]; then
    FILE=mass95-ishiyama-snap240_vel140_orbit1
    CFG=nws/pickup/live/snap240_vel140_orbit1_mass95-ishiyama.cfg
fi
###############################################################
# collision between DM sub-halo (Ishiyama & Ando relation) and NW stream (t = 6.5 Gyr)
if [ $PROBLEM -eq 556 ]; then
    FILE=mass95-ishiyama-snap240_vel140_orbit2
    CFG=nws/pickup/live/snap240_vel140_orbit2_mass95-ishiyama.cfg
fi
###############################################################
# collision between DM sub-halo (Ishiyama & Ando relation) and NW stream (t = 6.5 Gyr)
if [ $PROBLEM -eq 557 ]; then
    FILE=mass95-ishiyama-snap240_vel140_orbit3
    CFG=nws/pickup/live/snap240_vel140_orbit3_mass95-ishiyama.cfg
fi
###############################################################
# collision between DM sub-halo (Gilman relation) and NW stream (t = 6.5 Gyr)
if [ $PROBLEM -eq 560 ]; then
    FILE=mass80-gilman-snap240_vel140_orbit0
    CFG=nws/pickup/live/snap240_vel140_orbit0_mass80-gilman.cfg
fi
###############################################################
# collision between DM sub-halo (Gilman relation) and NW stream (t = 6.5 Gyr)
if [ $PROBLEM -eq 561 ]; then
    FILE=mass80-gilman-snap240_vel140_orbit1
    CFG=nws/pickup/live/snap240_vel140_orbit1_mass80-gilman.cfg
fi
###############################################################
# collision between DM sub-halo (Gilman relation) and NW stream (t = 6.5 Gyr)
if [ $PROBLEM -eq 562 ]; then
    FILE=mass80-gilman-snap240_vel140_orbit2
    CFG=nws/pickup/live/snap240_vel140_orbit2_mass80-gilman.cfg
fi
###############################################################
# collision between DM sub-halo (Gilman relation) and NW stream (t = 6.5 Gyr)
if [ $PROBLEM -eq 563 ]; then
    FILE=mass80-gilman-snap240_vel140_orbit3
    CFG=nws/pickup/live/snap240_vel140_orbit3_mass80-gilman.cfg
fi
###############################################################
# collision between DM sub-halo (Gilman relation) and NW stream (t = 6.5 Gyr)
if [ $PROBLEM -eq 564 ]; then
    FILE=mass85-gilman-snap240_vel140_orbit0
    CFG=nws/pickup/live/snap240_vel140_orbit0_mass85-gilman.cfg
fi
###############################################################
# collision between DM sub-halo (Gilman relation) and NW stream (t = 6.5 Gyr)
if [ $PROBLEM -eq 565 ]; then
    FILE=mass85-gilman-snap240_vel140_orbit1
    CFG=nws/pickup/live/snap240_vel140_orbit1_mass85-gilman.cfg
fi
###############################################################
# collision between DM sub-halo (Gilman relation) and NW stream (t = 6.5 Gyr)
if [ $PROBLEM -eq 566 ]; then
    FILE=mass85-gilman-snap240_vel140_orbit2
    CFG=nws/pickup/live/snap240_vel140_orbit2_mass85-gilman.cfg
fi
###############################################################
# collision between DM sub-halo (Gilman relation) and NW stream (t = 6.5 Gyr)
if [ $PROBLEM -eq 567 ]; then
    FILE=mass85-gilman-snap240_vel140_orbit3
    CFG=nws/pickup/live/snap240_vel140_orbit3_mass85-gilman.cfg
fi
###############################################################
# collision between DM sub-halo (Gilman relation) and NW stream (t = 6.5 Gyr)
if [ $PROBLEM -eq 570 ]; then
    FILE=mass90-gilman-snap240_vel140_orbit0
    CFG=nws/pickup/live/snap240_vel140_orbit0_mass90-gilman.cfg
fi
###############################################################
# collision between DM sub-halo (Gilman relation) and NW stream (t = 6.5 Gyr)
if [ $PROBLEM -eq 571 ]; then
    FILE=mass90-gilman-snap240_vel140_orbit1
    CFG=nws/pickup/live/snap240_vel140_orbit1_mass90-gilman.cfg
fi
###############################################################
# collision between DM sub-halo (Gilman relation) and NW stream (t = 6.5 Gyr)
if [ $PROBLEM -eq 572 ]; then
    FILE=mass90-gilman-snap240_vel140_orbit2
    CFG=nws/pickup/live/snap240_vel140_orbit2_mass90-gilman.cfg
fi
###############################################################
# collision between DM sub-halo (Gilman relation) and NW stream (t = 6.5 Gyr)
if [ $PROBLEM -eq 573 ]; then
    FILE=mass90-gilman-snap240_vel140_orbit3
    CFG=nws/pickup/live/snap240_vel140_orbit3_mass90-gilman.cfg
fi
###############################################################
# collision between DM sub-halo (Gilman relation) and NW stream (t = 6.5 Gyr)
if [ $PROBLEM -eq 574 ]; then
    FILE=mass95-gilman-snap240_vel140_orbit0
    CFG=nws/pickup/live/snap240_vel140_orbit0_mass95-gilman.cfg
fi
###############################################################
# collision between DM sub-halo (Gilman relation) and NW stream (t = 6.5 Gyr)
if [ $PROBLEM -eq 575 ]; then
    FILE=mass95-gilman-snap240_vel140_orbit1
    CFG=nws/pickup/live/snap240_vel140_orbit1_mass95-gilman.cfg
fi
###############################################################
# collision between DM sub-halo (Gilman relation) and NW stream (t = 6.5 Gyr)
if [ $PROBLEM -eq 576 ]; then
    FILE=mass95-gilman-snap240_vel140_orbit2
    CFG=nws/pickup/live/snap240_vel140_orbit2_mass95-gilman.cfg
fi
###############################################################
# collision between DM sub-halo (Gilman relation) and NW stream (t = 6.5 Gyr)
if [ $PROBLEM -eq 577 ]; then
    FILE=mass95-gilman-snap240_vel140_orbit3
    CFG=nws/pickup/live/snap240_vel140_orbit3_mass95-gilman.cfg
fi
###############################################################
# collision between DM sub-halo (test-particle) and NW stream (t = 7.0 Gyr)
if [ $PROBLEM -eq 600 ]; then
    FILE=mass80-test-snap260_vel140_orbit0
    CFG=nws/pickup/test/snap260_vel140_orbit0_mass80-test.cfg
fi
###############################################################
# collision between DM sub-halo (test-particle) and NW stream (t = 7.0 Gyr)
if [ $PROBLEM -eq 601 ]; then
    FILE=mass80-test-snap260_vel140_orbit1
    CFG=nws/pickup/test/snap260_vel140_orbit1_mass80-test.cfg
fi
###############################################################
# collision between DM sub-halo (test-particle) and NW stream (t = 7.0 Gyr)
if [ $PROBLEM -eq 602 ]; then
    FILE=mass80-test-snap260_vel140_orbit2
    CFG=nws/pickup/test/snap260_vel140_orbit2_mass80-test.cfg
fi
###############################################################
# collision between DM sub-halo (test-particle) and NW stream (t = 7.0 Gyr)
if [ $PROBLEM -eq 603 ]; then
    FILE=mass80-test-snap260_vel140_orbit3
    CFG=nws/pickup/test/snap260_vel140_orbit3_mass80-test.cfg
fi
###############################################################
# collision between DM sub-halo (test-particle) and NW stream (t = 7.0 Gyr)
if [ $PROBLEM -eq 604 ]; then
    FILE=mass85-test-snap260_vel140_orbit0
    CFG=nws/pickup/test/snap260_vel140_orbit0_mass85-test.cfg
fi
###############################################################
# collision between DM sub-halo (test-particle) and NW stream (t = 7.0 Gyr)
if [ $PROBLEM -eq 605 ]; then
    FILE=mass85-test-snap260_vel140_orbit1
    CFG=nws/pickup/test/snap260_vel140_orbit1_mass85-test.cfg
fi
###############################################################
# collision between DM sub-halo (test-particle) and NW stream (t = 7.0 Gyr)
if [ $PROBLEM -eq 606 ]; then
    FILE=mass85-test-snap260_vel140_orbit2
    CFG=nws/pickup/test/snap260_vel140_orbit2_mass85-test.cfg
fi
###############################################################
# collision between DM sub-halo (test-particle) and NW stream (t = 7.0 Gyr)
if [ $PROBLEM -eq 607 ]; then
    FILE=mass85-test-snap260_vel140_orbit3
    CFG=nws/pickup/test/snap260_vel140_orbit3_mass85-test.cfg
fi
###############################################################
# collision between DM sub-halo (test-particle) and NW stream (t = 7.0 Gyr)
if [ $PROBLEM -eq 610 ]; then
    FILE=mass90-test-snap260_vel140_orbit0
    CFG=nws/pickup/test/snap260_vel140_orbit0_mass90-test.cfg
fi
###############################################################
# collision between DM sub-halo (test-particle) and NW stream (t = 7.0 Gyr)
if [ $PROBLEM -eq 611 ]; then
    FILE=mass90-test-snap260_vel140_orbit1
    CFG=nws/pickup/test/snap260_vel140_orbit1_mass90-test.cfg
fi
###############################################################
# collision between DM sub-halo (test-particle) and NW stream (t = 7.0 Gyr)
if [ $PROBLEM -eq 612 ]; then
    FILE=mass90-test-snap260_vel140_orbit2
    CFG=nws/pickup/test/snap260_vel140_orbit2_mass90-test.cfg
fi
###############################################################
# collision between DM sub-halo (test-particle) and NW stream (t = 7.0 Gyr)
if [ $PROBLEM -eq 613 ]; then
    FILE=mass90-test-snap260_vel140_orbit3
    CFG=nws/pickup/test/snap260_vel140_orbit3_mass90-test.cfg
fi
###############################################################
# collision between DM sub-halo (test-particle) and NW stream (t = 7.0 Gyr)
if [ $PROBLEM -eq 614 ]; then
    FILE=mass95-test-snap260_vel140_orbit0
    CFG=nws/pickup/test/snap260_vel140_orbit0_mass95-test.cfg
fi
###############################################################
# collision between DM sub-halo (test-particle) and NW stream (t = 7.0 Gyr)
if [ $PROBLEM -eq 615 ]; then
    FILE=mass95-test-snap260_vel140_orbit1
    CFG=nws/pickup/test/snap260_vel140_orbit1_mass95-test.cfg
fi
###############################################################
# collision between DM sub-halo (test-particle) and NW stream (t = 7.0 Gyr)
if [ $PROBLEM -eq 616 ]; then
    FILE=mass95-test-snap260_vel140_orbit2
    CFG=nws/pickup/test/snap260_vel140_orbit2_mass95-test.cfg
fi
###############################################################
# collision between DM sub-halo (test-particle) and NW stream (t = 7.0 Gyr)
if [ $PROBLEM -eq 617 ]; then
    FILE=mass95-test-snap260_vel140_orbit3
    CFG=nws/pickup/test/snap260_vel140_orbit3_mass95-test.cfg
fi
###############################################################
# collision between DM sub-halo (Prada relation) and NW stream (t = 7.0 Gyr)
if [ $PROBLEM -eq 620 ]; then
    FILE=mass80-prada-snap260_vel140_orbit0
    CFG=nws/pickup/live/snap260_vel140_orbit0_mass80-prada.cfg
fi
###############################################################
# collision between DM sub-halo (Prada relation) and NW stream (t = 7.0 Gyr)
if [ $PROBLEM -eq 621 ]; then
    FILE=mass80-prada-snap260_vel140_orbit1
    CFG=nws/pickup/live/snap260_vel140_orbit1_mass80-prada.cfg
fi
###############################################################
# collision between DM sub-halo (Prada relation) and NW stream (t = 7.0 Gyr)
if [ $PROBLEM -eq 622 ]; then
    FILE=mass80-prada-snap260_vel140_orbit2
    CFG=nws/pickup/live/snap260_vel140_orbit2_mass80-prada.cfg
fi
###############################################################
# collision between DM sub-halo (Prada relation) and NW stream (t = 7.0 Gyr)
if [ $PROBLEM -eq 623 ]; then
    FILE=mass80-prada-snap260_vel140_orbit3
    CFG=nws/pickup/live/snap260_vel140_orbit3_mass80-prada.cfg
fi
###############################################################
# collision between DM sub-halo (Prada relation) and NW stream (t = 7.0 Gyr)
if [ $PROBLEM -eq 624 ]; then
    FILE=mass85-prada-snap260_vel140_orbit0
    CFG=nws/pickup/live/snap260_vel140_orbit0_mass85-prada.cfg
fi
###############################################################
# collision between DM sub-halo (Prada relation) and NW stream (t = 7.0 Gyr)
if [ $PROBLEM -eq 625 ]; then
    FILE=mass85-prada-snap260_vel140_orbit1
    CFG=nws/pickup/live/snap260_vel140_orbit1_mass85-prada.cfg
fi
###############################################################
# collision between DM sub-halo (Prada relation) and NW stream (t = 7.0 Gyr)
if [ $PROBLEM -eq 626 ]; then
    FILE=mass85-prada-snap260_vel140_orbit2
    CFG=nws/pickup/live/snap260_vel140_orbit2_mass85-prada.cfg
fi
###############################################################
# collision between DM sub-halo (Prada relation) and NW stream (t = 7.0 Gyr)
if [ $PROBLEM -eq 627 ]; then
    FILE=mass85-prada-snap260_vel140_orbit3
    CFG=nws/pickup/live/snap260_vel140_orbit3_mass85-prada.cfg
fi
###############################################################
# collision between DM sub-halo (Prada relation) and NW stream (t = 7.0 Gyr)
if [ $PROBLEM -eq 630 ]; then
    FILE=mass90-prada-snap260_vel140_orbit0
    CFG=nws/pickup/live/snap260_vel140_orbit0_mass90-prada.cfg
fi
###############################################################
# collision between DM sub-halo (Prada relation) and NW stream (t = 7.0 Gyr)
if [ $PROBLEM -eq 631 ]; then
    FILE=mass90-prada-snap260_vel140_orbit1
    CFG=nws/pickup/live/snap260_vel140_orbit1_mass90-prada.cfg
fi
###############################################################
# collision between DM sub-halo (Prada relation) and NW stream (t = 7.0 Gyr)
if [ $PROBLEM -eq 632 ]; then
    FILE=mass90-prada-snap260_vel140_orbit2
    CFG=nws/pickup/live/snap260_vel140_orbit2_mass90-prada.cfg
fi
###############################################################
# collision between DM sub-halo (Prada relation) and NW stream (t = 7.0 Gyr)
if [ $PROBLEM -eq 633 ]; then
    FILE=mass90-prada-snap260_vel140_orbit3
    CFG=nws/pickup/live/snap260_vel140_orbit3_mass90-prada.cfg
fi
###############################################################
# collision between DM sub-halo (Prada relation) and NW stream (t = 7.0 Gyr)
if [ $PROBLEM -eq 634 ]; then
    FILE=mass95-prada-snap260_vel140_orbit0
    CFG=nws/pickup/live/snap260_vel140_orbit0_mass95-prada.cfg
fi
###############################################################
# collision between DM sub-halo (Prada relation) and NW stream (t = 7.0 Gyr)
if [ $PROBLEM -eq 635 ]; then
    FILE=mass95-prada-snap260_vel140_orbit1
    CFG=nws/pickup/live/snap260_vel140_orbit1_mass95-prada.cfg
fi
###############################################################
# collision between DM sub-halo (Prada relation) and NW stream (t = 7.0 Gyr)
if [ $PROBLEM -eq 636 ]; then
    FILE=mass95-prada-snap260_vel140_orbit2
    CFG=nws/pickup/live/snap260_vel140_orbit2_mass95-prada.cfg
fi
###############################################################
# collision between DM sub-halo (Prada relation) and NW stream (t = 7.0 Gyr)
if [ $PROBLEM -eq 637 ]; then
    FILE=mass95-prada-snap260_vel140_orbit3
    CFG=nws/pickup/live/snap260_vel140_orbit3_mass95-prada.cfg
fi
###############################################################
# collision between DM sub-halo (Ishiyama & Ando relation) and NW stream (t = 7.0 Gyr)
if [ $PROBLEM -eq 640 ]; then
    FILE=mass80-ishiyama-snap260_vel140_orbit0
    CFG=nws/pickup/live/snap260_vel140_orbit0_mass80-ishiyama.cfg
fi
###############################################################
# collision between DM sub-halo (Ishiyama & Ando relation) and NW stream (t = 7.0 Gyr)
if [ $PROBLEM -eq 641 ]; then
    FILE=mass80-ishiyama-snap260_vel140_orbit1
    CFG=nws/pickup/live/snap260_vel140_orbit1_mass80-ishiyama.cfg
fi
###############################################################
# collision between DM sub-halo (Ishiyama & Ando relation) and NW stream (t = 7.0 Gyr)
if [ $PROBLEM -eq 642 ]; then
    FILE=mass80-ishiyama-snap260_vel140_orbit2
    CFG=nws/pickup/live/snap260_vel140_orbit2_mass80-ishiyama.cfg
fi
###############################################################
# collision between DM sub-halo (Ishiyama & Ando relation) and NW stream (t = 7.0 Gyr)
if [ $PROBLEM -eq 643 ]; then
    FILE=mass80-ishiyama-snap260_vel140_orbit3
    CFG=nws/pickup/live/snap260_vel140_orbit3_mass80-ishiyama.cfg
fi
###############################################################
# collision between DM sub-halo (Ishiyama & Ando relation) and NW stream (t = 7.0 Gyr)
if [ $PROBLEM -eq 644 ]; then
    FILE=mass85-ishiyama-snap260_vel140_orbit0
    CFG=nws/pickup/live/snap260_vel140_orbit0_mass85-ishiyama.cfg
fi
###############################################################
# collision between DM sub-halo (Ishiyama & Ando relation) and NW stream (t = 7.0 Gyr)
if [ $PROBLEM -eq 645 ]; then
    FILE=mass85-ishiyama-snap260_vel140_orbit1
    CFG=nws/pickup/live/snap260_vel140_orbit1_mass85-ishiyama.cfg
fi
###############################################################
# collision between DM sub-halo (Ishiyama & Ando relation) and NW stream (t = 7.0 Gyr)
if [ $PROBLEM -eq 646 ]; then
    FILE=mass85-ishiyama-snap260_vel140_orbit2
    CFG=nws/pickup/live/snap260_vel140_orbit2_mass85-ishiyama.cfg
fi
###############################################################
# collision between DM sub-halo (Ishiyama & Ando relation) and NW stream (t = 7.0 Gyr)
if [ $PROBLEM -eq 647 ]; then
    FILE=mass85-ishiyama-snap260_vel140_orbit3
    CFG=nws/pickup/live/snap260_vel140_orbit3_mass85-ishiyama.cfg
fi
###############################################################
# collision between DM sub-halo (Ishiyama & Ando relation) and NW stream (t = 7.0 Gyr)
if [ $PROBLEM -eq 650 ]; then
    FILE=mass90-ishiyama-snap260_vel140_orbit0
    CFG=nws/pickup/live/snap260_vel140_orbit0_mass90-ishiyama.cfg
fi
###############################################################
# collision between DM sub-halo (Ishiyama & Ando relation) and NW stream (t = 7.0 Gyr)
if [ $PROBLEM -eq 651 ]; then
    FILE=mass90-ishiyama-snap260_vel140_orbit1
    CFG=nws/pickup/live/snap260_vel140_orbit1_mass90-ishiyama.cfg
fi
###############################################################
# collision between DM sub-halo (Ishiyama & Ando relation) and NW stream (t = 7.0 Gyr)
if [ $PROBLEM -eq 652 ]; then
    FILE=mass90-ishiyama-snap260_vel140_orbit2
    CFG=nws/pickup/live/snap260_vel140_orbit2_mass90-ishiyama.cfg
fi
###############################################################
# collision between DM sub-halo (Ishiyama & Ando relation) and NW stream (t = 7.0 Gyr)
if [ $PROBLEM -eq 653 ]; then
    FILE=mass90-ishiyama-snap260_vel140_orbit3
    CFG=nws/pickup/live/snap260_vel140_orbit3_mass90-ishiyama.cfg
fi
###############################################################
# collision between DM sub-halo (Ishiyama & Ando relation) and NW stream (t = 7.0 Gyr)
if [ $PROBLEM -eq 654 ]; then
    FILE=mass95-ishiyama-snap260_vel140_orbit0
    CFG=nws/pickup/live/snap260_vel140_orbit0_mass95-ishiyama.cfg
fi
###############################################################
# collision between DM sub-halo (Ishiyama & Ando relation) and NW stream (t = 7.0 Gyr)
if [ $PROBLEM -eq 655 ]; then
    FILE=mass95-ishiyama-snap260_vel140_orbit1
    CFG=nws/pickup/live/snap260_vel140_orbit1_mass95-ishiyama.cfg
fi
###############################################################
# collision between DM sub-halo (Ishiyama & Ando relation) and NW stream (t = 7.0 Gyr)
if [ $PROBLEM -eq 656 ]; then
    FILE=mass95-ishiyama-snap260_vel140_orbit2
    CFG=nws/pickup/live/snap260_vel140_orbit2_mass95-ishiyama.cfg
fi
###############################################################
# collision between DM sub-halo (Ishiyama & Ando relation) and NW stream (t = 7.0 Gyr)
if [ $PROBLEM -eq 657 ]; then
    FILE=mass95-ishiyama-snap260_vel140_orbit3
    CFG=nws/pickup/live/snap260_vel140_orbit3_mass95-ishiyama.cfg
fi
###############################################################
# collision between DM sub-halo (Gilman relation) and NW stream (t = 7.0 Gyr)
if [ $PROBLEM -eq 660 ]; then
    FILE=mass80-gilman-snap260_vel140_orbit0
    CFG=nws/pickup/live/snap260_vel140_orbit0_mass80-gilman.cfg
fi
###############################################################
# collision between DM sub-halo (Gilman relation) and NW stream (t = 7.0 Gyr)
if [ $PROBLEM -eq 661 ]; then
    FILE=mass80-gilman-snap260_vel140_orbit1
    CFG=nws/pickup/live/snap260_vel140_orbit1_mass80-gilman.cfg
fi
###############################################################
# collision between DM sub-halo (Gilman relation) and NW stream (t = 7.0 Gyr)
if [ $PROBLEM -eq 662 ]; then
    FILE=mass80-gilman-snap260_vel140_orbit2
    CFG=nws/pickup/live/snap260_vel140_orbit2_mass80-gilman.cfg
fi
###############################################################
# collision between DM sub-halo (Gilman relation) and NW stream (t = 7.0 Gyr)
if [ $PROBLEM -eq 663 ]; then
    FILE=mass80-gilman-snap260_vel140_orbit3
    CFG=nws/pickup/live/snap260_vel140_orbit3_mass80-gilman.cfg
fi
###############################################################
# collision between DM sub-halo (Gilman relation) and NW stream (t = 7.0 Gyr)
if [ $PROBLEM -eq 664 ]; then
    FILE=mass85-gilman-snap260_vel140_orbit0
    CFG=nws/pickup/live/snap260_vel140_orbit0_mass85-gilman.cfg
fi
###############################################################
# collision between DM sub-halo (Gilman relation) and NW stream (t = 7.0 Gyr)
if [ $PROBLEM -eq 665 ]; then
    FILE=mass85-gilman-snap260_vel140_orbit1
    CFG=nws/pickup/live/snap260_vel140_orbit1_mass85-gilman.cfg
fi
###############################################################
# collision between DM sub-halo (Gilman relation) and NW stream (t = 7.0 Gyr)
if [ $PROBLEM -eq 666 ]; then
    FILE=mass85-gilman-snap260_vel140_orbit2
    CFG=nws/pickup/live/snap260_vel140_orbit2_mass85-gilman.cfg
fi
###############################################################
# collision between DM sub-halo (Gilman relation) and NW stream (t = 7.0 Gyr)
if [ $PROBLEM -eq 667 ]; then
    FILE=mass85-gilman-snap260_vel140_orbit3
    CFG=nws/pickup/live/snap260_vel140_orbit3_mass85-gilman.cfg
fi
###############################################################
# collision between DM sub-halo (Gilman relation) and NW stream (t = 7.0 Gyr)
if [ $PROBLEM -eq 670 ]; then
    FILE=mass90-gilman-snap260_vel140_orbit0
    CFG=nws/pickup/live/snap260_vel140_orbit0_mass90-gilman.cfg
fi
###############################################################
# collision between DM sub-halo (Gilman relation) and NW stream (t = 7.0 Gyr)
if [ $PROBLEM -eq 671 ]; then
    FILE=mass90-gilman-snap260_vel140_orbit1
    CFG=nws/pickup/live/snap260_vel140_orbit1_mass90-gilman.cfg
fi
###############################################################
# collision between DM sub-halo (Gilman relation) and NW stream (t = 7.0 Gyr)
if [ $PROBLEM -eq 672 ]; then
    FILE=mass90-gilman-snap260_vel140_orbit2
    CFG=nws/pickup/live/snap260_vel140_orbit2_mass90-gilman.cfg
fi
###############################################################
# collision between DM sub-halo (Gilman relation) and NW stream (t = 7.0 Gyr)
if [ $PROBLEM -eq 673 ]; then
    FILE=mass90-gilman-snap260_vel140_orbit3
    CFG=nws/pickup/live/snap260_vel140_orbit3_mass90-gilman.cfg
fi
###############################################################
# collision between DM sub-halo (Gilman relation) and NW stream (t = 7.0 Gyr)
if [ $PROBLEM -eq 674 ]; then
    FILE=mass95-gilman-snap260_vel140_orbit0
    CFG=nws/pickup/live/snap260_vel140_orbit0_mass95-gilman.cfg
fi
###############################################################
# collision between DM sub-halo (Gilman relation) and NW stream (t = 7.0 Gyr)
if [ $PROBLEM -eq 675 ]; then
    FILE=mass95-gilman-snap260_vel140_orbit1
    CFG=nws/pickup/live/snap260_vel140_orbit1_mass95-gilman.cfg
fi
###############################################################
# collision between DM sub-halo (Gilman relation) and NW stream (t = 7.0 Gyr)
if [ $PROBLEM -eq 676 ]; then
    FILE=mass95-gilman-snap260_vel140_orbit2
    CFG=nws/pickup/live/snap260_vel140_orbit2_mass95-gilman.cfg
fi
###############################################################
# collision between DM sub-halo (Gilman relation) and NW stream (t = 7.0 Gyr)
if [ $PROBLEM -eq 677 ]; then
    FILE=mass95-gilman-snap260_vel140_orbit3
    CFG=nws/pickup/live/snap260_vel140_orbit3_mass95-gilman.cfg
fi
###############################################################
# collision between DM sub-halo (test-particle) and NW stream (t = 7.5 Gyr)
if [ $PROBLEM -eq 700 ]; then
    FILE=mass80-test-snap280_vel140_orbit0
    CFG=nws/pickup/test/snap280_vel140_orbit0_mass80-test.cfg
fi
###############################################################
# collision between DM sub-halo (test-particle) and NW stream (t = 7.5 Gyr)
if [ $PROBLEM -eq 701 ]; then
    FILE=mass80-test-snap280_vel140_orbit1
    CFG=nws/pickup/test/snap280_vel140_orbit1_mass80-test.cfg
fi
###############################################################
# collision between DM sub-halo (test-particle) and NW stream (t = 7.5 Gyr)
if [ $PROBLEM -eq 702 ]; then
    FILE=mass80-test-snap280_vel140_orbit2
    CFG=nws/pickup/test/snap280_vel140_orbit2_mass80-test.cfg
fi
###############################################################
# collision between DM sub-halo (test-particle) and NW stream (t = 7.5 Gyr)
if [ $PROBLEM -eq 703 ]; then
    FILE=mass80-test-snap280_vel140_orbit3
    CFG=nws/pickup/test/snap280_vel140_orbit3_mass80-test.cfg
fi
###############################################################
# collision between DM sub-halo (test-particle) and NW stream (t = 7.5 Gyr)
if [ $PROBLEM -eq 704 ]; then
    FILE=mass85-test-snap280_vel140_orbit0
    CFG=nws/pickup/test/snap280_vel140_orbit0_mass85-test.cfg
fi
###############################################################
# collision between DM sub-halo (test-particle) and NW stream (t = 7.5 Gyr)
if [ $PROBLEM -eq 705 ]; then
    FILE=mass85-test-snap280_vel140_orbit1
    CFG=nws/pickup/test/snap280_vel140_orbit1_mass85-test.cfg
fi
###############################################################
# collision between DM sub-halo (test-particle) and NW stream (t = 7.5 Gyr)
if [ $PROBLEM -eq 706 ]; then
    FILE=mass85-test-snap280_vel140_orbit2
    CFG=nws/pickup/test/snap280_vel140_orbit2_mass85-test.cfg
fi
###############################################################
# collision between DM sub-halo (test-particle) and NW stream (t = 7.5 Gyr)
if [ $PROBLEM -eq 707 ]; then
    FILE=mass85-test-snap280_vel140_orbit3
    CFG=nws/pickup/test/snap280_vel140_orbit3_mass85-test.cfg
fi
###############################################################
# collision between DM sub-halo (test-particle) and NW stream (t = 7.5 Gyr)
if [ $PROBLEM -eq 710 ]; then
    FILE=mass90-test-snap280_vel140_orbit0
    CFG=nws/pickup/test/snap280_vel140_orbit0_mass90-test.cfg
fi
###############################################################
# collision between DM sub-halo (test-particle) and NW stream (t = 7.5 Gyr)
if [ $PROBLEM -eq 711 ]; then
    FILE=mass90-test-snap280_vel140_orbit1
    CFG=nws/pickup/test/snap280_vel140_orbit1_mass90-test.cfg
fi
###############################################################
# collision between DM sub-halo (test-particle) and NW stream (t = 7.5 Gyr)
if [ $PROBLEM -eq 712 ]; then
    FILE=mass90-test-snap280_vel140_orbit2
    CFG=nws/pickup/test/snap280_vel140_orbit2_mass90-test.cfg
fi
###############################################################
# collision between DM sub-halo (test-particle) and NW stream (t = 7.5 Gyr)
if [ $PROBLEM -eq 713 ]; then
    FILE=mass90-test-snap280_vel140_orbit3
    CFG=nws/pickup/test/snap280_vel140_orbit3_mass90-test.cfg
fi
###############################################################
# collision between DM sub-halo (test-particle) and NW stream (t = 7.5 Gyr)
if [ $PROBLEM -eq 714 ]; then
    FILE=mass95-test-snap280_vel140_orbit0
    CFG=nws/pickup/test/snap280_vel140_orbit0_mass95-test.cfg
fi
###############################################################
# collision between DM sub-halo (test-particle) and NW stream (t = 7.5 Gyr)
if [ $PROBLEM -eq 715 ]; then
    FILE=mass95-test-snap280_vel140_orbit1
    CFG=nws/pickup/test/snap280_vel140_orbit1_mass95-test.cfg
fi
###############################################################
# collision between DM sub-halo (test-particle) and NW stream (t = 7.5 Gyr)
if [ $PROBLEM -eq 716 ]; then
    FILE=mass95-test-snap280_vel140_orbit2
    CFG=nws/pickup/test/snap280_vel140_orbit2_mass95-test.cfg
fi
###############################################################
# collision between DM sub-halo (test-particle) and NW stream (t = 7.5 Gyr)
if [ $PROBLEM -eq 717 ]; then
    FILE=mass95-test-snap280_vel140_orbit3
    CFG=nws/pickup/test/snap280_vel140_orbit3_mass95-test.cfg
fi
###############################################################
# collision between DM sub-halo (Prada relation) and NW stream (t = 7.5 Gyr)
if [ $PROBLEM -eq 720 ]; then
    FILE=mass80-prada-snap280_vel140_orbit0
    CFG=nws/pickup/live/snap280_vel140_orbit0_mass80-prada.cfg
fi
###############################################################
# collision between DM sub-halo (Prada relation) and NW stream (t = 7.5 Gyr)
if [ $PROBLEM -eq 721 ]; then
    FILE=mass80-prada-snap280_vel140_orbit1
    CFG=nws/pickup/live/snap280_vel140_orbit1_mass80-prada.cfg
fi
###############################################################
# collision between DM sub-halo (Prada relation) and NW stream (t = 7.5 Gyr)
if [ $PROBLEM -eq 722 ]; then
    FILE=mass80-prada-snap280_vel140_orbit2
    CFG=nws/pickup/live/snap280_vel140_orbit2_mass80-prada.cfg
fi
###############################################################
# collision between DM sub-halo (Prada relation) and NW stream (t = 7.5 Gyr)
if [ $PROBLEM -eq 723 ]; then
    FILE=mass80-prada-snap280_vel140_orbit3
    CFG=nws/pickup/live/snap280_vel140_orbit3_mass80-prada.cfg
fi
###############################################################
# collision between DM sub-halo (Prada relation) and NW stream (t = 7.5 Gyr)
if [ $PROBLEM -eq 724 ]; then
    FILE=mass85-prada-snap280_vel140_orbit0
    CFG=nws/pickup/live/snap280_vel140_orbit0_mass85-prada.cfg
fi
###############################################################
# collision between DM sub-halo (Prada relation) and NW stream (t = 7.5 Gyr)
if [ $PROBLEM -eq 725 ]; then
    FILE=mass85-prada-snap280_vel140_orbit1
    CFG=nws/pickup/live/snap280_vel140_orbit1_mass85-prada.cfg
fi
###############################################################
# collision between DM sub-halo (Prada relation) and NW stream (t = 7.5 Gyr)
if [ $PROBLEM -eq 726 ]; then
    FILE=mass85-prada-snap280_vel140_orbit2
    CFG=nws/pickup/live/snap280_vel140_orbit2_mass85-prada.cfg
fi
###############################################################
# collision between DM sub-halo (Prada relation) and NW stream (t = 7.5 Gyr)
if [ $PROBLEM -eq 727 ]; then
    FILE=mass85-prada-snap280_vel140_orbit3
    CFG=nws/pickup/live/snap280_vel140_orbit3_mass85-prada.cfg
fi
###############################################################
# collision between DM sub-halo (Prada relation) and NW stream (t = 7.5 Gyr)
if [ $PROBLEM -eq 730 ]; then
    FILE=mass90-prada-snap280_vel140_orbit0
    CFG=nws/pickup/live/snap280_vel140_orbit0_mass90-prada.cfg
fi
###############################################################
# collision between DM sub-halo (Prada relation) and NW stream (t = 7.5 Gyr)
if [ $PROBLEM -eq 731 ]; then
    FILE=mass90-prada-snap280_vel140_orbit1
    CFG=nws/pickup/live/snap280_vel140_orbit1_mass90-prada.cfg
fi
###############################################################
# collision between DM sub-halo (Prada relation) and NW stream (t = 7.5 Gyr)
if [ $PROBLEM -eq 732 ]; then
    FILE=mass90-prada-snap280_vel140_orbit2
    CFG=nws/pickup/live/snap280_vel140_orbit2_mass90-prada.cfg
fi
###############################################################
# collision between DM sub-halo (Prada relation) and NW stream (t = 7.5 Gyr)
if [ $PROBLEM -eq 733 ]; then
    FILE=mass90-prada-snap280_vel140_orbit3
    CFG=nws/pickup/live/snap280_vel140_orbit3_mass90-prada.cfg
fi
###############################################################
# collision between DM sub-halo (Prada relation) and NW stream (t = 7.5 Gyr)
if [ $PROBLEM -eq 734 ]; then
    FILE=mass95-prada-snap280_vel140_orbit0
    CFG=nws/pickup/live/snap280_vel140_orbit0_mass95-prada.cfg
fi
###############################################################
# collision between DM sub-halo (Prada relation) and NW stream (t = 7.5 Gyr)
if [ $PROBLEM -eq 735 ]; then
    FILE=mass95-prada-snap280_vel140_orbit1
    CFG=nws/pickup/live/snap280_vel140_orbit1_mass95-prada.cfg
fi
###############################################################
# collision between DM sub-halo (Prada relation) and NW stream (t = 7.5 Gyr)
if [ $PROBLEM -eq 736 ]; then
    FILE=mass95-prada-snap280_vel140_orbit2
    CFG=nws/pickup/live/snap280_vel140_orbit2_mass95-prada.cfg
fi
###############################################################
# collision between DM sub-halo (Prada relation) and NW stream (t = 7.5 Gyr)
if [ $PROBLEM -eq 737 ]; then
    FILE=mass95-prada-snap280_vel140_orbit3
    CFG=nws/pickup/live/snap280_vel140_orbit3_mass95-prada.cfg
fi
###############################################################
# collision between DM sub-halo (Ishiyama & Ando relation) and NW stream (t = 7.5 Gyr)
if [ $PROBLEM -eq 740 ]; then
    FILE=mass80-ishiyama-snap280_vel140_orbit0
    CFG=nws/pickup/live/snap280_vel140_orbit0_mass80-ishiyama.cfg
fi
###############################################################
# collision between DM sub-halo (Ishiyama & Ando relation) and NW stream (t = 7.5 Gyr)
if [ $PROBLEM -eq 741 ]; then
    FILE=mass80-ishiyama-snap280_vel140_orbit1
    CFG=nws/pickup/live/snap280_vel140_orbit1_mass80-ishiyama.cfg
fi
###############################################################
# collision between DM sub-halo (Ishiyama & Ando relation) and NW stream (t = 7.5 Gyr)
if [ $PROBLEM -eq 742 ]; then
    FILE=mass80-ishiyama-snap280_vel140_orbit2
    CFG=nws/pickup/live/snap280_vel140_orbit2_mass80-ishiyama.cfg
fi
###############################################################
# collision between DM sub-halo (Ishiyama & Ando relation) and NW stream (t = 7.5 Gyr)
if [ $PROBLEM -eq 743 ]; then
    FILE=mass80-ishiyama-snap280_vel140_orbit3
    CFG=nws/pickup/live/snap280_vel140_orbit3_mass80-ishiyama.cfg
fi
###############################################################
# collision between DM sub-halo (Ishiyama & Ando relation) and NW stream (t = 7.5 Gyr)
if [ $PROBLEM -eq 744 ]; then
    FILE=mass85-ishiyama-snap280_vel140_orbit0
    CFG=nws/pickup/live/snap280_vel140_orbit0_mass85-ishiyama.cfg
fi
###############################################################
# collision between DM sub-halo (Ishiyama & Ando relation) and NW stream (t = 7.5 Gyr)
if [ $PROBLEM -eq 745 ]; then
    FILE=mass85-ishiyama-snap280_vel140_orbit1
    CFG=nws/pickup/live/snap280_vel140_orbit1_mass85-ishiyama.cfg
fi
###############################################################
# collision between DM sub-halo (Ishiyama & Ando relation) and NW stream (t = 7.5 Gyr)
if [ $PROBLEM -eq 746 ]; then
    FILE=mass85-ishiyama-snap280_vel140_orbit2
    CFG=nws/pickup/live/snap280_vel140_orbit2_mass85-ishiyama.cfg
fi
###############################################################
# collision between DM sub-halo (Ishiyama & Ando relation) and NW stream (t = 7.5 Gyr)
if [ $PROBLEM -eq 747 ]; then
    FILE=mass85-ishiyama-snap280_vel140_orbit3
    CFG=nws/pickup/live/snap280_vel140_orbit3_mass85-ishiyama.cfg
fi
###############################################################
# collision between DM sub-halo (Ishiyama & Ando relation) and NW stream (t = 7.5 Gyr)
if [ $PROBLEM -eq 750 ]; then
    FILE=mass90-ishiyama-snap280_vel140_orbit0
    CFG=nws/pickup/live/snap280_vel140_orbit0_mass90-ishiyama.cfg
fi
###############################################################
# collision between DM sub-halo (Ishiyama & Ando relation) and NW stream (t = 7.5 Gyr)
if [ $PROBLEM -eq 751 ]; then
    FILE=mass90-ishiyama-snap280_vel140_orbit1
    CFG=nws/pickup/live/snap280_vel140_orbit1_mass90-ishiyama.cfg
fi
###############################################################
# collision between DM sub-halo (Ishiyama & Ando relation) and NW stream (t = 7.5 Gyr)
if [ $PROBLEM -eq 752 ]; then
    FILE=mass90-ishiyama-snap280_vel140_orbit2
    CFG=nws/pickup/live/snap280_vel140_orbit2_mass90-ishiyama.cfg
fi
###############################################################
# collision between DM sub-halo (Ishiyama & Ando relation) and NW stream (t = 7.5 Gyr)
if [ $PROBLEM -eq 753 ]; then
    FILE=mass90-ishiyama-snap280_vel140_orbit3
    CFG=nws/pickup/live/snap280_vel140_orbit3_mass90-ishiyama.cfg
fi
###############################################################
# collision between DM sub-halo (Ishiyama & Ando relation) and NW stream (t = 7.5 Gyr)
if [ $PROBLEM -eq 754 ]; then
    FILE=mass95-ishiyama-snap280_vel140_orbit0
    CFG=nws/pickup/live/snap280_vel140_orbit0_mass95-ishiyama.cfg
fi
###############################################################
# collision between DM sub-halo (Ishiyama & Ando relation) and NW stream (t = 7.5 Gyr)
if [ $PROBLEM -eq 755 ]; then
    FILE=mass95-ishiyama-snap280_vel140_orbit1
    CFG=nws/pickup/live/snap280_vel140_orbit1_mass95-ishiyama.cfg
fi
###############################################################
# collision between DM sub-halo (Ishiyama & Ando relation) and NW stream (t = 7.5 Gyr)
if [ $PROBLEM -eq 756 ]; then
    FILE=mass95-ishiyama-snap280_vel140_orbit2
    CFG=nws/pickup/live/snap280_vel140_orbit2_mass95-ishiyama.cfg
fi
###############################################################
# collision between DM sub-halo (Ishiyama & Ando relation) and NW stream (t = 7.5 Gyr)
if [ $PROBLEM -eq 757 ]; then
    FILE=mass95-ishiyama-snap280_vel140_orbit3
    CFG=nws/pickup/live/snap280_vel140_orbit3_mass95-ishiyama.cfg
fi
###############################################################
# collision between DM sub-halo (Gilman relation) and NW stream (t = 7.5 Gyr)
if [ $PROBLEM -eq 760 ]; then
    FILE=mass80-gilman-snap280_vel140_orbit0
    CFG=nws/pickup/live/snap280_vel140_orbit0_mass80-gilman.cfg
fi
###############################################################
# collision between DM sub-halo (Gilman relation) and NW stream (t = 7.5 Gyr)
if [ $PROBLEM -eq 761 ]; then
    FILE=mass80-gilman-snap280_vel140_orbit1
    CFG=nws/pickup/live/snap280_vel140_orbit1_mass80-gilman.cfg
fi
###############################################################
# collision between DM sub-halo (Gilman relation) and NW stream (t = 7.5 Gyr)
if [ $PROBLEM -eq 762 ]; then
    FILE=mass80-gilman-snap280_vel140_orbit2
    CFG=nws/pickup/live/snap280_vel140_orbit2_mass80-gilman.cfg
fi
###############################################################
# collision between DM sub-halo (Gilman relation) and NW stream (t = 7.5 Gyr)
if [ $PROBLEM -eq 763 ]; then
    FILE=mass80-gilman-snap280_vel140_orbit3
    CFG=nws/pickup/live/snap280_vel140_orbit3_mass80-gilman.cfg
fi
###############################################################
# collision between DM sub-halo (Gilman relation) and NW stream (t = 7.5 Gyr)
if [ $PROBLEM -eq 764 ]; then
    FILE=mass85-gilman-snap280_vel140_orbit0
    CFG=nws/pickup/live/snap280_vel140_orbit0_mass85-gilman.cfg
fi
###############################################################
# collision between DM sub-halo (Gilman relation) and NW stream (t = 7.5 Gyr)
if [ $PROBLEM -eq 765 ]; then
    FILE=mass85-gilman-snap280_vel140_orbit1
    CFG=nws/pickup/live/snap280_vel140_orbit1_mass85-gilman.cfg
fi
###############################################################
# collision between DM sub-halo (Gilman relation) and NW stream (t = 7.5 Gyr)
if [ $PROBLEM -eq 766 ]; then
    FILE=mass85-gilman-snap280_vel140_orbit2
    CFG=nws/pickup/live/snap280_vel140_orbit2_mass85-gilman.cfg
fi
###############################################################
# collision between DM sub-halo (Gilman relation) and NW stream (t = 7.5 Gyr)
if [ $PROBLEM -eq 767 ]; then
    FILE=mass85-gilman-snap280_vel140_orbit3
    CFG=nws/pickup/live/snap280_vel140_orbit3_mass85-gilman.cfg
fi
###############################################################
# collision between DM sub-halo (Gilman relation) and NW stream (t = 7.5 Gyr)
if [ $PROBLEM -eq 770 ]; then
    FILE=mass90-gilman-snap280_vel140_orbit0
    CFG=nws/pickup/live/snap280_vel140_orbit0_mass90-gilman.cfg
fi
###############################################################
# collision between DM sub-halo (Gilman relation) and NW stream (t = 7.5 Gyr)
if [ $PROBLEM -eq 771 ]; then
    FILE=mass90-gilman-snap280_vel140_orbit1
    CFG=nws/pickup/live/snap280_vel140_orbit1_mass90-gilman.cfg
fi
###############################################################
# collision between DM sub-halo (Gilman relation) and NW stream (t = 7.5 Gyr)
if [ $PROBLEM -eq 772 ]; then
    FILE=mass90-gilman-snap280_vel140_orbit2
    CFG=nws/pickup/live/snap280_vel140_orbit2_mass90-gilman.cfg
fi
###############################################################
# collision between DM sub-halo (Gilman relation) and NW stream (t = 7.5 Gyr)
if [ $PROBLEM -eq 773 ]; then
    FILE=mass90-gilman-snap280_vel140_orbit3
    CFG=nws/pickup/live/snap280_vel140_orbit3_mass90-gilman.cfg
fi
###############################################################
# collision between DM sub-halo (Gilman relation) and NW stream (t = 7.5 Gyr)
if [ $PROBLEM -eq 774 ]; then
    FILE=mass95-gilman-snap280_vel140_orbit0
    CFG=nws/pickup/live/snap280_vel140_orbit0_mass95-gilman.cfg
fi
###############################################################
# collision between DM sub-halo (Gilman relation) and NW stream (t = 7.5 Gyr)
if [ $PROBLEM -eq 775 ]; then
    FILE=mass95-gilman-snap280_vel140_orbit1
    CFG=nws/pickup/live/snap280_vel140_orbit1_mass95-gilman.cfg
fi
###############################################################
# collision between DM sub-halo (Gilman relation) and NW stream (t = 7.5 Gyr)
if [ $PROBLEM -eq 776 ]; then
    FILE=mass95-gilman-snap280_vel140_orbit2
    CFG=nws/pickup/live/snap280_vel140_orbit2_mass95-gilman.cfg
fi
###############################################################
# collision between DM sub-halo (Gilman relation) and NW stream (t = 7.5 Gyr)
if [ $PROBLEM -eq 777 ]; then
    FILE=mass95-gilman-snap280_vel140_orbit3
    CFG=nws/pickup/live/snap280_vel140_orbit3_mass95-gilman.cfg
fi


###############################################################
# collision between DM sub-halo (test-particle) and NW stream (t = 8.0 Gyr)
if [ $PROBLEM -eq 800 ]; then
    FILE=mass80-test-snap300_vel140_orbit0
    CFG=nws/pickup/test/snap300_vel140_orbit0_mass80-test.cfg
fi
###############################################################
# collision between DM sub-halo (test-particle) and NW stream (t = 8.0 Gyr)
if [ $PROBLEM -eq 801 ]; then
    FILE=mass80-test-snap300_vel140_orbit1
    CFG=nws/pickup/test/snap300_vel140_orbit1_mass80-test.cfg
fi
###############################################################
# collision between DM sub-halo (test-particle) and NW stream (t = 8.0 Gyr)
if [ $PROBLEM -eq 802 ]; then
    FILE=mass80-test-snap300_vel140_orbit2
    CFG=nws/pickup/test/snap300_vel140_orbit2_mass80-test.cfg
fi
###############################################################
# collision between DM sub-halo (test-particle) and NW stream (t = 8.0 Gyr)
if [ $PROBLEM -eq 803 ]; then
    FILE=mass80-test-snap300_vel140_orbit3
    CFG=nws/pickup/test/snap300_vel140_orbit3_mass80-test.cfg
fi
###############################################################
# collision between DM sub-halo (test-particle) and NW stream (t = 8.0 Gyr)
if [ $PROBLEM -eq 804 ]; then
    FILE=mass85-test-snap300_vel140_orbit0
    CFG=nws/pickup/test/snap300_vel140_orbit0_mass85-test.cfg
fi
###############################################################
# collision between DM sub-halo (test-particle) and NW stream (t = 8.0 Gyr)
if [ $PROBLEM -eq 805 ]; then
    FILE=mass85-test-snap300_vel140_orbit1
    CFG=nws/pickup/test/snap300_vel140_orbit1_mass85-test.cfg
fi
###############################################################
# collision between DM sub-halo (test-particle) and NW stream (t = 8.0 Gyr)
if [ $PROBLEM -eq 806 ]; then
    FILE=mass85-test-snap300_vel140_orbit2
    CFG=nws/pickup/test/snap300_vel140_orbit2_mass85-test.cfg
fi
###############################################################
# collision between DM sub-halo (test-particle) and NW stream (t = 8.0 Gyr)
if [ $PROBLEM -eq 807 ]; then
    FILE=mass85-test-snap300_vel140_orbit3
    CFG=nws/pickup/test/snap300_vel140_orbit3_mass85-test.cfg
fi
###############################################################
# collision between DM sub-halo (test-particle) and NW stream (t = 8.0 Gyr)
if [ $PROBLEM -eq 810 ]; then
    FILE=mass90-test-snap300_vel140_orbit0
    CFG=nws/pickup/test/snap300_vel140_orbit0_mass90-test.cfg
fi
###############################################################
# collision between DM sub-halo (test-particle) and NW stream (t = 8.0 Gyr)
if [ $PROBLEM -eq 811 ]; then
    FILE=mass90-test-snap300_vel140_orbit1
    CFG=nws/pickup/test/snap300_vel140_orbit1_mass90-test.cfg
fi
###############################################################
# collision between DM sub-halo (test-particle) and NW stream (t = 8.0 Gyr)
if [ $PROBLEM -eq 812 ]; then
    FILE=mass90-test-snap300_vel140_orbit2
    CFG=nws/pickup/test/snap300_vel140_orbit2_mass90-test.cfg
fi
###############################################################
# collision between DM sub-halo (test-particle) and NW stream (t = 8.0 Gyr)
if [ $PROBLEM -eq 813 ]; then
    FILE=mass90-test-snap300_vel140_orbit3
    CFG=nws/pickup/test/snap300_vel140_orbit3_mass90-test.cfg
fi
###############################################################
# collision between DM sub-halo (test-particle) and NW stream (t = 8.0 Gyr)
if [ $PROBLEM -eq 814 ]; then
    FILE=mass95-test-snap300_vel140_orbit0
    CFG=nws/pickup/test/snap300_vel140_orbit0_mass95-test.cfg
fi
###############################################################
# collision between DM sub-halo (test-particle) and NW stream (t = 8.0 Gyr)
if [ $PROBLEM -eq 815 ]; then
    FILE=mass95-test-snap300_vel140_orbit1
    CFG=nws/pickup/test/snap300_vel140_orbit1_mass95-test.cfg
fi
###############################################################
# collision between DM sub-halo (test-particle) and NW stream (t = 8.0 Gyr)
if [ $PROBLEM -eq 816 ]; then
    FILE=mass95-test-snap300_vel140_orbit2
    CFG=nws/pickup/test/snap300_vel140_orbit2_mass95-test.cfg
fi
###############################################################
# collision between DM sub-halo (test-particle) and NW stream (t = 8.0 Gyr)
if [ $PROBLEM -eq 817 ]; then
    FILE=mass95-test-snap300_vel140_orbit3
    CFG=nws/pickup/test/snap300_vel140_orbit3_mass95-test.cfg
fi
###############################################################
# collision between DM sub-halo (Prada relation) and NW stream (t = 8.0 Gyr)
if [ $PROBLEM -eq 820 ]; then
    FILE=mass80-prada-snap300_vel140_orbit0
    CFG=nws/pickup/live/snap300_vel140_orbit0_mass80-prada.cfg
fi
###############################################################
# collision between DM sub-halo (Prada relation) and NW stream (t = 8.0 Gyr)
if [ $PROBLEM -eq 821 ]; then
    FILE=mass80-prada-snap300_vel140_orbit1
    CFG=nws/pickup/live/snap300_vel140_orbit1_mass80-prada.cfg
fi
###############################################################
# collision between DM sub-halo (Prada relation) and NW stream (t = 8.0 Gyr)
if [ $PROBLEM -eq 822 ]; then
    FILE=mass80-prada-snap300_vel140_orbit2
    CFG=nws/pickup/live/snap300_vel140_orbit2_mass80-prada.cfg
fi
###############################################################
# collision between DM sub-halo (Prada relation) and NW stream (t = 8.0 Gyr)
if [ $PROBLEM -eq 823 ]; then
    FILE=mass80-prada-snap300_vel140_orbit3
    CFG=nws/pickup/live/snap300_vel140_orbit3_mass80-prada.cfg
fi
###############################################################
# collision between DM sub-halo (Prada relation) and NW stream (t = 8.0 Gyr)
if [ $PROBLEM -eq 824 ]; then
    FILE=mass85-prada-snap300_vel140_orbit0
    CFG=nws/pickup/live/snap300_vel140_orbit0_mass85-prada.cfg
fi
###############################################################
# collision between DM sub-halo (Prada relation) and NW stream (t = 8.0 Gyr)
if [ $PROBLEM -eq 825 ]; then
    FILE=mass85-prada-snap300_vel140_orbit1
    CFG=nws/pickup/live/snap300_vel140_orbit1_mass85-prada.cfg
fi
###############################################################
# collision between DM sub-halo (Prada relation) and NW stream (t = 8.0 Gyr)
if [ $PROBLEM -eq 826 ]; then
    FILE=mass85-prada-snap300_vel140_orbit2
    CFG=nws/pickup/live/snap300_vel140_orbit2_mass85-prada.cfg
fi
###############################################################
# collision between DM sub-halo (Prada relation) and NW stream (t = 8.0 Gyr)
if [ $PROBLEM -eq 827 ]; then
    FILE=mass85-prada-snap300_vel140_orbit3
    CFG=nws/pickup/live/snap300_vel140_orbit3_mass85-prada.cfg
fi
###############################################################
# collision between DM sub-halo (Prada relation) and NW stream (t = 8.0 Gyr)
if [ $PROBLEM -eq 830 ]; then
    FILE=mass90-prada-snap300_vel140_orbit0
    CFG=nws/pickup/live/snap300_vel140_orbit0_mass90-prada.cfg
fi
###############################################################
# collision between DM sub-halo (Prada relation) and NW stream (t = 8.0 Gyr)
if [ $PROBLEM -eq 831 ]; then
    FILE=mass90-prada-snap300_vel140_orbit1
    CFG=nws/pickup/live/snap300_vel140_orbit1_mass90-prada.cfg
fi
###############################################################
# collision between DM sub-halo (Prada relation) and NW stream (t = 8.0 Gyr)
if [ $PROBLEM -eq 832 ]; then
    FILE=mass90-prada-snap300_vel140_orbit2
    CFG=nws/pickup/live/snap300_vel140_orbit2_mass90-prada.cfg
fi
###############################################################
# collision between DM sub-halo (Prada relation) and NW stream (t = 8.0 Gyr)
if [ $PROBLEM -eq 833 ]; then
    FILE=mass90-prada-snap300_vel140_orbit3
    CFG=nws/pickup/live/snap300_vel140_orbit3_mass90-prada.cfg
fi
###############################################################
# collision between DM sub-halo (Prada relation) and NW stream (t = 8.0 Gyr)
if [ $PROBLEM -eq 834 ]; then
    FILE=mass95-prada-snap300_vel140_orbit0
    CFG=nws/pickup/live/snap300_vel140_orbit0_mass95-prada.cfg
fi
###############################################################
# collision between DM sub-halo (Prada relation) and NW stream (t = 8.0 Gyr)
if [ $PROBLEM -eq 835 ]; then
    FILE=mass95-prada-snap300_vel140_orbit1
    CFG=nws/pickup/live/snap300_vel140_orbit1_mass95-prada.cfg
fi
###############################################################
# collision between DM sub-halo (Prada relation) and NW stream (t = 8.0 Gyr)
if [ $PROBLEM -eq 836 ]; then
    FILE=mass95-prada-snap300_vel140_orbit2
    CFG=nws/pickup/live/snap300_vel140_orbit2_mass95-prada.cfg
fi
###############################################################
# collision between DM sub-halo (Prada relation) and NW stream (t = 8.0 Gyr)
if [ $PROBLEM -eq 837 ]; then
    FILE=mass95-prada-snap300_vel140_orbit3
    CFG=nws/pickup/live/snap300_vel140_orbit3_mass95-prada.cfg
fi
###############################################################
# collision between DM sub-halo (Ishiyama & Ando relation) and NW stream (t = 8.0 Gyr)
if [ $PROBLEM -eq 840 ]; then
    FILE=mass80-ishiyama-snap300_vel140_orbit0
    CFG=nws/pickup/live/snap300_vel140_orbit0_mass80-ishiyama.cfg
fi
###############################################################
# collision between DM sub-halo (Ishiyama & Ando relation) and NW stream (t = 8.0 Gyr)
if [ $PROBLEM -eq 841 ]; then
    FILE=mass80-ishiyama-snap300_vel140_orbit1
    CFG=nws/pickup/live/snap300_vel140_orbit1_mass80-ishiyama.cfg
fi
###############################################################
# collision between DM sub-halo (Ishiyama & Ando relation) and NW stream (t = 8.0 Gyr)
if [ $PROBLEM -eq 842 ]; then
    FILE=mass80-ishiyama-snap300_vel140_orbit2
    CFG=nws/pickup/live/snap300_vel140_orbit2_mass80-ishiyama.cfg
fi
###############################################################
# collision between DM sub-halo (Ishiyama & Ando relation) and NW stream (t = 8.0 Gyr)
if [ $PROBLEM -eq 843 ]; then
    FILE=mass80-ishiyama-snap300_vel140_orbit3
    CFG=nws/pickup/live/snap300_vel140_orbit3_mass80-ishiyama.cfg
fi
###############################################################
# collision between DM sub-halo (Ishiyama & Ando relation) and NW stream (t = 8.0 Gyr)
if [ $PROBLEM -eq 844 ]; then
    FILE=mass85-ishiyama-snap300_vel140_orbit0
    CFG=nws/pickup/live/snap300_vel140_orbit0_mass85-ishiyama.cfg
fi
###############################################################
# collision between DM sub-halo (Ishiyama & Ando relation) and NW stream (t = 8.0 Gyr)
if [ $PROBLEM -eq 845 ]; then
    FILE=mass85-ishiyama-snap300_vel140_orbit1
    CFG=nws/pickup/live/snap300_vel140_orbit1_mass85-ishiyama.cfg
fi
###############################################################
# collision between DM sub-halo (Ishiyama & Ando relation) and NW stream (t = 8.0 Gyr)
if [ $PROBLEM -eq 846 ]; then
    FILE=mass85-ishiyama-snap300_vel140_orbit2
    CFG=nws/pickup/live/snap300_vel140_orbit2_mass85-ishiyama.cfg
fi
###############################################################
# collision between DM sub-halo (Ishiyama & Ando relation) and NW stream (t = 8.0 Gyr)
if [ $PROBLEM -eq 847 ]; then
    FILE=mass85-ishiyama-snap300_vel140_orbit3
    CFG=nws/pickup/live/snap300_vel140_orbit3_mass85-ishiyama.cfg
fi
###############################################################
# collision between DM sub-halo (Ishiyama & Ando relation) and NW stream (t = 8.0 Gyr)
if [ $PROBLEM -eq 850 ]; then
    FILE=mass90-ishiyama-snap300_vel140_orbit0
    CFG=nws/pickup/live/snap300_vel140_orbit0_mass90-ishiyama.cfg
fi
###############################################################
# collision between DM sub-halo (Ishiyama & Ando relation) and NW stream (t = 8.0 Gyr)
if [ $PROBLEM -eq 851 ]; then
    FILE=mass90-ishiyama-snap300_vel140_orbit1
    CFG=nws/pickup/live/snap300_vel140_orbit1_mass90-ishiyama.cfg
fi
###############################################################
# collision between DM sub-halo (Ishiyama & Ando relation) and NW stream (t = 8.0 Gyr)
if [ $PROBLEM -eq 852 ]; then
    FILE=mass90-ishiyama-snap300_vel140_orbit2
    CFG=nws/pickup/live/snap300_vel140_orbit2_mass90-ishiyama.cfg
fi
###############################################################
# collision between DM sub-halo (Ishiyama & Ando relation) and NW stream (t = 8.0 Gyr)
if [ $PROBLEM -eq 853 ]; then
    FILE=mass90-ishiyama-snap300_vel140_orbit3
    CFG=nws/pickup/live/snap300_vel140_orbit3_mass90-ishiyama.cfg
fi
###############################################################
# collision between DM sub-halo (Ishiyama & Ando relation) and NW stream (t = 8.0 Gyr)
if [ $PROBLEM -eq 854 ]; then
    FILE=mass95-ishiyama-snap300_vel140_orbit0
    CFG=nws/pickup/live/snap300_vel140_orbit0_mass95-ishiyama.cfg
fi
###############################################################
# collision between DM sub-halo (Ishiyama & Ando relation) and NW stream (t = 8.0 Gyr)
if [ $PROBLEM -eq 855 ]; then
    FILE=mass95-ishiyama-snap300_vel140_orbit1
    CFG=nws/pickup/live/snap300_vel140_orbit1_mass95-ishiyama.cfg
fi
###############################################################
# collision between DM sub-halo (Ishiyama & Ando relation) and NW stream (t = 8.0 Gyr)
if [ $PROBLEM -eq 856 ]; then
    FILE=mass95-ishiyama-snap300_vel140_orbit2
    CFG=nws/pickup/live/snap300_vel140_orbit2_mass95-ishiyama.cfg
fi
###############################################################
# collision between DM sub-halo (Ishiyama & Ando relation) and NW stream (t = 8.0 Gyr)
if [ $PROBLEM -eq 857 ]; then
    FILE=mass95-ishiyama-snap300_vel140_orbit3
    CFG=nws/pickup/live/snap300_vel140_orbit3_mass95-ishiyama.cfg
fi
###############################################################
# collision between DM sub-halo (Gilman relation) and NW stream (t = 8.0 Gyr)
if [ $PROBLEM -eq 860 ]; then
    FILE=mass80-gilman-snap300_vel140_orbit0
    CFG=nws/pickup/live/snap300_vel140_orbit0_mass80-gilman.cfg
fi
###############################################################
# collision between DM sub-halo (Gilman relation) and NW stream (t = 8.0 Gyr)
if [ $PROBLEM -eq 861 ]; then
    FILE=mass80-gilman-snap300_vel140_orbit1
    CFG=nws/pickup/live/snap300_vel140_orbit1_mass80-gilman.cfg
fi
###############################################################
# collision between DM sub-halo (Gilman relation) and NW stream (t = 8.0 Gyr)
if [ $PROBLEM -eq 862 ]; then
    FILE=mass80-gilman-snap300_vel140_orbit2
    CFG=nws/pickup/live/snap300_vel140_orbit2_mass80-gilman.cfg
fi
###############################################################
# collision between DM sub-halo (Gilman relation) and NW stream (t = 8.0 Gyr)
if [ $PROBLEM -eq 863 ]; then
    FILE=mass80-gilman-snap300_vel140_orbit3
    CFG=nws/pickup/live/snap300_vel140_orbit3_mass80-gilman.cfg
fi
###############################################################
# collision between DM sub-halo (Gilman relation) and NW stream (t = 8.0 Gyr)
if [ $PROBLEM -eq 864 ]; then
    FILE=mass85-gilman-snap300_vel140_orbit0
    CFG=nws/pickup/live/snap300_vel140_orbit0_mass85-gilman.cfg
fi
###############################################################
# collision between DM sub-halo (Gilman relation) and NW stream (t = 8.0 Gyr)
if [ $PROBLEM -eq 865 ]; then
    FILE=mass85-gilman-snap300_vel140_orbit1
    CFG=nws/pickup/live/snap300_vel140_orbit1_mass85-gilman.cfg
fi
###############################################################
# collision between DM sub-halo (Gilman relation) and NW stream (t = 8.0 Gyr)
if [ $PROBLEM -eq 866 ]; then
    FILE=mass85-gilman-snap300_vel140_orbit2
    CFG=nws/pickup/live/snap300_vel140_orbit2_mass85-gilman.cfg
fi
###############################################################
# collision between DM sub-halo (Gilman relation) and NW stream (t = 8.0 Gyr)
if [ $PROBLEM -eq 867 ]; then
    FILE=mass85-gilman-snap300_vel140_orbit3
    CFG=nws/pickup/live/snap300_vel140_orbit3_mass85-gilman.cfg
fi
###############################################################
# collision between DM sub-halo (Gilman relation) and NW stream (t = 8.0 Gyr)
if [ $PROBLEM -eq 870 ]; then
    FILE=mass90-gilman-snap300_vel140_orbit0
    CFG=nws/pickup/live/snap300_vel140_orbit0_mass90-gilman.cfg
fi
###############################################################
# collision between DM sub-halo (Gilman relation) and NW stream (t = 8.0 Gyr)
if [ $PROBLEM -eq 871 ]; then
    FILE=mass90-gilman-snap300_vel140_orbit1
    CFG=nws/pickup/live/snap300_vel140_orbit1_mass90-gilman.cfg
fi
###############################################################
# collision between DM sub-halo (Gilman relation) and NW stream (t = 8.0 Gyr)
if [ $PROBLEM -eq 872 ]; then
    FILE=mass90-gilman-snap300_vel140_orbit2
    CFG=nws/pickup/live/snap300_vel140_orbit2_mass90-gilman.cfg
fi
###############################################################
# collision between DM sub-halo (Gilman relation) and NW stream (t = 8.0 Gyr)
if [ $PROBLEM -eq 873 ]; then
    FILE=mass90-gilman-snap300_vel140_orbit3
    CFG=nws/pickup/live/snap300_vel140_orbit3_mass90-gilman.cfg
fi
###############################################################
# collision between DM sub-halo (Gilman relation) and NW stream (t = 8.0 Gyr)
if [ $PROBLEM -eq 874 ]; then
    FILE=mass95-gilman-snap300_vel140_orbit0
    CFG=nws/pickup/live/snap300_vel140_orbit0_mass95-gilman.cfg
fi
###############################################################
# collision between DM sub-halo (Gilman relation) and NW stream (t = 8.0 Gyr)
if [ $PROBLEM -eq 875 ]; then
    FILE=mass95-gilman-snap300_vel140_orbit1
    CFG=nws/pickup/live/snap300_vel140_orbit1_mass95-gilman.cfg
fi
###############################################################
# collision between DM sub-halo (Gilman relation) and NW stream (t = 8.0 Gyr)
if [ $PROBLEM -eq 876 ]; then
    FILE=mass95-gilman-snap300_vel140_orbit2
    CFG=nws/pickup/live/snap300_vel140_orbit2_mass95-gilman.cfg
fi
###############################################################
# collision between DM sub-halo (Gilman relation) and NW stream (t = 8.0 Gyr)
if [ $PROBLEM -eq 877 ]; then
    FILE=mass95-gilman-snap300_vel140_orbit3
    CFG=nws/pickup/live/snap300_vel140_orbit3_mass95-gilman.cfg
fi




###############################################################
# set input arguments
OPTION="-file=$FILE -list=$CFG -eps=$EPS -ft=$FINISH -eta=$ETA -snapshotInterval=$INTERVAL -saveInterval=$SAVE"
###############################################################


###############################################################
# job execution via UNIVA Grid Engine
###############################################################
# set stdout and stderr
STDOUT=log/${FILE}_$REQUEST.o$JOB_ID
STDERR=log/${FILE}_$REQUEST.e$JOB_ID
###############################################################
# load modules
. /etc/profile.d/modules.sh
export MODULEPATH=$MODULEPATH:/gs/hs1/jh180045/share/opt/Modules
module load intel/19.0.0.117 cuda/9.2.148 openmpi/2.1.2-opa10.9
module load phdf5
module list 1>>$STDOUT 2>>$STDERR
###############################################################
cat $PE_HOSTFILE 1>>$STDOUT 2>>$STDERR
TIME=`date`
echo "start: $TIME" 1>>$STDOUT 2>>$STDERR
###############################################################
# execute the job
if [ `which numactl` ]; then
    # run with numactl
    echo "numactl --localalloc $EXEC $OPTION 1>>$STDOUT 2>>$STDERR" 1>>$STDOUT 2>>$STDERR
    numactl --localalloc $EXEC $OPTION 1>>$STDOUT 2>>$STDERR
else
    # run without numactl
    echo "$EXEC $OPTION 1>>$STDOUT 2>>$STDERR" 1>>$STDOUT 2>>$STDERR
    $EXEC $OPTION 1>>$STDOUT 2>>$STDERR
fi
###############################################################
# finish logging
TIME=`date`
echo "finish: $TIME" 1>>$STDOUT 2>>$STDERR
###############################################################
exit 0
###############################################################
