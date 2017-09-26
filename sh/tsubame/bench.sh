#!/bin/sh
#########
qsub -g tge-17IJ0016 -N job01 -v ABSERR=5.0000000000000e-1 -hold_jid magi  sh/tsubame/run.sh
qsub -g tge-17IJ0016 -N job02 -v ABSERR=2.5000000000000e-1 -hold_jid job01 sh/tsubame/run.sh
qsub -g tge-17IJ0016 -N job03 -v ABSERR=1.2500000000000e-1 -hold_jid job02 sh/tsubame/run.sh
qsub -g tge-17IJ0016 -N job04 -v ABSERR=6.2500000000000e-2 -hold_jid job03 sh/tsubame/run.sh
qsub -g tge-17IJ0016 -N job05 -v ABSERR=3.1250000000000e-2 -hold_jid job04 sh/tsubame/run.sh
qsub -g tge-17IJ0016 -N job06 -v ABSERR=1.5625000000000e-2 -hold_jid job05 sh/tsubame/run.sh
qsub -g tge-17IJ0016 -N job07 -v ABSERR=7.8125000000000e-3 -hold_jid job06 sh/tsubame/run.sh
qsub -g tge-17IJ0016 -N job08 -v ABSERR=3.9062500000000e-3 -hold_jid job07 sh/tsubame/run.sh
qsub -g tge-17IJ0016 -N job09 -v ABSERR=1.9531250000000e-3 -hold_jid job08 sh/tsubame/run.sh
qsub -g tge-17IJ0016 -N job10 -v ABSERR=9.7656250000000e-4 -hold_jid job09 sh/tsubame/run.sh
qsub -g tge-17IJ0016 -N job11 -v ABSERR=4.8828125000000e-4 -hold_jid job10 sh/tsubame/run.sh
qsub -g tge-17IJ0016 -N job12 -v ABSERR=2.4414062500000e-4 -hold_jid job11 sh/tsubame/run.sh
qsub -g tge-17IJ0016 -N job13 -v ABSERR=1.2207031250000e-4 -hold_jid job12 sh/tsubame/run.sh
qsub -g tge-17IJ0016 -N job14 -v ABSERR=6.1035156250000e-5 -hold_jid job13 sh/tsubame/run.sh
qsub -g tge-17IJ0016 -N job15 -v ABSERR=3.0517578125000e-5 -hold_jid job14 sh/tsubame/run.sh
qsub -g tge-17IJ0016 -N job16 -v ABSERR=1.5258789062500e-5 -hold_jid job15 sh/tsubame/run.sh
qsub -g tge-17IJ0016 -N job17 -v ABSERR=7.6293945312500e-6 -hold_jid job16 sh/tsubame/run.sh
qsub -g tge-17IJ0016 -N job18 -v ABSERR=3.8146972656250e-6 -hold_jid job17 sh/tsubame/run.sh
qsub -g tge-17IJ0016 -N job19 -v ABSERR=1.9073486328125e-6 -hold_jid job18 sh/tsubame/run.sh
qsub -g tge-17IJ0016 -N job20 -v ABSERR=9.5367431640625e-7 -hold_jid job19 sh/tsubame/run.sh
#########
exit 0
#########
