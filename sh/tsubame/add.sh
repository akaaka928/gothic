#!/bin/sh
#########
qsub -g tge-17IJ0016 -N job06 -v ABSERR=1.5625000000000e-2                 sh/tsubame/run.sh
qsub -g tge-17IJ0016 -N job06 -v ABSERR=1.5625000000000e-2 -hold_jid job06 sh/tsubame/run.sh
qsub -g tge-17IJ0016 -N job06 -v ABSERR=1.5625000000000e-2 -hold_jid job06 sh/tsubame/run.sh
qsub -g tge-17IJ0016 -N job06 -v ABSERR=1.5625000000000e-2 -hold_jid job06 sh/tsubame/run.sh
qsub -g tge-17IJ0016 -N job06 -v ABSERR=1.5625000000000e-2 -hold_jid job06 sh/tsubame/run.sh
#########
qsub -g tge-17IJ0016 -N job05 -v ABSERR=3.1250000000000e-2 -hold_jid job06 sh/tsubame/run.sh
qsub -g tge-17IJ0016 -N job05 -v ABSERR=3.1250000000000e-2 -hold_jid job05 sh/tsubame/run.sh
qsub -g tge-17IJ0016 -N job05 -v ABSERR=3.1250000000000e-2 -hold_jid job05 sh/tsubame/run.sh
qsub -g tge-17IJ0016 -N job05 -v ABSERR=3.1250000000000e-2 -hold_jid job05 sh/tsubame/run.sh
qsub -g tge-17IJ0016 -N job05 -v ABSERR=3.1250000000000e-2 -hold_jid job05 sh/tsubame/run.sh
#########
qsub -g tge-17IJ0016 -N job04 -v ABSERR=6.2500000000000e-2 -hold_jid job05 sh/tsubame/run.sh
qsub -g tge-17IJ0016 -N job04 -v ABSERR=6.2500000000000e-2 -hold_jid job04 sh/tsubame/run.sh
qsub -g tge-17IJ0016 -N job04 -v ABSERR=6.2500000000000e-2 -hold_jid job04 sh/tsubame/run.sh
qsub -g tge-17IJ0016 -N job04 -v ABSERR=6.2500000000000e-2 -hold_jid job04 sh/tsubame/run.sh
qsub -g tge-17IJ0016 -N job04 -v ABSERR=6.2500000000000e-2 -hold_jid job04 sh/tsubame/run.sh
#########
qsub -g tge-17IJ0016 -N job03 -v ABSERR=1.2500000000000e-1 -hold_jid job04 sh/tsubame/run.sh
qsub -g tge-17IJ0016 -N job03 -v ABSERR=1.2500000000000e-1 -hold_jid job03 sh/tsubame/run.sh
qsub -g tge-17IJ0016 -N job03 -v ABSERR=1.2500000000000e-1 -hold_jid job03 sh/tsubame/run.sh
qsub -g tge-17IJ0016 -N job03 -v ABSERR=1.2500000000000e-1 -hold_jid job03 sh/tsubame/run.sh
qsub -g tge-17IJ0016 -N job03 -v ABSERR=1.2500000000000e-1 -hold_jid job03 sh/tsubame/run.sh
#########
qsub -g tge-17IJ0016 -N job01 -v ABSERR=5.0000000000000e-1 -hold_jid job03 sh/tsubame/run.sh
qsub -g tge-17IJ0016 -N job01 -v ABSERR=5.0000000000000e-1 -hold_jid job01 sh/tsubame/run.sh
qsub -g tge-17IJ0016 -N job01 -v ABSERR=5.0000000000000e-1 -hold_jid job01 sh/tsubame/run.sh
qsub -g tge-17IJ0016 -N job01 -v ABSERR=5.0000000000000e-1 -hold_jid job01 sh/tsubame/run.sh
qsub -g tge-17IJ0016 -N job01 -v ABSERR=5.0000000000000e-1 -hold_jid job01 sh/tsubame/run.sh
#########
exit 0
#########
