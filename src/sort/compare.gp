##############################################################
##############################################################
# GNUPLOT plot file
#     written by Yohei MIKI
#            last updated on 2015/07/22(Wed) 20:50:18
##############################################################
##############################################################
# set line styles
##############################################################
# lt 1: full line
# lt 3: dotted line
# lt 0: dashed line
# lt 5: dot-dashed line
se st l  1 lt 1 lw 3 lc rgb "black"
se st l  2 lt 3 lw 3 lc rgb "black"
se st l  3 lt 0 lw 3 lc rgb "black"
se st l  4 lt 5 lw 3 lc rgb "black"
se st l  5 lt 1 lw 3 lc rgb "red"
se st l  6 lt 3 lw 3 lc rgb "red"
se st l  7 lt 0 lw 3 lc rgb "red"
se st l  8 lt 5 lw 3 lc rgb "red"
se st l  9 lt 1 lw 3 lc rgb "blue"
se st l 10 lt 3 lw 3 lc rgb "blue"
se st l 11 lt 0 lw 3 lc rgb "blue"
se st l 12 lt 5 lw 3 lc rgb "blue"
se st l 13 lt 1 lw 3 lc rgb "magenta"
se st l 14 lt 3 lw 3 lc rgb "magenta"
se st l 15 lt 0 lw 3 lc rgb "magenta"
se st l 16 lt 5 lw 3 lc rgb "magenta"
##############################################################
# color of points
# lt 7: black
# lt 1: red
# lt 3: blue
# lt 4: magenta
##############################################################


##############################################################
# set global settings
##############################################################
se term post eps enh col 25
# se si sq
se st da l
se g
se key bottom
##############################################################
se log x
se format x "10^{%L}"
se xl '# of keys'
##############################################################
se log y
se format y "10^{%L}"
se yl 'Sorting rate (keys / sec)'
##############################################################


##############################################################
se out "rate.eps"
p 'log/bench.merge.dat' u 1:3 ls 15 ti 'thrust::sort',\
  'log/bench.merge.dat' u 1:3 w p pt 2 lt 4 ps 2 noti,\
  'log/bench.merge.dat' u 1:5 ls 12 ti 'thrust::stable\_sort',\
  'log/bench.merge.dat' u 1:5 w p pt 6 lt 3 ps 2 noti,\
  'log/bench.split.dat' u 1:7 ls  2 ti 'original (2 kernels mode)',\
  'log/bench.split.dat' u 1:7 w p pt 1 lt 0 ps 2 noti,\
  'log/bench.merge.dat' u 1:7 ls  5 ti 'original (1 kernel mode)',\
  'log/bench.merge.dat' u 1:7 w p pt 4 lt 1 ps 2 noti
se out
##############################################################
se yl 'Elapsed time (sec)'
se out "time.eps"
p 'log/bench.merge.dat' u 1:2 ls 15 ti 'thrust::sort',\
  'log/bench.merge.dat' u 1:2 w p pt 2 lt 4 ps 2 noti,\
  'log/bench.merge.dat' u 1:4 ls 12 ti 'thrust::stable\_sort',\
  'log/bench.merge.dat' u 1:4 w p pt 6 lt 3 ps 2 noti,\
  'log/bench.split.dat' u 1:6 ls  2 ti 'original (2 kernels mode)',\
  'log/bench.split.dat' u 1:6 w p pt 1 lt 0 ps 2 noti,\
  'log/bench.merge.dat' u 1:6 ls  5 ti 'original (1 kernel mode)',\
  'log/bench.merge.dat' u 1:6 w p pt 4 lt 1 ps 2 noti
se out
##############################################################
