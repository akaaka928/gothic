##############################################################
##############################################################
# GNUPLOT plot file
#     written by Yohei MIKI
#            last updated on 2015/02/17(Tue) 12:05:43
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
se st da l
se si sq
se g
se nokey
# se key outside
# se key bottom
##############################################################
# se xr [1.0e-3:1.0e+1]
# se log x
# se format x "10^{%L}"
# se xl 'B/F'
##############################################################
# se yr [20:3000]
# se log y
# se format y "10^{%L}"
# se yl 'performance (GFlop/s)'
##############################################################
se xl 'x'
se yl 'y'
# se out 'xy-plane.eps'
se out filename.".xy-plane.eps"
p filename ev :::0::1 u 2:3:8:9 w vec,\
  filename ev :::1::2 u 2:3:8:9 w vec,\
  filename ev :::0::1 u 2:3:5:6 w vec,\
  filename ev :::1::2 u 2:3:5:6 w vec
se out
##############################################################
