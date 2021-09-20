##############################################################
##############################################################
# GNUPLOT plot file
#     written by Yohei MIKI
#            last updated on 2015/01/15(Thu) 10:12:09
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
se st da p
se si sq
se g
#se nokey
##############################################################
se xr [-1:1]
se yr [-1:1]
##############################################################
se xl 'x'
se yl 'y'
se out filename.".xy-plane.eps"
p filename ev :::0::1 u 3:4 ti 'ful',\
  filename ev :::1::2 u 3:4 ti 'let'
se out
##############################################################
se xl 'y'
se yl 'z'
se out filename.".yz-plane.eps"
p filename ev :::0::1 u 4:5 ti 'ful',\
  filename ev :::1::2 u 4:5 ti 'let'
se out
##############################################################
se xl 'z'
se yl 'x'
se out filename.".zx-plane.eps"
p filename ev :::0::1 u 6:3 ti 'ful',\
  filename ev :::1::2 u 6:3 ti 'let'
se out
##############################################################
se xr [1.0e-4:1]
se log x
##############################################################
se xl 'm'
se yl 'x'
se out filename.".mx-plane.eps"
p filename ev :::0::1 u 2:3 ti 'ful',\
  filename ev :::1::2 u 2:3 ti 'let'
se out
##############################################################
se xl 'm'
se yl 'y'
se out filename.".my-plane.eps"
p filename ev :::0::1 u 2:4 ti 'ful',\
  filename ev :::1::2 u 2:4 ti 'let'
se out
##############################################################
se xl 'm'
se yl 'z'
se out filename.".mz-plane.eps"
p filename ev :::0::1 u 2:5 ti 'ful',\
  filename ev :::1::2 u 2:5 ti 'let'
se out
##############################################################
