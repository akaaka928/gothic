##############################################################
##############################################################
# GNUPLOT plot file
#     written by Yohei MIKI
#            last updated on 2015/01/20(Tue) 16:31:49
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
se key outside
##############################################################
se xl 'x'
se yl 'y'
se out filename.".xy-plane.eps"
p filename u 4:5:7:8   w vec ti 'tree',\
  filename u 4:5:11:12 w vec ti 'direct'
se out
se out filename.".xy-dist.eps"
p filename u 4:5 w d noti
se out
##############################################################
se xl 'y'
se yl 'z'
se out filename.".yz-plane.eps"
p filename u 5:6:8:9  w vec ti 'tree',\
  filename u 5:6:12:13 w vec ti 'direct'
se out
se out filename.".yz-dist.eps"
p filename u 5:6 w d noti
se out
##############################################################
se xl 'z'
se yl 'x'
se out filename.".zx-plane.eps"
p filename u 6:4:9:7   w vec ti 'tree',\
  filename u 6:4:13:11 w vec ti 'direct'
se out
se out filename.".zx-dist.eps"
p filename u 6:4 w d noti
se out
##############################################################
se xl 'r'
se yl 'pot'
se out filename.".pot.eps"
p filename u (($4**2+$5**2+$6**2)**0.5):10 w p ti 'tree',\
  filename u (($4**2+$5**2+$6**2)**0.5):14 w p ti 'direct'
se out
##############################################################
