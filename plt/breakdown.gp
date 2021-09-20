##############################################################
##############################################################
# GNUPLOT plot file
#     written by Yohei MIKI
#            last updated on 2016/03/09(Wed) 16:52:42
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
# se si sq
# se si 1.5,1
se si 2,1
se g
# se nokey
se key outside
# se key bottom
##############################################################
se xr [0:1024]
# se log x
# se format x "10^{%L}"
se xl 'Step'
##############################################################
#se yr [1.0e-2:1.0e+3]
se log y
se format y "10^{%L}"
se yl 'Time (s)'
##############################################################
se out filename.".breakdown_all.eps"
if( block_time_step == 1 ) \
p filename u 1:2  w p pt 65 lt 1 ps 1 ti 'walk tree',\
  filename u 1:4  w p pt  7 lt 4 ps 1 ti 'make tree',\
  filename u 1:6  w p pt 13 lt 3 ps 1 ti 'PH-key',\
  filename u 1:3  w p pt  1 lt 7 ps 1 ti 'calc MAC',\
  filename u 1:16 w p pt  5 lt 3 ps 1 ti 'split groups',\
  filename u 1:11 w p pt  1 lt 4 ps 1 ti 'predict',\
  filename u 1:12 w p pt 66 lt 7 ps 1 ti 'correct';\
else \
p filename u 1:2  w p pt  5 lt 1 ps 1 ti 'walk tree',\
  filename u 1:4  w p pt 64 lt 7 ps 1 ti 'make tree',\
  filename u 1:3  w p pt  1 lt 7 ps 1 ti 'calc MAC',\
  filename u 1:6  w p pt  7 lt 3 ps 1 ti 'PH-key',\
  filename u 1:13 w p pt  7 lt 3 ps 1 ti 'calc separation',\
  filename u 1:14 w p pt  7 lt 3 ps 1 ti 'split groups',\
  filename u 1:11 w p pt  2 lt 3 ps 1 ti 'adv pos',\
  filename u 1:12 w p pt  2 lt 4 ps 1 ti 'adv vel',\
  filename u 1:5  w p pt  0 lt 1 ps 1 ti 'set dt',\
  filename u 1:7  w p pt 68 lt 4 ps 1 ti 'copy body (d2h)',\
  filename u 1:8  w p pt 68 lt 3 ps 1 ti 'copy body (h2d)',\
  filename u 1:9  w p pt 66 lt 7 ps 1 ti 'copy node (h2d)',\
  filename u 1:10 w p pt 67 lt 7 ps 1 ti 'copy cell (h2d)';\
se out
##############################################################
se si 2,1
##############################################################
se yr [1.0e-3:*]
se out filename.".breakdown.eps"
p filename u 1:2  w p pt 65 lt 1 ps 1 ti 'walk tree',\
  filename u 1:4  w p pt  7 lt 4 ps 1 ti 'make tree',\
  filename u 1:6  w p pt 13 lt 3 ps 1 ti 'PH-key',\
  filename u 1:3  w p pt  1 lt 7 ps 1 ti 'calc MAC',\
  filename u 1:16 w p pt  5 lt 3 ps 1 ti 'split groups',\
  filename u 1:11 w p pt  1 lt 4 ps 1 ti 'predict',\
  filename u 1:12 w p pt 66 lt 7 ps 1 ti 'correct'
se out
##############################################################
se xr [0:255]
se out filename.".breakdown_0_255.eps"
p filename u 1:2  w p pt 65 lt 1 ps 1 ti 'walk tree',\
  filename u 1:4  w p pt  7 lt 4 ps 1 ti 'make tree',\
  filename u 1:6  w p pt 13 lt 3 ps 1 ti 'PH-key',\
  filename u 1:3  w p pt  1 lt 7 ps 1 ti 'calc MAC',\
  filename u 1:16 w p pt  5 lt 3 ps 1 ti 'split groups',\
  filename u 1:11 w p pt  1 lt 4 ps 1 ti 'predict',\
  filename u 1:12 w p pt 66 lt 7 ps 1 ti 'correct'
se out
##############################################################
se xr [256:511]
se out filename.".breakdown_256_511.eps"
p filename u 1:2  w p pt 65 lt 1 ps 1 ti 'walk tree',\
  filename u 1:4  w p pt  7 lt 4 ps 1 ti 'make tree',\
  filename u 1:6  w p pt 13 lt 3 ps 1 ti 'PH-key',\
  filename u 1:3  w p pt  1 lt 7 ps 1 ti 'calc MAC',\
  filename u 1:16 w p pt  5 lt 3 ps 1 ti 'split groups',\
  filename u 1:11 w p pt  1 lt 4 ps 1 ti 'predict',\
  filename u 1:12 w p pt 66 lt 7 ps 1 ti 'correct'
se out
##############################################################
se xr [512:767]
se out filename.".breakdown_512_767.eps"
p filename u 1:2  w p pt 65 lt 1 ps 1 ti 'walk tree',\
  filename u 1:4  w p pt  7 lt 4 ps 1 ti 'make tree',\
  filename u 1:6  w p pt 13 lt 3 ps 1 ti 'PH-key',\
  filename u 1:3  w p pt  1 lt 7 ps 1 ti 'calc MAC',\
  filename u 1:16 w p pt  5 lt 3 ps 1 ti 'split groups',\
  filename u 1:11 w p pt  1 lt 4 ps 1 ti 'predict',\
  filename u 1:12 w p pt 66 lt 7 ps 1 ti 'correct'
se out
##############################################################
se xr [768:1023]
se out filename.".breakdown_768_1023.eps"
p filename u 1:2  w p pt 65 lt 1 ps 1 ti 'walk tree',\
  filename u 1:4  w p pt  7 lt 4 ps 1 ti 'make tree',\
  filename u 1:6  w p pt 13 lt 3 ps 1 ti 'PH-key',\
  filename u 1:3  w p pt  1 lt 7 ps 1 ti 'calc MAC',\
  filename u 1:16 w p pt  5 lt 3 ps 1 ti 'split groups',\
  filename u 1:11 w p pt  1 lt 4 ps 1 ti 'predict',\
  filename u 1:12 w p pt 66 lt 7 ps 1 ti 'correct'
se out
##############################################################
se xr [0:511]
se out filename.".breakdown_0_511.eps"
p filename u 1:2  w p pt 65 lt 1 ps 1 ti 'walk tree',\
  filename u 1:4  w p pt  7 lt 4 ps 1 ti 'make tree',\
  filename u 1:6  w p pt 13 lt 3 ps 1 ti 'PH-key',\
  filename u 1:3  w p pt  1 lt 7 ps 1 ti 'calc MAC',\
  filename u 1:16 w p pt  5 lt 3 ps 1 ti 'split groups',\
  filename u 1:11 w p pt  1 lt 4 ps 1 ti 'predict',\
  filename u 1:12 w p pt 66 lt 7 ps 1 ti 'correct'
se out
##############################################################
se xr [512:1023]
se out filename.".breakdown_512_1023.eps"
p filename u 1:2  w p pt 65 lt 1 ps 1 ti 'walk tree',\
  filename u 1:4  w p pt  7 lt 4 ps 1 ti 'make tree',\
  filename u 1:6  w p pt 13 lt 3 ps 1 ti 'PH-key',\
  filename u 1:3  w p pt  1 lt 7 ps 1 ti 'calc MAC',\
  filename u 1:16 w p pt  5 lt 3 ps 1 ti 'split groups',\
  filename u 1:11 w p pt  1 lt 4 ps 1 ti 'predict',\
  filename u 1:12 w p pt 66 lt 7 ps 1 ti 'correct'
se out
##############################################################
se xr [0:300]
se out filename.".breakdown_300.eps"
p filename u 1:2  w p pt 65 lt 1 ps 1 ti 'walk tree',\
  filename u 1:4  w p pt  7 lt 4 ps 1 ti 'make tree',\
  filename u 1:6  w p pt 13 lt 3 ps 1 ti 'PH-key',\
  filename u 1:3  w p pt  1 lt 7 ps 1 ti 'calc MAC',\
  filename u 1:16 w p pt  5 lt 3 ps 1 ti 'split groups',\
  filename u 1:11 w p pt  1 lt 4 ps 1 ti 'predict',\
  filename u 1:12 w p pt 66 lt 7 ps 1 ti 'correct'
se out
##############################################################
