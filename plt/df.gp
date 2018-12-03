##############################################################
##############################################################
# GNUPLOT plot file
#     written by Yohei MIKI
#            last updated on 2018/11/22 (Thu) 15:28:28
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
# se term post eps enh col 25
# se st da p
se st da l
se si sq
se g
#se nokey
##############################################################
se log y
##############################################################

##############################################################
Mtot=1.0
rs=1.0
G=1.0
ra=0.5
vg=sqrt(G*Mtot/rs)
##############################################################
q(x)=sqrt(x*rs/(G*Mtot))
org(x)=Mtot*(3*asin(q(x))+q(x)*sqrt(1-q(x)*q(x))*(1-2*q(x)*q(x))*(8*q(x)**4-8*q(x)**2-3))/(((1-q(x)*q(x))**2.5)*8.0*sqrt(2)*pi*pi*pi*rs*rs*rs*vg*vg*vg)
add(x)=Mtot*q(x)*(1.0-2.0*q(x)*q(x))*(rs*rs/(ra*ra))/(sqrt(2.0)*pi*pi*pi*rs*rs*rs*vg*vg*vg)
fin(x)=org(x)+add(x)
##############################################################
# p 'log/hernquist_magi.e2270' u 1:2,\
#  fin(x)
##############################################################
p 'log/hernquist_magi.e2270' u ($1):(abs($2/fin($1)-1.0))
##############################################################
