##############################################################
##############################################################
# GNUPLOT plot file
#     written by Yohei MIKI
#            last updated on 2015/10/16(Fri) 16:09:04
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
# set parameter sets of Hernquist sphere
##############################################################
Mtot=1.0
Ntot=131072
Nng=16
rs=1.0
##############################################################
rho(x)=Mtot/(2.0*pi*(rs**3)*(x/rs)*((1.0+(x/rs))**3))
enc(x)=Mtot*(((x/rs)/(1.0+(x/rs)))**2)
##############################################################
# x = M(r)/Mtot
rad(x)=(1.0+sqrt(2.0-x))/(1.0-x)
##############################################################

##############################################################
neighbor(x)=(3.0*Mtot*Nng/(4.0*pi*Ntot*rho(x)))**(1.0/3.0)
##############################################################

##############################################################
# set global settings
##############################################################
se term post eps enh col 15
se st da l
se si sq
se g
#se g lt -1 lw 0, lt 0 lw 0
se nokey
##############################################################
se xl 'radius'
se yl 'neighbor length'
##############################################################
se xr [0:16]
se out 'arm.eps'
p neighbor(x) ls 1
se out
se xr [*:*]
##############################################################
se param
se tr [0:Ntot-1]
se xl 'particle index'
##############################################################
#se format x "%.0t*10^%T"
se format x "%.1te%+-T"
se yr [*:10]
se out 'near.eps'
p t,neighbor(rad(t/Ntot)) ls 1
se out
##############################################################
se tr [1:Ntot-1]
#se xr [0:10.0*rs]
se log x
se format x "10^{%L}"
#se g x mx
##############################################################
se yr [*:*]
se log y
se format y "10^{%L}"
# se g y my
##############################################################
se out 'neighbor.eps'
p t,neighbor(rad(t/Ntot)) ls 1
se out
##############################################################
