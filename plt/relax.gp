##############################################################
##############################################################
# GNUPLOT plot file
#     written by Yohei MIKI
#            last updated on 2015/09/29(Tue) 14:04:03
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
# set astrophysical constants (in CGS units)
##############################################################
newton=6.67428e-08
Gyr=3.15576e+16
Msun=1.9884e+33
pc=3.0856775975e+18
##############################################################

##############################################################
# set parameter sets of Hernquist sphere
##############################################################
Mtot=3.24e+10*Msun
rs=610.0*pc
eps=15.625*pc
##############################################################
rho(x)=Mtot/(2.0*pi*(rs**3)*(x/rs)*((1.0+(x/rs))**3))
enc(x)=Mtot*(((x/rs)/(1.0+(x/rs)))**2)
##############################################################

##############################################################
# set equation to calculate relaxation time by 2-body relaxation
# unit of x-axis:  pc
# unit of y-axis: Gyr
##############################################################
relax(x)=(enc(x*pc)/Mtot)*sqrt(3.0/(2.0*pi*newton*rho(x*pc)))/(32.0*log(x*pc/eps))/Gyr
##############################################################

##############################################################
# set global settings
##############################################################
se term post eps enh col 25
se st da l
#se si sq
#se g
se g lt -1 lw 0, lt 0 lw 0
se key outside
##############################################################
se xr [0.1*rs/pc:10.0*rs/pc]
se log x
se format x "10^{%L}"
se xl 'r (pc)'
se g x mx
##############################################################
se log y
se format y "10^{%L}"
se yl 't_{relax} (Gyr)'
# se g y my
##############################################################
# se param
# se tr [*:*]
# const=rs/pc
se out 'relax.eps'
p relax(x)*65536   ls  1 ti 'N =  64k',\
  relax(x)*131072  ls  6 ti 'N = 128k',\
  relax(x)*262144  ls 11 ti 'N = 256k',\
  relax(x)*524288  ls 16 ti 'N = 512k',\
  relax(x)*1048576 ls  2 ti 'N =   1M',\
  relax(x)*2097152 ls  7 ti 'N =   2M',\
  relax(x)*4194304 ls 12 ti 'N =   4M',\
  relax(x)*8388608 ls 13 ti 'N =   8M'
# rep const,t
se out
##############################################################
