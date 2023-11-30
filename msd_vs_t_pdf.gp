set terminal epslatex size 2.2,2.2 standalone color colortext 8
set output 'msd_vs_T.tex'

# set multiplot

# set lmargin at screen 0.15
# set rmargin at screen 1.0
# set tmargin at screen 0.9
# set bmargin at screen 0.2

set size square

unset key

set xtics offset 0,-0.3
set xtics ( '\fontsize{8}{60}\selectfont$0.0$' 0.0, '\fontsize{8}{60}\selectfont$5.0$' 50000, '\fontsize{8}{60}\selectfont$10.0$' 100000.0)
set xlabel '\fontsize{10}{60}\selectfont$t/\tau\,(10^4)$'
set xlabel offset 0.0,-0.5

# set format y "%.1f"


set ylabel '\fontsize{10}{60}\selectfont MSD (nm)'
set ylabel offset 0.4,0.0


plot 'len_50/seed_14/dt_0.01/T_1.0/D_0.1/com_pos.dat' u 1:5 w l

# unset multiplot

