set terminal postscript eps enhanced solid "Helvetica" 14
set output "demf.eps"
set tmargin 0
set bmargin 0
set lmargin 1
set rmargin 1

set multiplot layout 4,1 margins 0.15,0.95,.1,.95 spacing 0,0.02
#set ylabel "{/Symbol t}

set ylabel "S^{2}_{slow}"
set ytics (0.8, 0.9, 1.0)
set xtics ("" 10, "" 20, "" 30, "" 40, "" 50)
set grid xtics
set xrange [0:56]
set yrange [0.85:1]

plot 'demf/errors.dat' u 1:($6 > 0 ? $6 : NaN):($7>0?$7:NaN) w yerrorbars  notitle

set ylabel "{/Symbol t}_{slow} (ns)"
#set ytics (0.8, 0.9, 1.0)
unset xtics
unset ytics
set xtics ("" 10, "" 20, "" 30, "" 40, "" 50)
set grid xtics

set ytics ("  10^1" 1e1, "  10^3" 1e3, "  10^5" 1e5)
set logscale y
set xrange [0:56]
set yrange [1:1e5]

plot 'demf/errors.dat' u 1:($4 > 0 ? $4 : NaN):($5>0?$5:NaN) w yerrorbars  notitle


set ylabel "S^{2}_{fast}"
set ytics (0.4, 0.6, 0.8, 1.0)
set xtics ("" 10, "" 20, "" 30, "" 40, "" 50)
set grid xtics

set xrange [0:56]
set yrange [0.5:1]
plot 'demf/errors.dat' u 1:($10 > 0 ? $10 : NaN):($11>0?$11:NaN) w yerrorbars  notitle
# check errors on this one - S2f isn't directly calculated.

set ylabel "{/Symbol t}_{fast} (ps)"
#set ytics (0.8, 0.9, 1.0)
set xtics (10, 20, 30, 40, 50)
unset ytics
set ytics (0.1, 1, 10, 100, 1000)
set logscale y
set grid xtics

set xrange [0:56]
set yrange [0.1:1e3]
set xlabel "Peptide plane number (^{15}N_{i})"
plot 'demf/errors.dat' u 1:($8 > 1e-4 ? $8 * 1e3 : NaN):($9>0?$9*1e3:NaN) w yerrorbars  notitle

