set terminal postscript eps enhanced solid "Helvetica" 14
set output "smft.eps"
set tmargin 0
set bmargin 0
set lmargin 1
set rmargin 1

set multiplot layout 3,1 margins 0.15,0.95,.1,.95 spacing 0,0.02
#set ylabel "{/Symbol t}


set ylabel "S^{2}"
set ytics (0.4, 0.6, 0.8, 1.0)
set xtics ("" 10, "" 20, "" 30, "" 40, "" 50)
set grid xtics

set xrange [0:56]
set yrange [0.7:1]
plot 'smft/errors.dat' u 1:($6 > 0 ? $6 : NaN):($7>0?$7:NaN) w yerrorbars notitle

set ylabel "{/Symbol t}"
#set ytics (0.8, 0.9, 1.0)
set xtics ("" 10, "" 20, "" 30, "" 40, "" 50)
unset ytics
set ytics (0.1, 1, 10, 100, 1000)
set logscale y
set grid xtics
unset ytics
set ytics
set autoscale
set xrange [0:56]

plot 'smft/errors.dat' u 1:($4 > 0 ? $4 * 1e-9 : NaN):($5>0?$5 * 1e-9:NaN) w yerrorbars notitle


set ylabel "Ea / kJ/mol"
set ytics ("0" 0, "20" 20000, "40" 40000, "60" 60000)
set xtics ("" 10, "" 20, "" 30, "" 40, "" 50)
set grid xtics
set autoscale
unset logscale y
set xrange [0:56]

set xlabel "Peptide plane number (^{15}N_{i})"
set xtics (10, 20, 30, 40, 50)
plot 'smft/errors.dat' u 1:($8 > 0 ? $8 : NaN):($9>0?$9:NaN) w yerrorbars notitle

