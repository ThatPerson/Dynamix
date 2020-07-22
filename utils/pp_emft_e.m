set terminal postscript eps enhanced solid "Helvetica" 14
set output "emft.eps"
set tmargin 0
set bmargin 0
set lmargin 1
set rmargin 1

set multiplot layout 6,1 margins 0.15,0.95,.1,.95 spacing 0,0.02
#set ylabel "{/Symbol t}

set ylabel "S^{2}_{slow}"
set ytics (0.8, 0.9, 1.0)
set xtics ("" 10, "" 20, "" 30, "" 40, "" 50)
set grid xtics
set xrange [0:56]
set yrange [0.80:1]

plot 'emft/errors.dat' u 1:($6 > 0 ? $6 : NaN):($7>0?$7:NaN) w yerrorbars notitle

set ylabel "{/Symbol t}_{slow}"
#set ytics (0.8, 0.9, 1.0)
unset xtics
unset ytics
set xtics ("" 10, "" 20, "" 30, "" 40, "" 50)
set grid xtics

#set ytics ("  10^1" 1e1, "  10^3" 1e3, "  10^5" 1e5)
unset ytics
set ytics
set logscale y
set autoscale

set xrange [0:56]
#set yrange [1:1e5]
set ytics ("-16" 1e-16, "-12" 1e-12, "-8" 1e-8)

plot 'emft/errors.dat' u 1:($4 > 0 ? $4 : NaN):($5>0?$5:NaN) w yerrorbars notitle


set ylabel "S^{2}_{fast}"
set ytics (0.4, 0.6, 0.8, 1.0)
set xtics ("" 10, "" 20, "" 30, "" 40, "" 50)
set grid xtics

set xrange [0:56]
set yrange [0.5:1]
plot 'emft/errors.dat' u 1:($6 > 0 ? $2/$6 : NaN):($7>0?$7:NaN) w yerrorbars notitle

set ylabel "{/Symbol t}_{fast}"
#set ytics (0.8, 0.9, 1.0)

#unset ytics
#set ytics (0.1, 1, 10, 100, 1000)
set logscale y
set grid xtics
unset ytics
set ytics ("-20" 1e-20, "-16" 1e-16, "-12" 1e-12)
set autoscale

set xrange [0:56]
#set yrange [0.1:1e3]
plot 'emft/errors.dat' u 1:($8 > 0 ? $8 : NaN):($9>0?$9:NaN) w yerrorbars notitle

unset logscale y
set ylabel "Ea_{slow} / kJ/mol"
set ytics ("0" 0, "20" 20000, "40" 40000, "60" 60000)
set xtics ("" 10, "" 20, "" 30, "" 40, "" 50)
set grid xtics

set xrange [0:56]
set yrange [0:60000]
plot 'emft/errors.dat' u 1:($10 > 0 ? $10 : NaN):($11>0?$11:NaN) w yerrorbars notitle

set ylabel "Ea_{fast} / kJ/mol"
set ytics ("0" 0, "20" 20000, "40" 40000, "60" 60000)
set xtics ("" 10, "" 20, "" 30, "" 40, "" 50)
set grid xtics

set xrange [0:56]
set yrange [0:60000]

set xlabel "Peptide plane number (^{15}N_{i})"
set xtics (10, 20, 30, 40, 50)
plot 'emft/errors.dat' u 1:($12 > 0 ? $12 : NaN):($13>0?$13:NaN) w yerrorbars notitle





