set terminal postscript eps enhanced solid "Helvetica" 14
set output "egaf.eps"
set tmargin 0
set bmargin 0
set lmargin 1
set rmargin 1

set multiplot layout 2,1 margins 0.15,0.95,.1,.95 spacing 0,0.02
#set ylabel "{/Symbol t}


set ylabel "{/Symbol t}_{slow} (ns)"
set ytics ("1" 1e-9, "10" 1e-8, "100" 1e-7, "1000" 1e-6)
set yrange [1e-9:1e-6]
set xtics ("" 10, "" 20, "" 30, "" 40, "" 50)
set grid xtics

set xrange [0:56]
set logscale y
plot 'egaf/final.dat' u 1:($4 > 1e-1 ? $4 * 1e-9 : NaN):($3>0?$3:NaN) w points palette pt 7 notitle

set ylabel "{/Symbol t}_{fast} (ps)"
set xtics (10, 20, 30, 40, 50)
unset ytics
set ytics ("1" 1e-12, "100" 1e-10, "10,000" 1e-8)
set yrange [1e-12:1e-8]
set logscale y
set grid xtics

set xrange [0:56]
set xlabel "Peptide plane number (^{15}N_{i})"
plot 'egaf/final.dat' u 1:($5 > 1e-4 ? $5 * 1e-9 : NaN):($3>0?$3:NaN) w points palette pt 7 notitle

unset multiplot
set terminal postscript eps enhanced solid "Helvetica" 14
set output "egaf_slow.eps"
set tmargin 0
set bmargin 0
set lmargin 1
set rmargin 1

set multiplot layout 3,1 margins 0.15,0.95,.1,.95 spacing 0,0.02
#set ylabel "{/Symbol t}

unset xlabel
set ylabel "{/Symbol s}_{A}^{slow}"
set xtics ("" 10, "" 20, "" 30, "" 40, "" 50)
set grid xtics
set ytics (10, 20, 30)
unset logscale y
set xrange [0:56]
set yrange [0:30]
unset xlabel
plot 'egaf/final.dat' u 1:($6 > 0 ? $6 * (180./3.141593) : NaN):($3>0?$3:NaN) w points palette pt 7 notitle


set ylabel "{/Symbol s}_{B}^{slow}"
set xtics ("" 10, "" 20, "" 30, "" 40, "" 50)
set grid xtics
unset logscale y
set xrange [0:56]
set yrange [0:30]

plot 'egaf/final.dat' u 1:($7 > 0 ? $7 * (180./3.141593) : NaN):($3>0?$3:NaN) w points palette pt 7 notitle


set ylabel "{/Symbol s}_{G}^{slow}"
set xtics (10, 20, 30, 40, 50)
set grid xtics
unset logscale y
set xrange [0:56]
set yrange [0:30]
set xlabel "Peptide plane number (^{15}N_{i})"

plot 'egaf/final.dat' u 1:($8 > 0 ? $8 * (180./3.141593) : NaN):($3>0?$3:NaN) w points palette pt 7 notitle

unset multiplot
set terminal postscript eps enhanced solid "Helvetica" 14
set output "egaf_fast.eps"
set tmargin 0
set bmargin 0
set lmargin 1
set rmargin 1

set multiplot layout 1,1 margins 0.15,0.95,.1,.95 spacing 0,0.02
#set ylabel "{/Symbol t}
unset xlabel



set ylabel "S^{2}_{fast}"
set xtics (10, 20, 30, 40, 50)
set grid xtics
unset logscale y
set xrange [0:56]
set yrange [0.5:1]
set xlabel "Peptide plane number (^{15}N_{i})"
set ytics (0.6, 0.8, 1.0)

plot 'egaf/final.dat' u 1:($9 > 0 ? $9 : NaN):($3>0?$3:NaN) w points palette pt 7 notitle
