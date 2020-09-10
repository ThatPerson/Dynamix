set terminal postscript eps enhanced solid "Helvetica" 14
set output "egaft.eps"
set tmargin 0
set bmargin 0
set lmargin 1
set rmargin 1

set multiplot layout 4,1 margins 0.15,0.95,.1,.95 spacing 0,0.02
#set ylabel "{/Symbol t}


set ylabel "{/Symbol t}_{slow} (s)"
#set ytics ("-12" 1e-12, "-10" 1e-10, "-8" 1e-8, "-6" 1e-6)
set xtics ("" 10, "" 20, "" 30, "" 40, "" 50)
set grid xtics
set ytics
set xrange [0:56]
set logscale y
plot 'egaft/final.dat' u 1:($4 > 0 ? $4 * 1e-9 : NaN):($3>0?$3 * 1e-9:NaN) w points palette pt 7 notitle

set ylabel "{/Symbol t}_{fast} (s)"
set ytics ("-12" 1e-12, "-10" 1e-10, "-8" 1e-8, "-6" 1e-6)
set xtics ("" 10, "" 20, "" 30, "" 40, "" 50)
unset ytics
#set ytics ("-12" 1e-12, "-10" 1e-10, "-8" 1e-8, "-6" 1e-6)
set logscale y
set grid xtics
set ytics
set xrange [0:56]
plot 'egaft/final.dat' u 1:($5 > 0 ? $5 *1e-9: NaN):($3>0?$3*1e-9:NaN) w points palette pt 7 notitle

unset logscale y

set ylabel "Ea slow (kJ/mol)"
set ytics ("0" 0, "20" 20000, "40" 40000, "60" 60000)
set xtics ("" 10, "" 20, "" 30, "" 40, "" 50)
unset ytics
set ytics
set grid xtics

set xrange [0:56]
plot 'egaft/final.dat' u 1:($10 > 0 ? $10 : NaN):($3>0?$3:NaN) w points palette pt 7 notitle

set ylabel "Ea fast (kJ/mol)"
set ytics ("0" 0, "20" 20000, "40" 40000, "60" 60000)
set xtics (10, 20, 30, 40, 50)
unset ytics
set grid xtics
set ytics
set xrange [0:56]
set xlabel "Peptide plane number (^{15}N_{i})"
plot 'egaft/final.dat' u 1:($11 > 0 ? $11 : NaN):($3>0?$3:NaN) w points palette pt 7 notitle




unset multiplot

set terminal postscript eps enhanced solid "Helvetica" 14
set output "egaft_slow.eps"
set tmargin 0
set bmargin 0
set lmargin 1
set rmargin 1

set multiplot layout 3,1 margins 0.15,0.95,.1,.95 spacing 0,0.02
#set ylabel "{/Symbol t}


set ylabel "{/Symbol s}_{A}^{slow}"
set xtics ("" 10, "" 20, "" 30, "" 40, "" 50)
set grid xtics
unset logscale y
set xrange [0:56]
set yrange [0:30]

plot 'egaft/final.dat' u 1:($6 > 0 ? $6 * (180./3.141593) : NaN):($3>0?$3:NaN) w points palette pt 7 notitle


set ylabel "{/Symbol s}_{B}^{slow}"
set xtics ("" 10, "" 20, "" 30, "" 40, "" 50)
set grid xtics
unset logscale y
set xrange [0:56]
set yrange [0:30]

plot 'egaft/final.dat' u 1:($7 > 0 ? $7 * (180./3.141593) : NaN):($3>0?$3:NaN) w points palette pt 7 notitle


set ylabel "{/Symbol s}_{G}^{slow}"
set xtics (10,20,30,40,50)
set grid xtics
unset logscale y
set xrange [0:56]
set yrange [0:30]
set xlabel "Peptide plane number (^{15}N_{i})"

plot 'egaft/final.dat' u 1:($8 > 0 ? $8 * (180./3.141593) : NaN):($3>0?$3:NaN) w points palette pt 7 notitle


unset multiplot
set terminal postscript eps enhanced solid "Helvetica" 14
set output "egaft_fast.eps"
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

plot 'egaft/final.dat' u 1:($9 > 0 ? $9 : NaN):($3>0?$3:NaN) w points palette pt 7 notitle
