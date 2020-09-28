set terminal postscript eps enhanced solid "Helvetica" 14
set output "vgaft.eps"
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
plot 'vgaft/errors.dat' u 1:($4 > 0 ? $4 * 1e-9 : NaN):($5>0?$5 * 1e-9:NaN) w yerrorbars notitle

set ylabel "{/Symbol t}_{fast} (s)"
set ytics ("-12" 1e-12, "-10" 1e-10, "-8" 1e-8, "-6" 1e-6)
set xtics ("" 10, "" 20, "" 30, "" 40, "" 50)
unset ytics
#set ytics ("-12" 1e-12, "-10" 1e-10, "-8" 1e-8, "-6" 1e-6)
set logscale y
set grid xtics
set ytics
set xrange [0:56]
plot 'vgaft/errors.dat' u 1:($6 > 0 ? $6 *1e-9: NaN):($7>0?$7*1e-9:NaN) w yerrorbars notitle

unset logscale y

set ylabel "Ea slow (kJ/mol)"
set ytics ("0" 0, "20" 20000, "40" 40000, "60" 60000)
set xtics ("" 10, "" 20, "" 30, "" 40, "" 50)
unset ytics
set ytics
set grid xtics

set xrange [0:56]
plot 'vgaft/errors.dat' u 1:($20 > 0 ? $20 : NaN):($21>0?$21:NaN) w yerrorbars notitle

set ylabel "Ea fast (kJ/mol)"
set ytics ("0" 0, "20" 20000, "40" 40000, "60" 60000)
set xtics (10, 20, 30, 40, 50)
unset ytics
set grid xtics
set ytics
set xrange [0:56]
set xlabel "Peptide plane number (^{15}N_{i})"
plot 'vgaft/errors.dat' u 1:($22 > 0 ? $22 : NaN):($23>0?$23:NaN) w yerrorbars notitle




unset multiplot

set terminal postscript eps enhanced solid "Helvetica" 14
set output "vgaft_slow.eps"
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

plot 'vgaft/errors.dat' u 1:($8 > 0 ? $8 * (180./3.141593) : NaN):($9>0?$9* (180./3.141593):NaN) w yerrorbars notitle


set ylabel "{/Symbol s}_{B}^{slow}"
set xtics ("" 10, "" 20, "" 30, "" 40, "" 50)
set grid xtics
unset logscale y
set xrange [0:56]
set yrange [0:30]

plot 'vgaft/errors.dat' u 1:($10 > 0 ? $10 * (180./3.141593) : NaN):($11>0?$11* (180./3.141593):NaN) w yerrorbars notitle


set ylabel "{/Symbol s}_{G}^{slow}"
set xtics (10,20,30,40,50)
set grid xtics
unset logscale y
set xrange [0:56]
set yrange [0:30]
set xlabel "Peptide plane number (^{15}N_{i})"

plot 'vgaft/errors.dat' u 1:($12 > 0 ? $12 * (180./3.141593) : NaN):($13>0?$13* (180./3.141593):NaN) w yerrorbars notitle

unset multiplot

set terminal postscript eps enhanced solid "Helvetica" 14
set output "vgaft_fast.eps"
set tmargin 0
set bmargin 0
set lmargin 1
set rmargin 1

set multiplot layout 3,1 margins 0.15,0.95,.1,.95 spacing 0,0.02
#set ylabel "{/Symbol t}


set ylabel "{/Symbol s}_{A}^{fast}"
set xtics ("" 10, "" 20, "" 30, "" 40, "" 50)
set grid xtics
unset logscale y
set xrange [0:56]
set yrange [0:30]

plot 'vgaft/errors.dat' u 1:($14 > 0 ? $14 * (180./3.141593) : NaN):($15>0?$15* (180./3.141593):NaN) w yerrorbars notitle


set ylabel "{/Symbol s}_{B}^{fast}"
set xtics ("" 10, "" 20, "" 30, "" 40, "" 50)
set grid xtics
unset logscale y
set xrange [0:56]
set yrange [0:30]

plot 'vgaft/errors.dat' u 1:($16 > 0 ? $16 * (180./3.141593) : NaN):($17>0?$17* (180./3.141593):NaN) w yerrorbars notitle


set ylabel "{/Symbol s}_{G}^{fast}"
set xtics (10,20,30,40,50)
set grid xtics
unset logscale y
set xrange [0:56]
set yrange [0:30]
set xlabel "Peptide plane number (^{15}N_{i})"

plot 'vgaft/errors.dat' u 1:($18 > 0 ? $18 * (180./3.141593) : NaN):($19>0?$19* (180./3.141593):NaN) w yerrorbars notitle
