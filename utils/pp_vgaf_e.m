set terminal postscript eps enhanced solid "Helvetica" 14
set output "vgaf.eps"
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
plot 'vgaf/errors.dat' u 1:($4 > 1e-1 ? $4 * 1e-9 : NaN):($5>0?$5 * 1e-9:NaN) w yerrorbars notitle

set ylabel "{/Symbol t}_{fast} (ps)"
set xtics (10, 20, 30, 40, 50)
unset ytics
set ytics ("1" 1e-12, "100" 1e-10, "10,000" 1e-8)
set yrange [1e-12:1e-8]
set logscale y
set grid xtics

set xrange [0:56]
set xlabel "Peptide plane number (^{15}N_{i})"
plot 'vgaf/errors.dat' u 1:($6 > 1e-4 ? $6 * 1e-9 : NaN):($7>0?$7 * 1e-9:NaN) w yerrorbars notitle

unset multiplot
set terminal postscript eps enhanced solid "Helvetica" 14
set output "vgaf_slow.eps"
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
plot 'vgaf/errors.dat' u 1:($8 > 0 ? $8 * (180./3.141593) : NaN):($9>0?$9* (180./3.141593):NaN) w yerrorbars notitle


set ylabel "{/Symbol s}_{B}^{slow}"
set xtics ("" 10, "" 20, "" 30, "" 40, "" 50)
set grid xtics
unset logscale y
set xrange [0:56]
set yrange [0:30]

plot 'vgaf/errors.dat' u 1:($10 > 0 ? $10 * (180./3.141593) : NaN):($11>0?$11* (180./3.141593):NaN) w yerrorbars notitle


set ylabel "{/Symbol s}_{G}^{slow}"
set xtics (10, 20, 30, 40, 50)
set grid xtics
unset logscale y
set xrange [0:56]
set yrange [0:30]
set xlabel "Peptide plane number (^{15}N_{i})"

plot 'vgaf/errors.dat' u 1:($12 > 0 ? $12 * (180./3.141593) : NaN):($13>0?$13* (180./3.141593):NaN) w yerrorbars notitle

unset multiplot
set terminal postscript eps enhanced solid "Helvetica" 14
set output "vgaf_fast.eps"
set tmargin 0
set bmargin 0
set lmargin 1
set rmargin 1

set multiplot layout 3,1 margins 0.15,0.95,.1,.95 spacing 0,0.02
#set ylabel "{/Symbol t}
unset xlabel

set ylabel "{/Symbol s}_{A}^{fast}"
set xtics ("" 10, "" 20, "" 30, "" 40, "" 50)
set grid xtics
unset logscale y
set xrange [0:56]
set yrange [0:30]
set ytics (10, 20, 30)
plot 'vgaf/errors.dat' u 1:($14 > 0 ? $14 * (180./3.141593) : NaN):($15>0?$15* (180./3.141593):NaN) w yerrorbars notitle


set ylabel "{/Symbol s}_{B}^{fast}"
set xtics ("" 10, "" 20, "" 30, "" 40, "" 50)
set grid xtics
unset logscale y
set xrange [0:56]
set yrange [0:30]
set ytics

plot 'vgaf/errors.dat' u 1:($16 > 0 ? $16 * (180./3.141593) : NaN):($17>0?$17* (180./3.141593):NaN) w yerrorbars notitle


set ylabel "{/Symbol s}_{G}^{fast}"
set xtics (10, 20, 30, 40, 50)
set grid xtics
unset logscale y
set xrange [0:56]
set yrange [0:30]
set xlabel "Peptide plane number (^{15}N_{i})"


plot 'vgaf/errors.dat' u 1:($18 > 0 ? $18 * (180./3.141593) : NaN):($19>0?$19* (180./3.141593):NaN) w yerrorbars notitle



unset multiplot
set terminal postscript eps enhanced solid "Helvetica" 14
set output "vgaf_orientation.eps"
set tmargin 0
set bmargin 0
set lmargin 1
set rmargin 1

set multiplot layout 3,1 margins 0.15,0.95,.1,.95 spacing 0,0.02
#set ylabel "{/Symbol t}

unset xlabel
set ylabel "{/Symbol a}"
set xtics ("" 10, "" 20, "" 30, "" 40, "" 50)
set grid xtics
set ytics
set ytics (-180, -120, -60, 0, 60, 120, 180)
unset logscale y
set xrange [0:56]
set yrange [-180:180]
unset xlabel
plot 'vgaf/errors.dat' u 1:($9 > 0 ? $20 * (180./3.141593) : NaN):($21>0?$21 * (180/3.141593):NaN) w yerrorbars  notitle

unset xlabel
set ylabel "{/Symbol b}"
set xtics ("" 10, "" 20, "" 30, "" 40, "" 50)
set grid xtics
set ytics
unset logscale y
set ytics (-180, -120, -60, 0, 60, 120, 180)
set xrange [0:56]
set yrange [-180:180]
unset xlabel
plot 'vgaf/errors.dat' u 1:($9 > 0 ? $22 * (180./3.141593) : NaN):($23>0?$23 * (180/3.141593):NaN) w yerrorbars  notitle

unset xlabel
set ylabel "{/Symbol g}"
set grid xtics
set ytics (-180, -120, -60, 0, 60, 120, 180)
unset logscale y
set xrange [0:56]
set yrange [-180:180]
set xtics (10, 20, 30, 40, 50)
set xlabel "Peptide plane number (^{15}N_{i})"
plot 'vgaf/errors.dat' u 1:($9 > 0 ? $24 * (180./3.141593) : NaN):($25>0?$25 * (180/3.141593):NaN) w yerrorbars notitle
