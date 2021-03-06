##
# Start of plot..
##
set terminal %(terminal)s
set output "%(filename)s"

set size 1,1
unset key
set tmargin 1
set notitle
set multiplot
#set key bottom
set size 1,0.57
set ylabel "%(ylabel1)s" 1
set lmargin 7
set rmargin 2
set key 16000, -0.24 samplen 1
set origin 0,0
set xlabel "Longitudinal location [m]"
set xrange [%(xmin)s:%(xmax)s]
set yrange [%(ymin)s:%(ymax)s]

plot \
   "< awk '$6<20{print }' %(dir1)s/get%(data1_p1)s.out" u %(xfunc)s:%(yfunc)s:%(errfunc)s  t "%(title1)s" w e pt 7 lw 2 ps 0.8 lt 1

set nolabel
set origin 0, 0.5
set size 1,0.525

l=%(yminl)s
%(xtics_plot2)s

set ylabel "%(ylabel2)s"
set label "%(title)s" at graph 1,0.95 right
set xlabel ""
set yrange [%(ymin2)s:%(ymax2)s]

plot  \
    "< awk '$6<20{print }' %(dir1)s/get%(data1_p2)s.out" u %(xfunc)s:%(yfunc2)s:%(errfunc2)s w e pt 7 lw 2 ps 0.8 lt 1

unset multiplot
reset
##
# End of plot
##


