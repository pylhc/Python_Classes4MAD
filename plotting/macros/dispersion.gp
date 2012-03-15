set terminal %(terminal)s
set output "disp.esp"

set multiplot

set size 1,0.5
set origin 0, 0.5
unset xtics

set xrange [%(xmin)f:%(xmax)f]
set yrange [%(ymin)f:%(ymax)f]

set ylabel "{/Symbol D}D_{x} [m]"

p "%(dir1)s/getDx.out" u 2:($4-$7):5 w e notitle lt 1 lw 2 ps 0.8 pt 7,\
  "%(dir2)s/getDx.out" u 2:($4-$7):5 w e notitle lt 3 lw 2 ps 0.8 pt 7

set origin 0,0
set xtics
set xlabel "Longitudinal position[m]"
set ylabel "{/Symbol D}D_{y} [m]"

set yrange [%(ymin2)f:%(ymax2)f]

p [:27000][-0.15:0.15] "%(dir1)s/getDy.out" u 2:($4-$7):5 w e notitle lt 1 lw 2 ps 0.8 pt 7,\
             "%(dir2)s/getDy.out" u 2:($4-$7):5 w e notitle lt 3 lw 2 ps 0.8 pt 7
