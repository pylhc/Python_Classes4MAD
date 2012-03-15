set terminal %(terminal)s
set output "%(filename)s"

set multiplot

set size 1,0.5
set origin 0, 0.5
unset xtics
l=-0.056
set label "IR3" at 3000, l
set label "IR4" at 6000, l
set label "IR5" at 9100, l
set label "IR6" at 12500, l
set label "IR7" at 15800, l
set label "IR2" at 0, l
set label "IR8" at 19600, l
set label "IR1" at 22500, l
set label "LHCB1 %(energy_in_gev)i GeV" at graph 1,1.05 right

set ylabel "{/Symbol D}D_{x}/{/Symbol b}_x^{1/2} [m^{1/2}]"
set key samplen 1
p "%(dir1)s/getNDx.out"  u 2:($4-$8):5 t "%(title1)s"  w e  lt 1 lw 2 ps 0.8 pt 7,\
  "%(dir2)s/getNDx.out"  u 2:($4-$8):5 t "%(title2)s"  w e lt 3 lw 2 ps 0.8 pt 7

set origin 0,0
set size 1,.54
set xtics 5000
set xlabel "Longitudinal position[m]"
set ylabel "{/Symbol D}D_{y} [m]"
unset label
p [:27000][-0.15:0.15] "%(dir1)s/getDy.out" u 2:($4-$7):5 w e notitle lt 1 lw 2 ps 0.8 pt 7,\
             "%(dir2)s/getDy.out" u 2:($4-$7):5 w e notitle lt 3 lw 2 ps 0.8 pt 7
