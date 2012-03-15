set terminal %(terminal)s
set output"%(filename)s"

set size 1,1
unset key
set tmargin 1
set notitle
set multiplot
#set key bottom
set size 1,0.57
set ylabel"{/Symbol Db/b}_y" 1
set lmargin 7
set rmargin 2
set key 16000, -0.24 samplen 1
set origin 0,0
set xlabel "Longitudinal location [m]"
set xrange [%(xmin)f:%(xmax)f]

plot \
   "< awk '$6<20{print }' %(dir1)s/getbetay.out" u 2:(($4-$10)/$10):(sqrt($6**2+$5**2)/$10)  t "%(title1)s" w e pt 7 lw 2 ps 0.8 lt 1, \
   "< awk '$6<20{print }' %(dir2)s/getbetay.out" u 2:(($4-$10)/$10):(sqrt($6**2+$5**2)/$10)  t "%(title2)s" w e pt 7 lw 2 ps 0.8 lt 3

set nolabel
set origin 0, 0.5
set size 1,0.525
set noxtics
unset key
l=-0.3
set label "IR3" at 3000, l
set label "IR4" at 6000, l
set label "IR5" at 9100, l
set label "IR6" at 12500, l
set label "IR7" at 15800, l
set label "IR2" at 0, l
set label "IR8" at 19600, l
set label "IR1" at 22500, l
set ylabel"{/Symbol Db/b}_x"
set label "LHCB1 450 GeV" at graph 1,0.95 right
set xlabel""

plot  \
    "< awk '$6<20{print }' %(dir1)s/getbetax.out" u 2:(($4-$10)/$10):(sqrt($6**2+$5**2)/$10) w e pt 7 lw 2 ps 0.8 lt 1, \
    "< awk '$6<20{print }' %(dir2)s/getbetax.out" u 2:(($4-$10)/$10):(sqrt($6**2+$5**2)/$10) w e pt 7 lw 2 ps 0.8 lt 3
