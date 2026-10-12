set term pdfcairo size 16,16 enhanced background rgb "white"
set output 'stone_89_exp4.pdf'
set grid
dir_stone="stone/"
dir_falcon_12="FALCON/output_12/"
dir_falcon_13="FALCON/output_13/"
set key reverse Left
set key spacing 1.5
set key at 0.22, 90.465
set multiplot layout 2,2 title "Stone 89 - Experiment 4" font 'Helvetic,26'
set xtics 0.1 offset 0,0 font 'Helvetic,22'
set ytics 0.1 offset 0,0 font 'Helvetic,22'
set xlabel offset 0,-0.5 '{/Helvetic=26 Time}'
set ylabel offset -0.5,0 '{/Helvetic=26 Angle (°)}'
set xrange [0:0.75]
set yrange [89.5:90.5]
plot sprintf("%sold.data",dir_falcon_12) using ($1):2 w l lw 5 lc rgb '#4393C3' title '{/Helvetic=22 FALCON - theta old}',\
     sprintf("%snew.data",dir_falcon_13) using ($1):2 w l lw 5 lc rgb '#E07B2E' title '{/Helvetic=22 FALCON - theta new}',\
     sprintf("%sold.ascii",dir_stone) using 1:2 w l lw 5 dt 2 lc rgb '#8B3A6B' title '{/Helvetic=22 stone 89 - theta old}',\
     sprintf("%snew.ascii",dir_stone) using 1:2 w l lw 5 dt 2 lc rgb '#3D3D3D' title '{/Helvetic=22 stone 89 - theta new}'

set key at 0.2, 1.4
set xtics 0.1 offset 0,0 font 'Helvetic,22'
set ytics 0.5 offset 0,0 font 'Helvetic,22'
set xlabel offset 0,-0.5 '{/Helvetic=26 Time}'
set ylabel offset -0.5,0 '{/Helvetic=26 Strain}'
set xrange [0:0.75]
set yrange [-1.5:1.5]
plot sprintf("%sold.data",dir_falcon_12) using ($1):3 w l lw 5 lc rgb '#4393C3' title '{/Helvetic=22 FALCON - e1 old}',\
     sprintf("%snew.data",dir_falcon_13) using ($1):3 w l lw 5 lc rgb '#E07B2E' title '{/Helvetic=22 FALCON - e1 new}',\
     sprintf("%sold.ascii",dir_stone) using 1:3 w l lw 5 dt 2 lc rgb '#8B3A6B' title '{/Helvetic=22 stone 89 - e1 old}',\
     sprintf("%snew.ascii",dir_stone) using 1:3 w l lw 5 dt 2 lc rgb '#3D3D3D' title '{/Helvetic=22 stone 89 - e1 new}',\
     sprintf("%sold.data",dir_falcon_12) using ($1):4 w l lw 5 lc rgb '#4393C3' title '{/Helvetic=22 FALCON - e2 old}',\
     sprintf("%snew.data",dir_falcon_13) using ($1):4 w l lw 5 lc rgb '#E07B2E' title '{/Helvetic=22 FALCON - e2 new}',\
     sprintf("%sold.ascii",dir_stone) using 1:4 w l lw 5 dt 2 lc rgb '#8B3A6B' title '{/Helvetic=22 stone 89 - e2 old}',\
     sprintf("%snew.ascii",dir_stone) using 1:4 w l lw 5 dt 2 lc rgb '#3D3D3D' title '{/Helvetic=22 stone 89 - e2 new}'

set key at 0.18, 0.95
set xtics 0.1 offset 0,0 font 'Helvetic,22'
set ytics 0.2 offset 0,0 font 'Helvetic,22'
set xlabel offset 0,-0.5 '{/Helvetic=26 Time}'
set ylabel offset -0.5,0 '{/Helvetic=26 R tensor}'
set xrange [0:0.75]
set yrange [-0.2:1]
plot sprintf("%snew.data",dir_falcon_13) using ($1):5 w l lw 5 lc rgb '#4393C3' title '{/Helvetic=22 FALCON - Rxx}',\
     sprintf("%snew.data",dir_falcon_13) using ($1):6 w l lw 5 lc rgb '#E07B2E' title '{/Helvetic=22 FALCON - Rxy}',\
     sprintf("%snew.data",dir_falcon_13) using 1:7 w l lw 5 lc rgb '#8B3A6B' title '{/Helvetic=22 FALCON - Ryx}',\
     sprintf("%snew.data",dir_falcon_13) using 1:8 w l lw 5 lc rgb '#3D3D3D' title '{/Helvetic=22 FALCON - Ryy}',\
     sprintf("%snew.ascii",dir_stone) using ($1):5 w l lw 5 dt 2 lc rgb '#4393C3' title '{/Helvetic=22 stone 89 - Rxx}',\
     sprintf("%snew.ascii",dir_stone) using ($1):6 w l lw 5 dt 2 lc rgb '#E07B2E' title '{/Helvetic=22 stone 89 - Rxy}',\
     sprintf("%snew.ascii",dir_stone) using 1:7 w l lw 5 dt 2 lc rgb '#8B3A6B' title '{/Helvetic=22 stone 89 - Ryx}',\
     sprintf("%snew.ascii",dir_stone) using 1:8 w l lw 5 dt 2 lc rgb '#3D3D3D' title '{/Helvetic=22 stone 89 - Ryy}'

set key at 0.18, 2.4
set xtics 0.1 offset 0,0 font 'Helvetic,22'
set ytics 0.5 offset 0,0 font 'Helvetic,22'
set xlabel offset 0,-0.5 '{/Helvetic=26 Time}'
set ylabel offset -0.5,0 '{/Helvetic=26 V tensor}'
set xrange [0:0.75]
set yrange [-0.1:2.5]
plot sprintf("%snew.data",dir_falcon_13) using ($1):9 w l lw 5 lc rgb '#4393C3' title '{/Helvetic=22 FALCON - Vxx}',\
     sprintf("%snew.data",dir_falcon_13) using ($1):10 w l lw 5 lc rgb '#E07B2E' title '{/Helvetic=22 FALCON - Vxy}',\
     sprintf("%snew.data",dir_falcon_13) using 1:11 w l lw 5 lc rgb '#8B3A6B' title '{/Helvetic=22 FALCON - Vyx}',\
     sprintf("%snew.data",dir_falcon_13) using 1:12 w l lw 5 lc rgb '#3D3D3D' title '{/Helvetic=22 FALCON - Vyy}',\
     sprintf("%snew.ascii",dir_stone) using ($1):9 w l lw 5 dt 2 lc rgb '#4393C3' title '{/Helvetic=22 stone 89 - Vxx}',\
     sprintf("%snew.ascii",dir_stone) using ($1):10 w l lw 5 dt 2 lc rgb '#E07B2E' title '{/Helvetic=22 stone 89 - Vxy}',\
     sprintf("%snew.ascii",dir_stone) using 1:11 w l lw 5 dt 2 lc rgb '#8B3A6B' title '{/Helvetic=22 stone 89 - Vyx}',\
     sprintf("%snew.ascii",dir_stone) using 1:12 w l lw 5 dt 2 lc rgb '#3D3D3D' title '{/Helvetic=22 stone 89 - Vyy}'

unset multiplot