#set terminal png size 800,600
#set terminal png size 800,800
#set terminal png size 791.9,520.3
set terminal png size 814.0,523.4
set output "temp.png"
#load '/home/group/code/gnuplot/colormap/gnuplot-colorbrewer-master/sequential/Reds.plt'
#load '/home/group/code/gnuplot/colormap/gnuplot-colorbrewer-master/diverging/Spectral.plt'
set palette defined ( 0 0 0 0, 0.1667 0 0 1, 0.5 0 1 0, 0.8333 1 0 0, 1 1 1 1 )
#set palette defined (1 1 1 1,  0.8333 1 0 0, 0.5 0 1 0,  0.1667 0 0 1, 0 0 0 0)
#set palette negative
set cbrange [0:25]
set xrange[-0.5:20]
set yrange[-0.5:20]
set xtics -0.5,50,300
set xtics font "Times-Roman,0"
set ytics -0.5,50,300
set ytics font "Times-Roman,0"
set tics front 
set xlabel "snapshot ID"
set ylabel "snapshot ID"
#set xtics 50
#set ytics 50
plot "rmsdmatrix_mirror.dat" matrix with image

