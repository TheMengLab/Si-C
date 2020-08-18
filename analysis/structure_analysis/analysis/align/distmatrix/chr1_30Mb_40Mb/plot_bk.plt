#set terminal png size 800,600
#set terminal png size 800,800
#set terminal png size 791.9,520.3
#set terminal png size 814.0,523.4
#set terminal png size 814.0,610.9
set terminal png size 814.0,666.0
set output "output.png"
#load '/home/group/code/gnuplot/colormap/gnuplot-colorbrewer-master/sequential/Reds.plt'
#load '/home/group/code/gnuplot/colormap/gnuplot-colorbrewer-master/diverging/Spectral.plt'
#set palette defined ( 0 0 0 0, 0.1667 0 0 1, 0.5 0 1 0, 0.8333 1 0 0, 1 1 1 1 )
#set palette defined ( 0 0.204 0.278 0.62, 0.1 0.271 0.624 0.835, 0.2 0.451 0.784 0.733, 0.3 0.624 0.8 0.592, 0.4 0.89 0.894 0.094 ,0.5 1 0.756 0.043 ,0.6 0.961 0.486 0.086 ,0.7 0.937 0.231 0.141 ,0.8 0.922 0.11 0.145 ,0.9 0.733 0.125 0.145 ,1 0.502 0.09 0.09)
#set palette defined (0 1 0 0 ,0.3 0.737 0.231 0.141,  0.5 1 1 1, 0.7 0.271 0.624 0.835, 1 0 0 1)
set palette defined (0 1 0 0 ,  0.5 1 1 1,  1 0 0 1)
set palette negative
set cbrange [-1.0:1.0]
set xrange[-0.5:500.5]
set yrange[-0.5:500.5]
set xtics 0,100,500
set xtics font "Times-Roman,0"
set ytics 0,100,500
set ytics font "Times-Roman,0"
set tics front
plot "distmatrix.dat" matrix with image

