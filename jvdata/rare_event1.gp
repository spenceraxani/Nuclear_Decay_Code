binwidth=0.00002
bin(x,width)=width*floor(x/width)

set title 'very rough histogram of binned averages'
set xlabel 'average cps'
set ylabel 'number'

threshold=0.99964235
ymax=15

set arrow from threshold,0 to threshold,ymax nohead lt 2

plot[0.9995:1.0005][0:ymax] 'rare_event1.txt' using (bin($1,binwidth)):(1.0) smooth freq with boxes

pause -1
