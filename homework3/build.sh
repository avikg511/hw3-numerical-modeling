gnuplot -persist <<-EOFMarker
	reset session
	set xlabel 'time (dimensionless)'
	set ylabel 'concentration'
	set title 'Forward Euler, dt = 10^{-6}'
#    set logscale y
	set datafile separator ","
	plot 'FEData.csv' w l
	quit
EOFMarker

gnuplot -persist <<-EOFMarker
    reset session
    set xlabel 'time (dimensionless)'
    set ylabel 'concentration (log)'
    set title 'Forward Euler with LeapFrog'
    set logscale y
    set datafile separator ","
    plot 'FELFrogMix10-6.csv' w l, 'FELFrogMix10-4.csv' w l
    quit
EOFMarker
