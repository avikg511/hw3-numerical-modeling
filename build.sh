cmake -B ./out -S . -DCMAKE_BUILD_TYPE=Debug;
cd out;
rm 'FEData.csv';

make;
./DiffusionAdvectionWork;

gnuplot -persist <<-EOFMarker
	reset session
	set datafile separator ","
	plot 'FEData.csv'  
	quit
EOFMarker

cd ..;
