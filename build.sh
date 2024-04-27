cmake -B ./out -S . -DCMAKE_BUILD_TYPE=Debug;
cd ./out
make;
./DiffusionAdvectionWork;
cd ..;
