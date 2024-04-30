#include <iostream>
#include "DiffusionModel.h"
#include "filesystem"

int main() {
	std::cout << "Hello world!" << std::endl;

	int pecletVal = 100;
	int spatialPoints = 40;
	double timeLength = 0.05; // How many units of time are ran for 
	double timeStep = 1 / 10000.0;
	auto method = pureForwardEuler;	

	auto model = DiffusionModel(pecletVal, spatialPoints, timeLength, timeStep, method);

	std::filesystem::path fpath("FEData.csv");
	std::vector<int> pointsOfInterest {0};
	model.outputData(fpath, pointsOfInterest);

	return 0;
}
