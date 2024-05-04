//
//  main.cpp
//  homework3
//
//  Created by Avik Ghosh on 4/30/24.
//

#include <iostream>
#include <chrono>
#include "DiffusionModel.hpp"

int main(int argc, const char * argv[]) {
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    
    int pecletVal = 100;
    int spatialPoints = 40;
    double timeLength = 0.05; // How many units of time are ran for
    double timeStep = 1 / 1000000.0;

    auto modelPureFE = DiffusionModel(pecletVal, spatialPoints, timeLength, timeStep, pureForwardEuler);

    std::filesystem::path fpathPureFE("FEData.csv");
    std::vector<int> pointsOfInterest {0};
    modelPureFE.outputData(fpathPureFE, pointsOfInterest);
    
    auto modelFELFMix = DiffusionModel(pecletVal, spatialPoints, timeLength, timeStep, forwardEulerIsh);

    std::filesystem::path fpathFELFMix("FELFrogMixData.csv");
    modelFELFMix.outputData(fpathFELFMix, pointsOfInterest);
    
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << "[µs]" << std::endl;

    return 0;
}
