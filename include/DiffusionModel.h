#include <vector>
#include <iostream>
#include <filesystem>
enum TimeSteppingApproach {
    pureForwardEuler = 0,
    forwardEulerIsh = 1,
};


class DiffusionModel {
    public:
        std::vector<double> data;
        std::vector<double> timeSeries;

        double pecletVal;
        int spatialPoints;
        double dt;
        TimeSteppingApproach method;

        DiffusionModel(double pecletVal, 
                        int spatialPoints, 
                        double timeLength, 
                        const double timeStep, 
                        TimeSteppingApproach method
                    );
        void RunPureForwardEuler();
        void RunForwardEulerIsh();
        void outputData(std::filesystem::path outputFile, std::vector<int> pointsOfInterest);
        inline int calcIndex(int xStep, int timeStep);
};