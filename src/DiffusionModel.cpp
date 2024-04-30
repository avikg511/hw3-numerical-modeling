#include "DiffusionModel.h"
#include <vector>
#include <numeric>
#include <fstream>

DiffusionModel::DiffusionModel(double pecletVal, int spatialPoints, 
                                double timeLength, const double timeStep, 
                                TimeSteppingApproach method) {

        this->dt = timeStep;
        this->spatialPoints = spatialPoints;
        this->method = method;

        // Spatial resolution (points) * time resolution (length / stepSize)
        this->data.resize(int(1.0 * spatialPoints * timeLength / timeStep));
        this->timeSeries.resize(int(timeLength / timeStep));

        // std::cout << "timeseries len: " << timeSeries.size() << std::endl;
        
        // std::iota here will start the vector at 0 and increment by 1 until the end.
        std::iota(this->timeSeries.begin(), this->timeSeries.end(), 0);
        for (auto it = timeSeries.begin(); it != timeSeries.end(); it++) {
            *it *= timeStep;
            if (std::distance(timeSeries.begin(), it) < 20) {
                std::cout << *it << ", ";
            }
        }

        // Set up the initial concentrations:
        // We have int spatialPoints divisions for our data, and we're forcing the x range to be from 0 to 1 (not inclusive), so 
        // 1 / spatialPoints gives us the change in x from one step to the next. The initial condition we're setting up is
        // that x < 1/2 ==> C(x,0) = 1, else C(x,0) = 0.

        // This is just using an np.linspace equivalent to put the index at each position. Probably inefficient but makes 
        // it simpler to write up, and it's not really a huge cost anyways.
        for (auto it = data.begin(); it != data.end(); it++) {
            if (std::distance(data.begin(), it) < spatialPoints / 2) {
                *it = 1.0;
            } else if (std::distance(data.begin(), it) < spatialPoints) {
                *it = 0.0;
            } else {
                break;
            }
        }

        std::cout << "Model setup complete!" << std::endl;
        if (this->method == pureForwardEuler) {
            std::cout << "Running forward euler!" << std::endl;
            RunPureForwardEuler();
        } else if (this->method == forwardEulerIsh) {
            std::cout << "Running kinda sorta forward euler ish with leapfrog!" << std::endl;
            RunForwardEulerIsh();
        }
}

inline int DiffusionModel::calcIndex(int xStep, int timeInd) {
    if (timeInd * this->spatialPoints + xStep > this->data.size()) {
        std::cout << "Messing up with xstep: " << xStep << " and timeStep: " << timeInd << std::endl;
        assert( 1 == 0 ); // Force stop for debugging help
    }
    return timeInd * this->spatialPoints + xStep;
}

void DiffusionModel::RunPureForwardEuler() {
    // The way we calculate Forward Euler is written below. Note that C_n here is the nth space step, and independent of time.
    // \vec{C}_n = \vec{C}_{n-1} + \Delta \vec{f}(\vec{C}), where we define \vec{f} as below:
    // \vec{f}(C_{k+1}, C_{k}, C_{k-1}) = \frac{ \partial \vec{C} }{ \partial t }  \\
    //                                  =  \frac{ -Pe * \Delta x + 1 }{ (\Delta x)^2 } * C_{k+1} + C_{k} * \frac{ -2 }{ (\Delta x)^2 } \\
    //                                        + \frac{ Pe * \Delta x + 1 }{ (\Delta x)^ 2} * C_{k-1} \\
    // (The above should theoretically compile in LaTeX :)

    // Let's define 3 constants, \alpha, \beta, \gamma. \alpha is the coefficient that multiples C_{k+1}, 
    // \beta multiplies C_{k} and \gamma multiplies C_{k-1}
    // Then \vec{f}(C_{k+1}, C_{k}, C_{k-1}) = \frac{ \partial \vec{C} }{ \partial t }  \\
    //                                       =  \alpha * C_{k+1} + \beta * C_{k} + \gamma * C_{k-1}

    double delX = 1.0 / spatialPoints;
    
    double alpha = ( -1 * this->pecletVal * delX + 1 ) / std::pow(delX, 2);
    double beta = -2 / std::pow(delX, 2);
    double gamma = ( this->pecletVal * delX + 1 ) / std::pow(delX, 2);

    std::cout << "delX: " << delX << " alpha: " << alpha << " beta: " << beta << " gamma: " << gamma << std::endl;  

    double nMinus1ValLastTime;
    double nValLastTime;
    double nPlus1ValLastTime;

    double del;
    int lastTimeIndex;
    std::vector<double>::iterator dataIt;
    for (auto it = this->timeSeries.begin() + 1, dataIt = this->data.begin();
                 it != this->timeSeries.end(); it++) {
            
        // We have alpha, beta, and gamma.
        // C_{n+1}^cur = C_{n} + Delta(t) * (alpha * C_{n+1}_last + beta * C_n_last + gamma * C_{n-1}_last)
        lastTimeIndex = std::distance(this->timeSeries.begin(), it);
        for (int xStep = 0; xStep < this->spatialPoints; xStep++, dataIt++) {
            if (xStep == 0) {
                nMinus1ValLastTime = this->data.at(calcIndex(spatialPoints - 1, lastTimeIndex - 1));
                nValLastTime = this->data.at(calcIndex(xStep, lastTimeIndex - 1));
                nPlus1ValLastTime = this->data.at(calcIndex(xStep + 1, lastTimeIndex - 1));
            } else if (xStep == spatialPoints - 1) {
                nMinus1ValLastTime = this->data.at(calcIndex(xStep - 1, lastTimeIndex - 1));
                nValLastTime = this->data.at(calcIndex(xStep, lastTimeIndex - 1));
                nPlus1ValLastTime = this->data.at(calcIndex(0, lastTimeIndex - 1));
            } else {
                nMinus1ValLastTime = this->data.at(calcIndex(xStep - 1, lastTimeIndex - 1));
                nValLastTime = this->data.at(calcIndex(xStep, lastTimeIndex - 1));
                nPlus1ValLastTime = this->data.at(calcIndex(xStep + 1, lastTimeIndex - 1));
            }
            
            del = dt * (alpha * nPlus1ValLastTime + beta * nValLastTime + gamma * nMinus1ValLastTime);
            *dataIt =  data.at(calcIndex(xStep, lastTimeIndex)) + del;
            
        }
    }
}

void DiffusionModel::RunForwardEulerIsh() {

}

void DiffusionModel::outputData(std::filesystem::path outputFile, std::vector<int> pointsOfInterest) {
    std::cout << "Outputting to file!" << std::endl;
    std::ofstream myfile;
    myfile.open(outputFile);

    if (myfile.is_open()) {
        // Here, the points of interest are basically the spatial point number, where time stretches over the entire 
        // array.
        for (auto ind : pointsOfInterest) {
            for (auto it = this->timeSeries.begin(); it != this->timeSeries.end(); it++) {
                myfile << *it << "," << this->data.at(calcIndex(ind, std::distance(this->timeSeries.begin(), it))) << std::endl;
            }

            myfile << std::endl;
        }
    }
    myfile.close();
}
