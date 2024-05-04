//
//  DiffusionModel.cpp
//  homework3
//
//  Created by Avik Ghosh on 4/30/24.
//

#include "DiffusionModel.hpp"
#include <vector>
#include <numeric>
#include <fstream>

DiffusionModel::DiffusionModel(double pecletVal, int spatialPoints,
                                double timeLength, const double timeStep,
                                TimeSteppingApproach method) {

        this->dt = timeStep;
        this->spatialPoints = spatialPoints;
        this->method = method;
        this->pecletVal = pecletVal;

        // Spatial resolution (points) * time resolution (length / stepSize)
        this->data.resize(int(1.0 * spatialPoints * timeLength / timeStep));
        this->timeSeries.resize(int(timeLength / timeStep));

        // std::cout << "timeseries len: " << timeSeries.size() << std::endl;
        
        // std::iota here will start the vector at 0 and increment by 1 until the end.
        std::iota(this->timeSeries.begin(), this->timeSeries.end(), 0);
        for (auto it = timeSeries.begin(); it != timeSeries.end(); it++) {
            *it *= timeStep;
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
            RunForwardEulerLeapFrogIsh();
        }
}

inline int DiffusionModel::calcIndex(int xStep, int timeInd) {
    return timeInd * this->spatialPoints + xStep;
}

void DiffusionModel::RunPureForwardEuler() {
    /*
     The way we calculate Forward Euler is written below. Note that C_k here is the nth
     step, and C^n as the nth timestep.

     \vec{C}_k = \vec{C}_k^{n-1} + \Delta t * \vec{f}(\vec{C}), where we define \vec{f} as below:
     \vec{f_k}(C_{k+1}^{n-1}, C_{k}^{n-1}, C_{k-1}^{n-1}) \\
                     = \frac{ \partial \vec{C_k} }{ \partial t }  \\
                     =   C_{k+1}^{n-1} * \frac{ -Pe * \Delta x + 1 }{ (\Delta x)^2 } \\
                      +    C_{k}^{n-1} * \frac{ -2 }{ (\Delta x)^2 } \\
                      +  C_{k-1}^{n-1} * \frac{ Pe * \Delta x + 1 }{ (\Delta x)^ 2}
   (The above snippets should theoretically compile in LaTeX)

    Let's define 3 constants, \alpha, \beta, \gamma.
         \alpha multiples  C_{k+1},
         \beta  multiplies C_{k}
         \gamma multiplies C_{k-1}
    Then \vec{f}(C_{k+1}, C_{k}, C_{k-1}) = \frac{ \partial \vec{C} }{ \partial t }  \\
                                          =  \alpha * C_{k+1} + \beta * C_{k} + \gamma * C_{k-1}
     */


    double delX = 1.0 / spatialPoints;
    
    double alpha = ( -1 * this->pecletVal * delX + 1 ) / std::pow(delX, 2);
    double beta = -2 / std::pow(delX, 2);
    double gamma = ( this->pecletVal * delX + 1 ) / std::pow(delX, 2);

    std::cout << "delX: " << delX << " alpha: " << alpha << " beta: " << beta << " gamma: " << gamma << std::endl;

    double cKMinus1TMinus1;
    double cKTMinus1;
    double cKPlus1TMinus1;

    double del;
    int tMinus1;
    
    typedef std::vector<double>::iterator doubleVecIter;
    for (std::pair<doubleVecIter, doubleVecIter> it(this->timeSeries.begin() + 1, this->data.begin() + spatialPoints) ; it.first != this->timeSeries.end(); it.first++) {
            
        // We have alpha, beta, and gamma.
        // C_{k+1}^{n+1} = C_{k}^{n} + \Delta(t) * (\alpha * C_{k+1}^{n} + \beta * C_{k}^{n} + \gamma * C_{k-1}^{n})
        tMinus1 = int(std::distance(this->timeSeries.begin(), it.first - 1));
        for (int xStep = 0; xStep < this->spatialPoints && it.second != this->data.end(); xStep++, it.second++) {
            if (xStep == 0) {
                cKMinus1TMinus1 = this->data.at(calcIndex(spatialPoints - 1, tMinus1));
                cKTMinus1 = this->data.at(calcIndex(xStep, tMinus1));
                cKPlus1TMinus1 = this->data.at(calcIndex(xStep + 1, tMinus1));
            } else if (xStep == spatialPoints - 1) {
                cKMinus1TMinus1 = this->data.at(calcIndex(xStep - 1, tMinus1));
                cKTMinus1 = this->data.at(calcIndex(xStep, tMinus1));
                cKPlus1TMinus1 = this->data.at(calcIndex(0, tMinus1));
            } else {
                cKMinus1TMinus1 = this->data.at(calcIndex(xStep - 1, tMinus1));
                cKTMinus1 = this->data.at(calcIndex(xStep, tMinus1));
                cKPlus1TMinus1 = this->data.at(calcIndex(xStep + 1, tMinus1));
            }
            
            del = dt * (alpha * cKPlus1TMinus1 + beta * cKTMinus1 + gamma * cKMinus1TMinus1);
            *it.second =  data.at(calcIndex(xStep, tMinus1)) + del;
            
        }
    }
}

void DiffusionModel::RunForwardEulerLeapFrogIsh() {
    /*
     This code is largely copied from RunPureForwardEuler() above because the coefficients are \\
         all that is changing. We now define a new function f' that represents our time derivative
     
      Again, note that C_k here is the kth space step, and C^n is the nth timestep
      \vec{C}^{n+1} = \vec{C}^{n-1} + 2\Delta t( -Pe \cdot \frac{ \partial \vec{C}^{n} }{ \partial x } \\
                   + \frac{ \partial^2 C^{n-1} }{ \partial x^2 } ) \\
    Or similarly:
        C_{k}^{n+1} = C_{k}^{n-1} + 2 \Delta t \cdot ( -Pe \cdot \frac{ C_{k+1}^n - C_{k-1}^n }{2 \Delta x} \\ + \frac{ C_{k+1}^{n-1} - 2 C_{k}^{n-1} + C_{k-1}^{n-1} }{ (\Delta x)^2 } ) \\
     
     With this form, there is no reason to create \alpha, \beta, \gamma or similar forms because we'd have
        to make 5 coefficients (using values from 2 time values). Instead, we'll just calculate the
        derivatives. Derivative 1 is going to be \frac{ C_{k+1}^n - C_{k-1}^n }{2 \Delta x}
        and derivative 2 is \frac{ C_{k+1}^{n-1} - 2 C_{k}^{n-1} + C_{k-1}^{n-1} }{ (\Delta x)^2 }
        
     */
    double delX = 1.0 / spatialPoints;

    double cKMinus1TMinus1;
    double cKTMinus1;
    double cKPlus1TMinus1;
    
    double cKPlus1T; // C(x = x_{k-1}, t = t_{n})
    double cKMinus1T;

//    Derivative One and Two defined above in the comment towards the end
    double derivOne;
    double derivTwo;
    double newContribution;
    
    int tMinus1;
    int t;
    
    
    typedef std::vector<double>::iterator doubleVecIter;
    for (std::pair<doubleVecIter, doubleVecIter> it(this->timeSeries.begin() + 1, this->data.begin() + spatialPoints) ; it.first != this->timeSeries.end(); it.first++) {
            
        tMinus1 = int(std::distance(this->timeSeries.begin(), it.first - 1));
        t = int(std::distance(this->timeSeries.begin(), it.first));

        for (int xStep = 0; xStep < this->spatialPoints && it.second != this->data.end(); xStep++, it.second++) {
            if (xStep == 0) {
                cKTMinus1 = this->data.at(calcIndex(xStep, tMinus1));
                cKPlus1TMinus1 = this->data.at(calcIndex(xStep + 1, tMinus1));
                cKMinus1TMinus1 = this->data.at(calcIndex(spatialPoints - 1, tMinus1));
                
                cKPlus1T = this->data.at(calcIndex(xStep + 1, t));
                cKMinus1T = this->data.at(calcIndex(spatialPoints - 1, t));
            } else if (xStep == spatialPoints - 1) {
                cKTMinus1 = this->data.at(calcIndex(xStep, tMinus1));
                cKPlus1TMinus1 = this->data.at(calcIndex(0, tMinus1));
                cKMinus1TMinus1 = this->data.at(calcIndex(xStep - 1, tMinus1));
                
                cKPlus1T = this->data.at(calcIndex(0, t));
                cKMinus1T = this->data.at(calcIndex(xStep - 1, t));
            } else {
                cKTMinus1 = this->data.at(calcIndex(xStep, tMinus1));
                cKPlus1TMinus1 = this->data.at(calcIndex(xStep + 1, tMinus1));
                cKMinus1TMinus1 = this->data.at(calcIndex(xStep - 1, tMinus1));
                
                cKPlus1T = this->data.at(calcIndex(xStep + 1, t));
                cKMinus1T = this->data.at(calcIndex(xStep - 1, t));
            }
            
            derivOne = (cKPlus1T - cKMinus1T) / (2 * delX);
            derivTwo = (cKPlus1TMinus1 - 2 * cKTMinus1 + cKMinus1TMinus1) / (std::pow(delX, 2));
            
            newContribution = 2 * dt * ( -1 * this->pecletVal * derivOne + derivTwo);
            *it.second =  data.at(calcIndex(xStep, tMinus1)) + newContribution;
            
        }
    }
}

void DiffusionModel::outputData(std::filesystem::path outputFile, std::vector<int> pointsOfInterest) {
    std::cout << "Outputting to file!" << std::endl;
    std::ofstream myfile;
    myfile.open(outputFile);

    if (myfile.is_open()) {
//         Here, the points of interest are basically the spatial point number, where \\
            time stretches over the entire array.
        for (auto ind : pointsOfInterest) {
            for (auto it = this->timeSeries.begin(); it != this->timeSeries.end(); it++) {
                myfile << *it << "," << this->data.at(calcIndex(ind, int(std::distance(this->timeSeries.begin(), it)))) << std::endl;
            }

            myfile << std::endl;
        }
    }
    myfile.close();
}
