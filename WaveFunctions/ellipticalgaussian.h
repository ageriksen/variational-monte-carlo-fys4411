#pragma once
#include "wavefunction.h"

class EllipticalGaussian : public WaveFunction {
public:
    EllipticalGaussian(class System* system, double alpha, double beta, double dt);
    double evaluate(std::vector<class Particle*> particles);
    double computeDerivative(std::vector<class Particle*> particles);
    double computeDoubleDerivative(std::vector<class Particle*> particles);
    // double computeDoubleDerivative(std::vector<class Particle*> particles, int particle);

private:
    double m_beta[3] = {1,1,1};
};
