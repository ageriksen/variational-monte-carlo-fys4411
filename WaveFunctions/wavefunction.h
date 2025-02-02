#pragma once
#include <vector>


class WaveFunction {
public:
    WaveFunction(class System* system);
    int     getNumberOfParameters() { return m_numberOfParameters; }
    std::vector<double> getParameters() { return m_parameters; }
    void updateParameters(double upd);
    virtual double evaluate(std::vector<class Particle*> particles) = 0;
    virtual double computeDerivative(std::vector<class Particle*> particles) = 0;
    // virtual double computeDoubleDerivative(double r2) = 0;
    virtual double computeDoubleDerivative(std::vector<class Particle*> particles) = 0;

protected:
    int     m_numberOfParameters = 0;
    std::vector<double> m_parameters = std::vector<double>();
    class System* m_system = nullptr;
};

