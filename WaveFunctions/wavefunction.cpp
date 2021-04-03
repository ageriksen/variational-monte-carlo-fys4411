#include "wavefunction.h"


WaveFunction::WaveFunction(System* system) {
    m_system = system;
}


void WaveFunction::updateParameters(double update) {
    m_parameters.push_back(update);
    m_numberOfParameters += 1;
}