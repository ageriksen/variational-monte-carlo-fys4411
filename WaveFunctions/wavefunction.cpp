#include "wavefunction.h"
#include <iostream>



WaveFunction::WaveFunction(System* system) {
    m_system = system;
}


void WaveFunction::updateParameters(double upd) {
    m_numberOfParameters = m_numberOfParameters + 1;
    // std::cout << "here 3 " << m_numberOfParameters << std::endl;
    m_parameters.push_back(upd);
}