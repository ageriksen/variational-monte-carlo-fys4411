#include "repulsiveinteraction.h"
#include <cassert>
#include <iostream>
#include <cmath>
#include "../system.h"
#include "../particle.h"
#include "../WaveFunctions/wavefunction.h"

using std::cout;
using std::endl;

RepulsiveInteraction::RepulsiveInteraction(System* system, double omega, double gamma, bool mode) :
        Hamiltonian(system) {
    assert(omega > 0);
	m_mode = mode;
    m_omega  = omega;
    m_gamma[2] = gamma*gamma; //
    // m_a = a;
    // cout << "here 1" << endl;
    // m_system->getWaveFunction()->updateParameters(gamma);
    // cout << "here 2" << endl;
    // system->getWaveFunction()->m_numberOfParameters += 1;
    // system->getWaveFunction()->m_parameters.push_back(gamma);
}

double RepulsiveInteraction::computeLocalEnergy(std::vector<Particle*> particles) {
    /* Here, you need to compute the kinetic and potential energies. Note that
     * when using numerical differentiation, the computation of the kinetic
     * energy becomes the same for all Hamiltonians, and thus the code for
     * doing this should be moved up to the super-class, Hamiltonian.
     *
     * You may access the wave function currently used through the
     * getWaveFunction method in the m_system object in the super-class, i.e.
     * m_system->getWaveFunction()...
     */

    double r2 = 0;
    for (int par = 0; par < m_system->getNumberOfParticles(); par++)
    {
        for( int dim=0; dim < m_system->getNumberOfDimensions(); dim++ )
        {            
            r2 += m_gamma[dim] * std::pow(particles[par]->getPosition()[dim], 2);
        }
    }
    
    double potentialEnergy	= 0.5*r2;
    // double potentialEnergy	= 0.5*m_omega*m_omega*r2;
    double kineticEnergy = 0;

    if( m_mode ){
        kineticEnergy	=
            -0.5*m_system->getWaveFunction()->computeDoubleDerivative(particles);
    } else
    {// mode here means numeric or analytic. TODO clarify variable
        kineticEnergy	= numeric();// -.5*numeric();
    }

    return kineticEnergy + potentialEnergy;
}

