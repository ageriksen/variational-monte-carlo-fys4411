#include "repulsiveinteraction.h"
#include <cassert>
#include <iostream>
#include <cmath>
#include "../system.h"
#include "../particle.h"
#include "../WaveFunctions/wavefunction.h"

using std::cout;
using std::endl;

RepulsiveInteraction::RepulsiveInteraction(System* system, double omega, double gamma, double a, bool mode) :
        Hamiltonian(system) {
    assert(omega > 0);
	m_mode = mode;
    m_omega  = omega;
    m_gamma[2] = gamma*gamma; //
    m_a = a;
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

    double V_int = 0;
    double r_ij2 = 0;

    //first we check if V_int is 0 or inf
    for (int i = 0; i < m_system->getNumberOfParticles(); i++)
    {
        for (int j = i+1; j < m_system->getNumberOfParticles(); j++)
        {
            r_ij2 = 0;
            for( int dim=0; dim < m_system->getNumberOfDimensions(); dim++ )
            {            
                r_ij2 += std::pow(particles[i]->getPosition()[dim] - particles[j]->getPosition()[dim], 2);
            }
            if (r_ij2 <= m_a)
            {
                V_int = 1e6; //some large number
                break;
            }
        }
    }
    
    if (V_int>0)
    {
        return V_int;
    }
    else
    {
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
                -0.5*m_system->getWaveFunction()->computeDoubleDerivative(r2);
        } else
        {// mode here means numeric or analytic. TODO clarify variable
            kineticEnergy	= numeric();// -.5*numeric();
        }

        return kineticEnergy + potentialEnergy;
    }
}

