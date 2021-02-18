#include "system.h"
#include <cassert>
#include <cmath>
#include "sampler.h"
#include "particle.h"
#include "WaveFunctions/wavefunction.h"
#include "Hamiltonians/hamiltonian.h"
#include "InitialStates/initialstate.h"
#include "Math/random.h"


System::System() {
    m_random = new Random();
}

System::System(int seed) {
    m_random = new Random(seed);
}

bool System::metropolisStep(int particle) {
    /* Perform the actual Metropolis step: Choose a particle at random and
     * change it's position by a random amount, and check if the step is
     * accepted by the Metropolis test (compare the wave function evaluated
     * at this new position with the one at the old position).
	 * double step = m_stepLength*(Random::nextDouble() - .5);
     */

	std::vector<double> step(m_numberOfDimensions);
	double wfold = m_waveFunction->evaluate(m_particles);

	for( int dim = 0; dim < m_numberOfDimensions; dim++ )
	{
		step[dim] = m_stepLength*2*(m_random->nextDouble() - .5);
		m_particles[particle]->adjustPosition(step[dim], dim);
	}

	double wfnew = m_waveFunction->evaluate(m_particles);

	if( m_random->nextDouble() <= 2*(wfnew - wfold) )
	{
		return true;
	}
	else
	{
		for( int dim = 0; dim < m_numberOfDimensions; dim++ )
		{
			m_particles[particle]->adjustPosition(-step[dim], dim);
		}
		return false;
	}

    //return false;
}

void System::runMetropolisSteps(int numberOfMetropolisSteps) {
    m_particles                 = m_initialState->getParticles();
    m_sampler                   = new Sampler(this);
    m_numberOfMetropolisSteps   = numberOfMetropolisSteps;
    m_sampler->setNumberOfMetropolisSteps(numberOfMetropolisSteps);

    for (int i=0; i < numberOfMetropolisSteps; i++) {
		for( int particle=0; particle < m_numberOfParticles; particle++ )
		{
			bool acceptedStep = metropolisStep(particle);

			/* Here you should sample the energy (and maybe other things using
			* the m_sampler instance of the Sampler class. Make sure, though,
			* to only begin sampling after you have let the system equilibrate
			* for a while. You may handle this using the fraction of steps which
			* are equilibration steps; m_equilibrationFraction.
			*/
			if( i > numberOfMetropolisSteps*m_equilibrationFraction )
			{
				m_sampler->sample(acceptedStep);
			}
		}
    }
    m_sampler->computeAverages();
    m_sampler->getOutput();
    m_sampler->printOutputToTerminal();
    m_sampler->printOutputToFile();
}

void System::setNumberOfParticles(int numberOfParticles) {
    m_numberOfParticles = numberOfParticles;
}

void System::setNumberOfDimensions(int numberOfDimensions) {
    m_numberOfDimensions = numberOfDimensions;
}

void System::setStepLength(double stepLength) {
    assert(stepLength >= 0);
    m_stepLength = stepLength;
}

void System::setEquilibrationFraction(double equilibrationFraction) {
    assert(equilibrationFraction >= 0);
    m_equilibrationFraction = equilibrationFraction;
}

void System::setHamiltonian(Hamiltonian* hamiltonian) {
    m_hamiltonian = hamiltonian;
}

void System::setWaveFunction(WaveFunction* waveFunction) {
    m_waveFunction = waveFunction;
}

void System::setInitialState(InitialState* initialState) {
    m_initialState = initialState;
}


