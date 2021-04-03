#pragma once
#include "hamiltonian.h"
#include <vector>

class RepulsiveInteraction : public Hamiltonian {
public:
    RepulsiveInteraction(class System* system, double omega, double gamma, bool mode);
    double computeLocalEnergy(std::vector<class Particle*> particles);

private:
    double m_omega = 0;
    double m_gamma[3] = {1,1,1};
	bool m_mode = true;
};

