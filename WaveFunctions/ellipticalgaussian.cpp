#include "ellipticalgaussian.h"
#include <cmath>
#include <cassert>
#include <iostream>
#include "wavefunction.h"
#include "../system.h"
#include "../particle.h"

EllipticalGaussian::EllipticalGaussian(System* system, double alpha, double beta, double dt) :
        WaveFunction(system) {
    assert(alpha >= 0);
    assert(beta >= 0);
    assert(dt >= 0);
    m_numberOfParameters = 3;
    m_parameters.reserve(6);
    m_parameters.push_back(alpha);
    m_parameters.push_back(dt);
    m_parameters.push_back(beta);
    m_beta[2] = beta;
}

double EllipticalGaussian::evaluate(std::vector<class Particle*> particles) {
    /* You need to implement a Gaussian wave function here. The positions of
     * the particles are accessible through the particle[i].getPosition()
     * function.
     *
     * For the actual expression, use exp(-alpha * r^2), with alpha being the
     * (only) variational parameter.
     */
	double r2 = 0;
	for (int particle = 0; particle < m_system->getNumberOfParticles(); particle++)
    {
        for( int dim = 0; dim < m_system->getNumberOfDimensions(); dim++ )
        {
            r2 += m_beta[dim] * std::pow(particles[particle]->getPosition()[dim], 2);
        }
    }

    double r;
    double f = 1;
    if (m_parameters[4]>0)
    {
        for (int i = 0; i < m_system->getNumberOfParticles(); i++)
        {
            for (int j = i+1; j < m_system->getNumberOfParticles(); j++)
            {
                r = 0;
                for( int dim = 0; dim < m_system->getNumberOfDimensions(); dim++ )
                {
                    r += std::pow( particles[i]->getPosition()[dim] - particles[j]->getPosition()[dim] , 2);
                }
                if (r > m_parameters[4])
                {
                    // f = 0; //some large number
                    // // std::cout << "f=0" << std::endl;
                    // break;
                // }
                // else
                // {
                    f *= 1 - ( m_parameters[4] / std::sqrt(r) );
                }
                // f *= 1 - ( m_parameters[4] / std::sqrt(r) );
            }
        }
    }
    

	//return -m_parameters[0]*r2;
	return std::exp(-m_parameters[0]*r2) * f;
}

double EllipticalGaussian::computeDerivative(std::vector<class Particle*> particles) {
    /* Computes the derivative of wave function ansatz as function of variational parameters.
     * (I think this is correct)
     */
	double r2 = 0;
    for (int particle = 0; particle < m_system->getNumberOfParticles(); particle++)
    {
        for( int dim = 0; dim < m_system->getNumberOfDimensions(); dim++ )
        {
            r2 += m_beta[dim] * std::pow(particles[particle]->getPosition()[dim], 2);
        }
    }

    double r;
    double f = 1;
    if (m_parameters[4]>0)
    {
        for (int i = 0; i < m_system->getNumberOfParticles(); i++)
        {
            for (int j = i+1; j < m_system->getNumberOfParticles(); j++)
            {
                r = 0;
                for( int dim = 0; dim < m_system->getNumberOfDimensions(); dim++ )
                {
                    r += std::pow( particles[i]->getPosition()[dim] - particles[j]->getPosition()[dim] , 2);
                }
                if (r > m_parameters[4])
                {
                    // f = 0; //some large number
                    // // std::cout << "f=0" << std::endl;
                    // break;
                // }
                // else
                // {
                    f *= 1 - ( m_parameters[4] / std::sqrt(r) );
                }
            }
        }
    } 
    // std::cout << f << std::endl;
    // if (f!=1)
    // {
    //     std::cout << f << std::endl;
    // }
    

	return -.5*r2*f;
}

double EllipticalGaussian::computeDoubleDerivative(std::vector<class Particle*> particles) {
    /* All wave functions need to implement this function, so you need to
     * find the double derivative analytically. Note that by double derivative,
     * we actually mean the sum of the Laplacians with respect to the
     * coordinates of each particle.
     *
     * This quantity is needed to compute the (local) energy (consider the
     * Schrödinger equation to see how the two are related).
     */

	double r2 = 0;
    int N = m_system->getNumberOfParticles();
	for (int particle = 0; particle < N; particle++)
    {
        for( int dim = 0; dim < m_system->getNumberOfDimensions(); dim++ )
        {
            r2 += m_beta[dim] * std::pow(particles[particle]->getPosition()[dim], 2);
        }
    }
    double res0 =  4*pow(m_parameters[0], 2)*r2;
    res0 +=  2*m_parameters[0]*N*(m_parameters[2] + 2);

    double a = m_parameters[4];
    double res1 = 0;
    double r_kj;
    double denm = 0;
    for (int k = 0; k < N; k++)
    {
        //rk_rj = {0.,0.,0.};
        double rk_rj[3] = {0.,0.,0.};
        for (int j = 0; j < k; j++)
        {
            r_kj = 0;
            for( int dim = 0; dim < m_system->getNumberOfDimensions(); dim++ )
            {
                r_kj += abs(particles[k]->getPosition()[dim] - particles[j]->getPosition()[dim]);
                rk_rj[dim] += particles[k]->getPosition()[dim] - particles[j]->getPosition()[dim];
            }
            denm += pow(r_kj,3) - pow(r_kj,2) * a;
        }
        for (int j = k+1; j < N; j++)
        {
            r_kj = 0;
            for( int dim = 0; dim < m_system->getNumberOfDimensions(); dim++ )
            {
                r_kj += abs(particles[k]->getPosition()[dim] - particles[j]->getPosition()[dim]);
                rk_rj[dim] += particles[k]->getPosition()[dim] - particles[j]->getPosition()[dim];
            }
            denm += pow(r_kj,3) - pow(r_kj,2) * a;
        }
        for( int dim = 0; dim < m_system->getNumberOfDimensions(); dim++ )
        {
            res1 += m_beta[dim] * particles[k]->getPosition()[dim] * rk_rj[dim];
        }
    }
    res1 *= 4*m_parameters[0]*a/denm;

    double res2 = 0;
    for (int k = 0; k < N; k++)
    {
        for (int j = 0; j < k; j++)
        {
            r_kj = 0;
            for( int dim = 0; dim < m_system->getNumberOfDimensions(); dim++ )
            {
                r_kj += abs(particles[k]->getPosition()[dim] - particles[j]->getPosition()[dim]);
            }
            res2 += pow(r_kj,2) * pow(r_kj - a, 2);
        }
        for (int j = k+1; j < N; j++)
        {
            r_kj = 0;
            for( int dim = 0; dim < m_system->getNumberOfDimensions(); dim++ )
            {
                r_kj += abs(particles[k]->getPosition()[dim] - particles[j]->getPosition()[dim]);
            }
            res2 += pow(r_kj,2) * pow(r_kj - a, 2);
        }
    }
    res2 = pow(a,2) / res2;

	return res0+res1+res2;
}



    // since HarmonicOscillator::computeLocalEnergy allready finds r2, we can just reuse
    // that one. This might have to change in the future though, if we need the double
    // derivative somewhere else


































