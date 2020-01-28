//#include <gsl/gsl_rng.h>
//#include <cmath>
#include <algorithm>
//#include <iostream>
//
//#include "simulationParameters.h"
#include "particle.h"
#include "main.h"
#include "externDefs.h"


void checkIfDomainsCanBeBuiltOnLangevinParticles()
{
    blankCells();
    sortParticles();

    std::vector<particle>::iterator it = particleList.begin();
    int particleId=0;
    while( it != particleList.end())
    {
        int cellID = findCellOfTheParticle(*it);
        std::pair<double, int> nearestParticle = returnNearestDistance (cellID, particleId, "neighbor");
        double nearestNeighborDistance = nearestParticle.first;
        int nnid  = nearestParticle.second;
        ///if there is only one particle draw domain of constant size at each step
        ///three particles form the D particle, so when they are bound, then
        ///it means there is only one particle in the system effective
        int pass=1;
        
        if((idleParticleList.size() == 2 && particleList.size()==3) ||
           (idleParticleList.size() == 1 && particleList.size()==2))
        {
            //std::cout << "hi" << std::endl;
            //std::cout << "single domain" << std::endl;
            nearestNeighborDistance = 20.0*prm.sigma;
            nnid  = 0;
        }
        else
        {
            //check if the current particle is in the idle particle list
            
            for(std::vector<int>::iterator ipl=idleParticleList.begin();ipl!=idleParticleList.end();++ipl)
            {
                if(particleId == *ipl)
                {
                    pass = 0;
                }
            }
        }
        //proceed only if the particle is not in the idleParticleList
        if(nearestNeighborDistance > 0 && it->particleType == 2 && pass == 1)
        {
            //std::cout << " domain built" << std::endl;
            
            it->buildADomain(nearestNeighborDistance);
            
            
        }
        
        ++it;
        ++particleId;
    }
}


bool wayToSort(int i, int j)
{
    return particleList[i].nextEventTime < particleList[j].nextEventTime;
}
///sorts the scheduler list based on the next event times of the particles from the particle list
/// the instruction to sort is given the the wayToSort function
void sortSchedulerList ()
{
    std::sort(schedulerList.begin(), schedulerList.end(), wayToSort);
}

