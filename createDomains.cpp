
#include "particle.h"
#include "main.h"

void burstAllPossibleDomainsLocally()
{
    blankCells();
    sortParticles();
    std::vector<particle>::iterator it = particleList.begin();
    int particleId = 0;
    while( it != particleList.end())
    {       
        int cellID = findCellOfTheParticle(*it);
        std::pair<double, int> nearestParticle = returnNearestDistance (cellID, particleId, "domain");
        double nearestDomainDistance = nearestParticle.first;
        int nnid  = nearestParticle.second;
        if( (nearestDomainDistance < prm.rc) && it->particleType == 2 && particleList[nnid].particleType == 1 && ( nnid != -1 && nearestDomainDistance != -1.0 ))
        {
            particleList[nnid].burstTheDomain("burst");
            it = particleList.begin();
            particleId = 0;
        }
        else
        {
            ++it;
            ++particleId;
        }
        
    }    
}

void createDomains()
{   
    burstAllPossibleDomainsLocally();
    checkIfDomainsCanBeBuiltOnLangevinParticles();
    sortSchedulerList();
}

