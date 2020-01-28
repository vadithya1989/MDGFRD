  #include "externDefs.h"
#include "particle.h"
#include <algorithm> // remove and remove_if


void dissociationBD()
{
    ///all checks are made for the first scheduler list entry
    ///check if the particle is a LD particle
    if (particleList[schedulerList[0]].particleType == 2     &&
        ///and if the next reaction time is greater than tSim
        tSim >= particleList[schedulerList[0]].nextReactionTime &&
        ///next event type is reaction
        particleList[schedulerList[0]].nextEventType == 1
       )
    {
        std::cout << "*****************************************BD DISSOCIATION****************************************************" << std::endl;
        ///If the dissociating particle is a D particle
        if(particleList[schedulerList[0]].complexList[0].first != -1 &&
           particleList[schedulerList[0]].complexList[1].first != -1)
        {
            
            std::cout << " D complex BD dissociation" << std::endl;
            ///Draw a random number 0 or 1 to randomly dissociate a particle
            int n=1;
            int patchIndex = gsl_rng_uniform_int (rng, n);
            int otherPatchIndex=0;
            if(patchIndex == 0)
                otherPatchIndex = 1;
            
            ///get the ids of the particles that are dissociating
            int dissociatingEnzyme = (particleList[schedulerList[0]].complexList[patchIndex]).first;
            int dissociatingSubstrate = schedulerList[0];
            
            idleParticleList.erase( std::remove(idleParticleList.begin(), idleParticleList.end(), dissociatingEnzyme), idleParticleList.end() );
            
            
            ///////////////////////////////////ENZYME CHANGES////////////////////////////////////////////////
            if (particleList[schedulerList[0]].nextReactionType == 1)
            {
                ///For now the position of the dissociated particle is put just at a distance of 2Sigma
                myVector patchUnitVector = particleList[schedulerList[0]].q.bodyToSpace(speciesList[particleList[schedulerList[0]].speciesType].bodyPatchVectors[patchIndex])/ particleList[schedulerList[0]].q.bodyToSpace(speciesList[particleList[schedulerList[0]].speciesType].bodyPatchVectors[patchIndex]).magnitude();
                myVector newPosition = patchUnitVector*(2.0*prm.sigma);
                
                particleList[dissociatingEnzyme].particlePosition = newPosition;
                ///BD particle till domain is built
                particleList[dissociatingEnzyme].particleType = 2;
            }
            else if(particleList[schedulerList[0]].nextReactionType == 2)
            {
                ///put enzyme to other patch
                myVector patchUnitVector = particleList[schedulerList[0]].q.bodyToSpace(speciesList[particleList[schedulerList[0]].speciesType].bodyPatchVectors[otherPatchIndex])/ particleList[schedulerList[0]].q.bodyToSpace(speciesList[particleList[schedulerList[0]].speciesType].bodyPatchVectors[patchIndex]).magnitude();
                myVector newPosition = patchUnitVector*(2.0*prm.sigma);
                
                particleList[dissociatingEnzyme].particlePosition = newPosition;
                ///BD particle till domain is built
                particleList[dissociatingEnzyme].particleType = 2;
            }
            ///if it phosphorylates, then enzyme remains in the same place
            
            ///////////////////////////////////SUBSTRATE CHANGES////////////////////////////////////////////////
            
            ///Change the state of the patch
            //If the unbinding enzyme is a kinase then the substrate is phosphorylated
            else if(particleList[schedulerList[0]].nextReactionType == 3)
            {
                if(particleList[dissociatingEnzyme].speciesType == 3)
                {
                    particleList[dissociatingSubstrate].patchState[patchIndex] = 1;
                    if( particleList[dissociatingSubstrate].patchState[otherPatchIndex] == 1)
                        particleList[dissociatingSubstrate].speciesType = 2;
                    else if (particleList[dissociatingSubstrate].patchState[otherPatchIndex] == 0)
                        particleList[dissociatingSubstrate].speciesType = 1;
                    
                    //if the other patch is also phosphorylated then turn ON D switch
                    if(particleList[dissociatingSubstrate].patchState[otherPatchIndex] == 1)
                    {
                        ///counter for doubly phosphorylated (D particle) D increases
                        tBoundEndD = tSim;
                        tBoundD = tBoundD + numberOfDParticles * (tBoundEndD - tBoundStartD);
                        probabilityOfBeingBoundD = tBoundD/(maxDparticles*tSim);
                        std::cout << "D " << numberOfDParticles << " " << (tBoundEndD - tBoundStartD) << " " << tSim << " " << probabilityOfBeingBoundD << std::endl;
                        tBoundStartD = tSim;
                        numberOfDParticles = numberOfDParticles + 1;
                        ///counter for singly phosphorylated (C particle) C decreases
                        tBoundEndC = tSim;
                        tBoundC = tBoundC + numberOfCParticles * (tBoundEndC - tBoundStartC);
                        probabilityOfBeingBoundC = tBoundC/(maxCparticles*tSim);
                        std::cout << "C "  << numberOfCParticles << " " << (tBoundEndC - tBoundStartC) << " " << tSim << " " << probabilityOfBeingBoundC << std::endl;
                        tBoundStartC = tSim;
                        numberOfCParticles = numberOfCParticles - 1;
                        
                    }
                    //if the other patch is de-phosphorylated then turn ON C switch
                    else if(particleList[dissociatingSubstrate].patchState[otherPatchIndex] == 0)
                    {
                        ///counter for singly phosphorylated (C particle) C increases
                        tBoundEndC = tSim;
                        tBoundC = tBoundC + numberOfCParticles * (tBoundEndC - tBoundStartC);
                        probabilityOfBeingBoundC = tBoundC/(maxCparticles*tSim);
                        std::cout << "C "<< numberOfCParticles << " " << (tBoundEndC - tBoundStartC) << " " << tSim << " " << probabilityOfBeingBoundC << std::endl;
                        tBoundStartC = tSim;
                        numberOfCParticles = numberOfCParticles + 1;
                    }
                }
                else if(particleList[dissociatingEnzyme].speciesType == 4)
                {
                    particleList[dissociatingSubstrate].patchState[patchIndex] = 0;
                    if( particleList[dissociatingSubstrate].patchState[otherPatchIndex] == 1)
                        particleList[dissociatingSubstrate].speciesType = 1;
                    else if (particleList[dissociatingSubstrate].patchState[otherPatchIndex] == 0)
                        particleList[dissociatingSubstrate].speciesType = 0;
                    
                    
                   //if the other patch is phosphorylated then turn ON D switch
                    if(particleList[dissociatingSubstrate].patchState[otherPatchIndex] == 1)
                    {
                        ///counter for doubly phosphorylated (D particle) D increases
                        tBoundEndD = tSim;
                        tBoundD = tBoundD + numberOfDParticles * (tBoundEndD - tBoundStartD);
                        probabilityOfBeingBoundD = tBoundD/(maxDparticles*tSim);
                        std::cout << " D " << numberOfDParticles << " " << (tBoundEndD - tBoundStartD) << " " << tSim << " " << probabilityOfBeingBoundD << std::endl;
                        tBoundStartD = tSim;
                        numberOfDParticles = numberOfDParticles - 1;
                        ///counter for singly phosphorylated (C particle) C decreases
                        tBoundEndC = tSim;
                        tBoundC = tBoundC + numberOfCParticles * (tBoundEndC - tBoundStartC);
                        probabilityOfBeingBoundC = tBoundC/(maxCparticles*tSim);
                        std::cout << " C " << numberOfCParticles << " " << (tBoundEndC - tBoundStartC) << " " << tSim << " " << probabilityOfBeingBoundC << std::endl;
                        tBoundStartC = tSim;
                        numberOfCParticles = numberOfCParticles + 1;
                        
                    }
                    //if the other patch is de-phosphorylated then turn ON C switch
                    else if(particleList[dissociatingSubstrate].patchState[otherPatchIndex] == 0)
                    {
                        ///counter for singly phosphorylated (C particle) C increases
                        tBoundEndC = tSim;
                        tBoundC = tBoundC + numberOfCParticles * (tBoundEndC - tBoundStartC);
                        probabilityOfBeingBoundC = tBoundC/(maxCparticles*tSim);
                        std::cout<< "C " << numberOfCParticles << " " << (tBoundEndC - tBoundStartC) << " " << tSim << " " << probabilityOfBeingBoundC << std::endl;
                        tBoundStartC = tSim;
                        numberOfCParticles = numberOfCParticles - 1;
                    }
                }
            }
            
            ///draw next event times (dissociation times)
            double nextDissociationTime     = (1.0/dissociationRate) * (-log(gsl_rng_uniform(rng)));
            double nextPhosphorylationTime  = (1.0/phosphorylationRate) * (-log(gsl_rng_uniform(rng)));
            double nextHoppingTime          = (1.0/hoppingRate) * (-log(gsl_rng_uniform(rng)));
            
            ///the minimum of the above times is the next evebt time
            if( (nextDissociationTime < nextPhosphorylationTime) && (nextDissociationTime < nextHoppingTime) )
            {
                particleList[dissociatingSubstrate].nextEventTime = nextDissociationTime;
                particleList[dissociatingSubstrate].nextEventTime += tSim;
                particleList[dissociatingSubstrate].nextReactionTime = nextDissociationTime;
                particleList[dissociatingSubstrate].nextReactionTime += tSim;
                particleList[dissociatingSubstrate].nextReactionType = 1;
                particleList[dissociatingSubstrate].nextEventType = 1;
            }
            else if( (nextPhosphorylationTime < nextDissociationTime) && (nextPhosphorylationTime < nextHoppingTime) )
            {
                particleList[dissociatingSubstrate].nextEventTime = nextPhosphorylationTime;
                particleList[dissociatingSubstrate].nextEventTime += tSim;
                particleList[dissociatingSubstrate].nextReactionTime = nextPhosphorylationTime;
                particleList[dissociatingSubstrate].nextReactionTime += tSim;
                particleList[dissociatingSubstrate].nextReactionType = 3;
                particleList[dissociatingSubstrate].nextEventType = 1;
            }
            else if ( (nextHoppingTime < nextPhosphorylationTime) && (nextHoppingTime < nextDissociationTime) )
            {
                particleList[dissociatingSubstrate].nextEventTime = nextHoppingTime;
                particleList[dissociatingSubstrate].nextEventTime += tSim;
                particleList[dissociatingSubstrate].nextReactionTime = nextHoppingTime;
                particleList[dissociatingSubstrate].nextReactionTime += tSim;
                particleList[dissociatingSubstrate].nextReactionType = 2;
                particleList[dissociatingSubstrate].nextEventType = 1;
            }
            
            ///Remove the dissociated enzyme from the complex list of the substrate and viceversa
            //particleList[dissociatingEnzyme].complexList.pop_back();
            //particleList[dissociatingSubstrate].complexList.erase(particleList[dissociatingSubstrate].complexList.begin()+patchIndex);
            ///Remove the substrate from the scheduler list
            //schedulerList.erase(schedulerList.begin());
            //Add the new substrate to the scheduler list
            //schedulerList.push_back(dissociatingSubstrate);
            //change complex list of the substrate and enzyme
            std::pair<int, int>  tempPair(-1,-1);
            particleList[dissociatingSubstrate].complexList[particleList[dissociatingEnzyme].complexList[0].second] = tempPair;
            particleList[dissociatingEnzyme].complexList[0] = tempPair;
        }
        
        ///If the dissociating particle is a C particle
        else if((particleList[schedulerList[0]].complexList[0].first == -1 &&
                 particleList[schedulerList[0]].complexList[1].first != -1)
                ||
                (particleList[schedulerList[0]].complexList[1].first == -1 &&
                 particleList[schedulerList[0]].complexList[0].first != -1)
                )
        {
            std::cout << " C complex BD dissociation" << std::endl;

            int dissociatingEnzyme=0;
            //get the ids of the dissociating particles
            if(particleList[schedulerList[0]].complexList[0].first == -1 &&
               particleList[schedulerList[0]].complexList[1].first != -1)
            {
                std::cout << "dfh" << std::endl;
                dissociatingEnzyme = (particleList[schedulerList[0]].complexList[1]).first;
            }
            else if(particleList[schedulerList[0]].complexList[1].first == -1 &&
               particleList[schedulerList[0]].complexList[0].first != -1)
            {
                std::cout << "dfgjf`" << std::endl;
                dissociatingEnzyme = (particleList[schedulerList[0]].complexList[0]).first;
            }
            std::cout << "sdfg" << std::endl;
            
            int dissociatingSubstrate = schedulerList[0];
            std::cout << "ag" << std::endl;
            std::cout <<dissociatingEnzyme << std::endl;
            int patchIndex = particleList[dissociatingEnzyme].complexList[0].second;
            std::cout << "sghsh" << std::endl;
            int otherPatchIndex=0;
            if(patchIndex == 0)
                otherPatchIndex = 1;
            idleParticleList.erase( std::remove(idleParticleList.begin(), idleParticleList.end(), dissociatingEnzyme), idleParticleList.end() );
            std::cout << "dfgj" << std::endl;
            
            ///////////////////////////////////ENZYME CHANGES////////////////////////////////////////////////
            if (particleList[schedulerList[0]].nextReactionType == 1)
            {
                ///For now the position of the dissociated particle is put just at a distance of 2Sigma
                myVector patchUnitVector = particleList[schedulerList[0]].q.bodyToSpace(speciesList[particleList[schedulerList[0]].speciesType].bodyPatchVectors[patchIndex])/ particleList[schedulerList[0]].q.bodyToSpace(speciesList[particleList[schedulerList[0]].speciesType].bodyPatchVectors[patchIndex]).magnitude();
                myVector newPosition = patchUnitVector*(2.0*prm.sigma);
                
                particleList[dissociatingEnzyme].particlePosition = newPosition;
                ///BD particle till domain is built
                particleList[dissociatingEnzyme].particleType = 2;
            }
            else if(particleList[schedulerList[0]].nextReactionType == 2)
            {
                ///put enzyme to other patch
                myVector patchUnitVector = particleList[schedulerList[0]].q.bodyToSpace(speciesList[particleList[schedulerList[0]].speciesType].bodyPatchVectors[otherPatchIndex])/ particleList[schedulerList[0]].q.bodyToSpace(speciesList[particleList[schedulerList[0]].speciesType].bodyPatchVectors[patchIndex]).magnitude();
                myVector newPosition = patchUnitVector*(2.0*prm.sigma);
                
                particleList[dissociatingEnzyme].particlePosition = newPosition;
                ///BD particle till domain is built
                particleList[dissociatingEnzyme].particleType = 2;
            }
            ///Change the state of the patch
            ///If the unbinding enzyme is a kinase then the substrate is phosphorylated
            else if(particleList[schedulerList[0]].nextReactionType == 3)
            {
                if(particleList[dissociatingEnzyme].speciesType == 3)
                {
                    
                    particleList[dissociatingSubstrate].patchState[patchIndex] = 1;
                    if( particleList[dissociatingSubstrate].patchState[otherPatchIndex] == 1)
                        particleList[dissociatingSubstrate].speciesType = 2;
                    else if (particleList[dissociatingSubstrate].patchState[otherPatchIndex] == 0)
                        particleList[dissociatingSubstrate].speciesType = 1;
                    //if the other patch is also phosphorylated then turn ON D switch
                    if(particleList[dissociatingSubstrate].patchState[otherPatchIndex] == 1)
                    {
                        ///counter for doubly phosphorylated (D particle) D increases
                        tBoundEndD = tSim;
                        tBoundD = tBoundD + numberOfDParticles * (tBoundEndD - tBoundStartD);
                        probabilityOfBeingBoundD = tBoundD/(maxDparticles*tSim);
                        std::cout << "D " << numberOfDParticles << " " << (tBoundEndD - tBoundStartD) << " " << tSim << " " << probabilityOfBeingBoundD << std::endl;
                        tBoundStartD = tSim;
                        numberOfDParticles = numberOfDParticles + 1;
                        ///counter for singly phosphorylated (C particle) C decreases
                        tBoundEndC = tSim;
                        tBoundC = tBoundC + numberOfCParticles * (tBoundEndC - tBoundStartC);
                        probabilityOfBeingBoundC = tBoundC/(maxCparticles*tSim);
                        std::cout << "C " << numberOfCParticles << " " << (tBoundEndC - tBoundStartC) << " " << tSim << " " << probabilityOfBeingBoundC << std::endl;
                        tBoundStartC = tSim;
                        numberOfCParticles = numberOfCParticles - 1;
                        
                    }
                    //if the other patch is de-phosphorylated then turn ON C switch
                    else if(particleList[dissociatingSubstrate].patchState[otherPatchIndex] == 0)
                    {
                        ///counter for singly phosphorylated (C particle) C increases
                        tBoundEndC = tSim;
                        tBoundC = tBoundC + numberOfCParticles * (tBoundEndC - tBoundStartC);
                        probabilityOfBeingBoundC = tBoundC/(maxCparticles*tSim);
                        std::cout << "C "<< numberOfCParticles << " " << (tBoundEndC - tBoundStartC) << " " << tSim << " " << probabilityOfBeingBoundC << std::endl;
                        tBoundStartC = tSim;
                        numberOfCParticles = numberOfCParticles + 1;
                    }
                }
                else if(particleList[dissociatingEnzyme].speciesType == 4)
                {
                    std::cout << "now here " << std::endl;
                    particleList[dissociatingSubstrate].patchState[patchIndex] = 0;
                    if( particleList[dissociatingSubstrate].patchState[otherPatchIndex] == 1)
                        particleList[dissociatingSubstrate].speciesType = 1;
                    else if (particleList[dissociatingSubstrate].patchState[otherPatchIndex] == 0)
                        particleList[dissociatingSubstrate].speciesType = 0;
                    
                    //if the other patch is phosphorylated then turn ON D switch
                    if(particleList[dissociatingSubstrate].patchState[otherPatchIndex] == 1)
                    {
                        ///counter for doubly phosphorylated (D particle) D increases
                        tBoundEndD = tSim;
                        tBoundD = tBoundD + numberOfDParticles * (tBoundEndD - tBoundStartD);
                        probabilityOfBeingBoundD = tBoundD/(maxDparticles*tSim);
                        std::cout << "D " << numberOfDParticles << " " << (tBoundEndD - tBoundStartD) << " " << tSim << " " << probabilityOfBeingBoundD << std::endl;
                        tBoundStartD = tSim;
                        numberOfDParticles = numberOfDParticles - 1;
                        ///counter for singly phosphorylated (C particle) C decreases
                        tBoundEndC = tSim;
                        tBoundC = tBoundC + numberOfCParticles * (tBoundEndC - tBoundStartC);
                        probabilityOfBeingBoundC = tBoundC/(maxCparticles*tSim);
                        std::cout << "C " << numberOfCParticles << " " << (tBoundEndC - tBoundStartC) << " " << tSim << " " << probabilityOfBeingBoundC << std::endl;
                        tBoundStartC = tSim;
                        numberOfCParticles = numberOfCParticles + 1;
                        
                    }
                    //if the other patch is de-phosphorylated then turn ON C switch
                    else if(particleList[dissociatingSubstrate].patchState[otherPatchIndex] == 0)
                    {
                        std::cout << "still here " << std::endl;
                        ///counter for singly phosphorylated (C particle) C increases
                        tBoundEndC = tSim;
                        tBoundC = tBoundC + numberOfCParticles * (tBoundEndC - tBoundStartC);
                        probabilityOfBeingBoundC = tBoundC/(maxCparticles*tSim);
                        std::cout << "C " << numberOfCParticles << " " << (tBoundEndC - tBoundStartC) << " " << tSim << " " << probabilityOfBeingBoundC << std::endl;
                        tBoundStartC = tSim;
                        numberOfCParticles = numberOfCParticles - 1;
                    }
                }
            }
            
            ///draw next event times (dissociation times)
            double nextDissociationTime     = (1.0/singleParticleDissociationRate) * (-log(gsl_rng_uniform(rng)));
            
            particleList[dissociatingSubstrate].nextEventTime = nextDissociationTime;
            particleList[dissociatingSubstrate].nextEventTime += tSim;
            particleList[dissociatingSubstrate].nextReactionTime = nextDissociationTime;
            particleList[dissociatingSubstrate].nextReactionTime += tSim;
            particleList[dissociatingSubstrate].nextReactionType = 0;
            particleList[dissociatingSubstrate].nextEventType = 1;
            

            ///Remove the dissociated enzyme from the complex list of the substrate and viceversa
            //particleList[dissociatingSubstrate].complexList.pop_back();
            //particleList[dissociatingEnzyme].complexList.pop_back();
            
            ///Remove the substrate from the scheduler list
            //schedulerList.erase(schedulerList.begin());
            //Add the new substrate to the scheduler list
            //schedulerList.push_back(dissociatingSubstrate);
            std::pair<int, int>  tempPair(-1,-1);
            particleList[dissociatingSubstrate].complexList[particleList[dissociatingEnzyme].complexList[0].second] = tempPair;
            particleList[dissociatingEnzyme].complexList[0] = tempPair;
            
        }
      std::cout << "*********************************************************************************************" << std::endl;
    }
    
}







































