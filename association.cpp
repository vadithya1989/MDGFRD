#include "externDefs.h"
#include "simulationParameters.h"
#include "particle.h"
#include "main.h"

void updateAssociationList( std::vector< std::tr1::tuple<int, int, int> > &associationList
                          )
{
    //Loop over all the particle pairs and calculate their respective energies
    // if the energy obeys the criteria then add them to the association list
    //if the list already has a particle that is one of the current particles
    //then the particle-pair with the least energy is added to the list and the other is discarded
    for(int i=0; i<particleList.size();++i)
    {
        
        for(int j=i+1; j<particleList.size();++j)
        {
            
            ///calculate interparticle and interpatch vectors
            myVector r;
            std::vector< std::pair<int, myVector> > rPatchVector;
            measureInterParticleDistace(r, particleList[i].particlePosition, particleList[j].particlePosition);
            measureInterPatchDistance  (r, rPatchVector, particleList[i], particleList[j]);
            ///determine which patch of the substrate is interacting
            int interactingPatch=0;

            if(rPatchVector[0].second.magnitude() < rPatchVector[1].second.magnitude())
                interactingPatch = rPatchVector[0].first;
            else
                interactingPatch = rPatchVector[1].first;
            
            /** check if the correct species are interacting and only kinases can interact with dephosphorylated patch and phosphotases can interact with phosphorylated patch and only if they are proceed*/
            if( ( (particleList[i].speciesType==0 && particleList[j].speciesType==3) ||
                  ( (particleList[i].speciesType==1 && particleList[j].speciesType==3) && (particleList[i].patchState[interactingPatch] == 0) ) ||
                  ( (particleList[i].speciesType==1 && particleList[j].speciesType==4) && (particleList[i].patchState[interactingPatch] == 1) ) ||
                  (particleList[i].speciesType==2 && particleList[j].speciesType==4)
                ) ||
                ( (particleList[i].speciesType==3 && particleList[j].speciesType==0) ||
                  ( (particleList[i].speciesType==3 && particleList[j].speciesType==1) && (particleList[j].patchState[interactingPatch] == 0) )||
                  ( (particleList[i].speciesType==4 && particleList[j].speciesType==1) && (particleList[j].patchState[interactingPatch] == 1) )||
                  (particleList[i].speciesType==4 && particleList[j].speciesType==2)
                )
              )
            {
                
                int pass = 1;///only if pass is 1 the function is executed for the particles
                int substrate=0;
                int enzyme=0;
                
                //if the particle i is a substrate and both the sites are occupied then pass=0
                
                if ( (particleList[i].speciesType==0) || (particleList[i].speciesType==1) || (particleList[i].speciesType==2) )
                {
                    
                    substrate = particleList[i].particleID;
                    enzyme = particleList[j].particleID;
                    
                    if (particleList[i].complexList[0].first != -1 && particleList[i].complexList[1].first != -1) 
                    {
                       
                        pass=0;
                    }
                    
                }
                
                //if the particle j is a substrate and both the sites are occupied then pass=0
                if( (particleList[j].speciesType==0) || (particleList[j].speciesType==1) || (particleList[j].speciesType==2) )
                {
                    
                    substrate = particleList[j].particleID;
                    enzyme = particleList[i].particleID;
                    //(particleList[j].complexList.size() == 2)
                    if( (particleList[j].complexList[0].first != -1 && particleList[j].complexList[1].first != -1))
                    {
                        pass=0;
                    }
                }
                
                for(std::vector<int>::iterator it = idleParticleList.begin(); it!=idleParticleList.end(); ++it)
                {
                    //if one of the particles is in the idle particle list, then it does not react with other particles
                    if(*it == i || *it ==j)
                    {
                        pass=0;
                    }
                }
                
                if(pass == 1)
                {
                    
                    double E=0.0;//energy bet pairs of particles
                    calcEnergy(r, rPatchVector, E, particleList[i], particleList[j]);
                    //std::cout << "energy of interaction between particles : "<<enzyme <<" and " << substrate << " is " << E/(prm.boltzmannConstant*prm.temperature) << std::endl;
                    ///tuple of particleID, particleID, patchID of the patch in the substrate
                    std::tr1::tuple<int, int, int> tempTuple = std::tr1::make_tuple(i, j, interactingPatch);
                    
                    //std::cout << "sherlock " << particleList[enzyme].speciesType << " " << particleList[substrate].patchState[interactingPatch] << " " << interactingPatch << std::endl;
                    //std::cout << "soole " <<particleList[substrate].patchState[1] << std::endl;
                    if( (E < -10.0*prm.boltzmannConstant*prm.temperature)  &&
                        ((particleList[enzyme].speciesType == 3 && particleList[substrate].patchState[interactingPatch] == 0) || (particleList[enzyme].speciesType == 4 && particleList[substrate].patchState[interactingPatch] == 1))
                       )
                    {
                        
                        std::cout << "Adding to association list " << interactingPatch << std::endl;
                        associationList.push_back(tempTuple);
                    }
                    
                }
                
                
            }
            
        }
    }      
}

void associationReaction(
                         std::vector<std::tr1::tuple<int, int,int> > &associationList
                        )
{
    
    ///Do this only if you have associating particles
    if(associationList.size() > 0)
    {
        std::cout << "******************************** PARTICLES ASSOCIATING*******************************   "<<std::endl;
        std::cout << "species : "<<particleList[std::tr1::get<0>(associationList[0])].speciesType << " and " << particleList[std::tr1::get<1>(associationList[0])].speciesType << " are reacting to the patch "<< std::tr1::get<2>(associationList[0])<< " of the substrate" << std::endl;
        for(int i=0; i<associationList.size();++i)
        {
            
//            std::cout << particleList[1].nextReactionTime << " " << particleList[1].nextReactionType << " " << particleList[1].nextEventTime << " " << particleList[1].nextEventType << std::endl;
            if( (particleList[std::tr1::get<0>(associationList[i])].speciesType == 0) ||
                (particleList[std::tr1::get<0>(associationList[i])].speciesType == 1) ||
                (particleList[std::tr1::get<0>(associationList[i])].speciesType == 2)
              )
            {
                idleParticleList.push_back(std::tr1::get<1>(associationList[i]));
                
                std::pair<int, int> temp1(std::tr1::get<1>(associationList[i]),0);
                particleList[std::tr1::get<0>(associationList[i])].complexList[std::tr1::get<2>(associationList[i])] = temp1;
                //particleList[std::tr1::get<0>(associationList[i])].complexList.push_back(temp1);
                
                std::pair<int, int> temp2(std::tr1::get<0>(associationList[i]),(std::tr1::get<2>(associationList[i])));
                //since enzyme has complexList.size()=1
                particleList[std::tr1::get<1>(associationList[i])].complexList[0] = temp2;
                //particleList[std::tr1::get<1>(associationList[i])].complexList.push_back(temp2);
                
                //counter
                //if the new complex formed is a D complex
                if(particleList[std::tr1::get<0>(associationList[i])].complexList[0].first != -1 &&
                   particleList[std::tr1::get<0>(associationList[i])].complexList[1].first != -1)
                {
                    std::cout << " double complex formed... " << std::endl;
                    std::cout <<particleList[std::tr1::get<1>(associationList[i])].complexList[0].first << " " << particleList[std::tr1::get<1>(associationList[i])].complexList[1].first << std::endl;
                    ///draw next event times (dissociation times or phosphorylation)
                    ///Can redraw time for the other bound particle and the current particle
                    ///Or can just draw once with twice the dissociation rate
                    double nextDissociationTime = (1.0/(2.0*blockedDissociationRate)) * (-log(gsl_rng_uniform(rng)));
                    double nextPhosphorylationTime =  (1.0/phosphorylationRate) * (-log(gsl_rng_uniform(rng)));
                    
                    /// get the next event and reaction times and types
                    if(nextDissociationTime < nextPhosphorylationTime )
                    {
                        particleList[std::tr1::get<0>(associationList[i])].nextReactionTime = nextDissociationTime;
                        particleList[std::tr1::get<0>(associationList[i])].nextReactionTime += tSim;
                        particleList[std::tr1::get<0>(associationList[i])].nextEventTime = nextDissociationTime;
                        particleList[std::tr1::get<0>(associationList[i])].nextEventTime += tSim;
                        particleList[std::tr1::get<0>(associationList[i])].nextReactionType = 1;
                        particleList[std::tr1::get<0>(associationList[i])].nextEventType = 1;
                    }
                    else if( nextPhosphorylationTime < nextDissociationTime )
                    {
                        particleList[std::tr1::get<0>(associationList[i])].nextReactionTime = nextPhosphorylationTime;
                        particleList[std::tr1::get<0>(associationList[i])].nextReactionTime += tSim;
                        particleList[std::tr1::get<0>(associationList[i])].nextEventTime = nextPhosphorylationTime;
                        particleList[std::tr1::get<0>(associationList[i])].nextEventTime += tSim;
                        particleList[std::tr1::get<0>(associationList[i])].nextReactionType = 3;
                        particleList[std::tr1::get<0>(associationList[i])].nextEventType = 1;
                    }
                    ///add the substrate particle to the scheduler list
                    //schedulerList.push_back(std::tr1::get<0>(associationList[i]));
                }
                //if the new complex formed is a C complex
                else if((particleList[std::tr1::get<0>(associationList[i])].complexList[0].first == -1 &&
                        particleList[std::tr1::get<0>(associationList[i])].complexList[1].first != -1)
                        ||
                        (particleList[std::tr1::get<0>(associationList[i])].complexList[1].first == -1 &&
                         particleList[std::tr1::get<0>(associationList[i])].complexList[0].first != -1)
                        )
                {
                    
                    std::cout << " single complex formed... " << std::endl;
                    //draw next event times (dissociation times)
                    double nextDissociationTime     = (1.0/dissociationRate) * (-log(gsl_rng_uniform(rng)));
                    double nextPhosphorylationTime  = (1.0/phosphorylationRate) * (-log(gsl_rng_uniform(rng)));
                    double nextHoppingTime          = (1.0/hoppingRate) * (-log(gsl_rng_uniform(rng)));
                    
                    ///the minimum of the above times is the next evebt time
                    if( (nextDissociationTime < nextPhosphorylationTime) && (nextDissociationTime < nextHoppingTime) )
                    {
                        particleList[std::tr1::get<0>(associationList[i])].nextEventTime = nextDissociationTime;
                        particleList[std::tr1::get<0>(associationList[i])].nextEventTime += tSim;
                        particleList[std::tr1::get<0>(associationList[i])].nextReactionTime = nextDissociationTime;
                        particleList[std::tr1::get<0>(associationList[i])].nextReactionTime += tSim;
                        particleList[std::tr1::get<0>(associationList[i])].nextReactionType = 1;
                        particleList[std::tr1::get<0>(associationList[i])].nextEventType = 1;
                    }
                    else if( (nextPhosphorylationTime < nextDissociationTime) && (nextPhosphorylationTime < nextHoppingTime) )
                    {
                        particleList[std::tr1::get<0>(associationList[i])].nextEventTime = nextPhosphorylationTime;
                        particleList[std::tr1::get<0>(associationList[i])].nextEventTime += tSim;
                        particleList[std::tr1::get<0>(associationList[i])].nextReactionTime = nextPhosphorylationTime;
                        particleList[std::tr1::get<0>(associationList[i])].nextReactionTime += tSim;
                        particleList[std::tr1::get<0>(associationList[i])].nextReactionType = 3;
                        particleList[std::tr1::get<0>(associationList[i])].nextEventType = 1;
                    }
                    else if ( (nextHoppingTime < nextPhosphorylationTime) && (nextHoppingTime < nextDissociationTime) )
                    {
                        particleList[std::tr1::get<0>(associationList[i])].nextEventTime = nextHoppingTime;
                        particleList[std::tr1::get<0>(associationList[i])].nextEventTime += tSim;
                        particleList[std::tr1::get<0>(associationList[i])].nextReactionTime = nextHoppingTime;
                        particleList[std::tr1::get<0>(associationList[i])].nextReactionTime += tSim;
                        particleList[std::tr1::get<0>(associationList[i])].nextReactionType = 2;
                        particleList[std::tr1::get<0>(associationList[i])].nextEventType = 1;
                    }
                    ///add the substrate particle to the scheduler list
                    //schedulerList.push_back(std::tr1::get<0>(associationList[i]));
                }
                //change the nextEventTime of the enzyme
                particleList[std::tr1::get<1>(associationList[i])].nextEventTime = particleList[std::tr1::get<1>(associationList[i])].nextReactionTime;
                particleList[std::tr1::get<1>(associationList[i])].nextEventType = 1;
                
                
                
                std::cout << " The next event time and type and the next reaction time and type : "<< particleList[std::tr1::get<0>(associationList[i])].nextEventTime << " " << particleList[std::tr1::get<0>(associationList[i])].nextEventType << " " << particleList[std::tr1::get<0>(associationList[i])].nextReactionTime << " " << particleList[std::tr1::get<0>(associationList[i])].nextReactionType << std::endl;
                
            }
            
            else if( (particleList[std::tr1::get<1>(associationList[i])].speciesType == 0) ||
                     (particleList[std::tr1::get<1>(associationList[i])].speciesType == 1) ||
                     (particleList[std::tr1::get<1>(associationList[i])].speciesType == 2)
                   )
            {
                idleParticleList.push_back(std::tr1::get<0>(associationList[i]));
                std::pair<int, int> temp1(std::tr1::get<0>(associationList[i]),0);
                 particleList[std::tr1::get<1>(associationList[i])].complexList[std::tr1::get<2>(associationList[i])] = temp1;
                //particleList[std::tr1::get<1>(associationList[i])].complexList.push_back(temp1);
                std::pair<int, int> temp2(std::tr1::get<1>(associationList[i]),std::tr1::get<2>(associationList[i]));
                particleList[std::tr1::get<0>(associationList[i])].complexList[0] = temp2;
                //particleList[std::tr1::get<0>(associationList[i])].complexList.push_back(temp2);
                
                //if D particle is formed
                if(particleList[std::tr1::get<1>(associationList[i])].complexList[0].first != -1 &&
                   particleList[std::tr1::get<1>(associationList[i])].complexList[1].first != -1)
                {
                    std::cout << "double complex formed... " << std::endl;
                    std::cout <<particleList[std::tr1::get<1>(associationList[i])].complexList[0].first << " " << particleList[std::tr1::get<1>(associationList[i])].complexList[1].first << std::endl;
                    //draw next event times (dissociation times)
                    double nextDissociationTime = (1.0/(2.0*blockedDissociationRate)) * (-log(gsl_rng_uniform(rng)));
                    double nextPhosphorylationTime =  (1.0/phosphorylationRate) * (-log(gsl_rng_uniform(rng)));
                    
                    if(nextDissociationTime < nextPhosphorylationTime )
                    {
                        particleList[std::tr1::get<1>(associationList[i])].nextEventTime = nextDissociationTime;
                        particleList[std::tr1::get<1>(associationList[i])].nextEventTime += tSim;
                        particleList[std::tr1::get<1>(associationList[i])].nextReactionTime = nextDissociationTime;
                        particleList[std::tr1::get<1>(associationList[i])].nextReactionTime += tSim;
                        particleList[std::tr1::get<1>(associationList[i])].nextReactionType = 1;
                        particleList[std::tr1::get<1>(associationList[i])].nextEventType = 1;
                    }
                    else if( nextPhosphorylationTime < nextDissociationTime )
                    {
                        particleList[std::tr1::get<1>(associationList[i])].nextEventTime = nextPhosphorylationTime;
                        particleList[std::tr1::get<1>(associationList[i])].nextEventTime += tSim;
                        particleList[std::tr1::get<1>(associationList[i])].nextReactionTime = nextPhosphorylationTime;
                        particleList[std::tr1::get<1>(associationList[i])].nextReactionTime += tSim;
                        particleList[std::tr1::get<1>(associationList[i])].nextReactionType = 3;
                        particleList[std::tr1::get<1>(associationList[i])].nextEventType = 1;
                    }
                    ///add the substrate particle to the scheduler list
                    //schedulerList.push_back(std::tr1::get<1>(associationList[i]));
                }
                ///if its a c particle
                else if((particleList[std::tr1::get<1>(associationList[i])].complexList[0].first == -1 &&
                         particleList[std::tr1::get<1>(associationList[i])].complexList[1].first != -1)
                        ||
                        (particleList[std::tr1::get<1>(associationList[i])].complexList[1].first == -1 &&
                         particleList[std::tr1::get<1>(associationList[i])].complexList[0].first != -1)
                       )
                {
                    std::cout << "single complex formed... " << std::endl;
                    //draw next event times (dissociation times)
                    double nextDissociationTime     = (1.0/dissociationRate) * (-log(gsl_rng_uniform(rng)));
                    double nextPhosphorylationTime  = (1.0/phosphorylationRate) * (-log(gsl_rng_uniform(rng)));
                    double nextHoppingTime          = (1.0/hoppingRate) * (-log(gsl_rng_uniform(rng)));
                    
                    
                    if( (nextDissociationTime < nextPhosphorylationTime) && (nextDissociationTime < nextHoppingTime) )
                    {
                        particleList[std::tr1::get<1>(associationList[i])].nextEventTime = nextDissociationTime;
                        particleList[std::tr1::get<1>(associationList[i])].nextEventTime += tSim;
                        particleList[std::tr1::get<1>(associationList[i])].nextReactionTime = nextDissociationTime;
                        particleList[std::tr1::get<1>(associationList[i])].nextReactionTime += tSim;
                        particleList[std::tr1::get<1>(associationList[i])].nextReactionType = 1;
                        particleList[std::tr1::get<1>(associationList[i])].nextEventType = 1;
                    }
                    else if( (nextPhosphorylationTime < nextDissociationTime) && (nextPhosphorylationTime < nextHoppingTime) )
                    {
                        particleList[std::tr1::get<1>(associationList[i])].nextEventTime = nextPhosphorylationTime;
                        particleList[std::tr1::get<1>(associationList[i])].nextEventTime += tSim;
                        particleList[std::tr1::get<1>(associationList[i])].nextEventType = 1;
                        particleList[std::tr1::get<1>(associationList[i])].nextReactionTime = nextPhosphorylationTime;
                        particleList[std::tr1::get<1>(associationList[i])].nextReactionTime += tSim;
                        particleList[std::tr1::get<1>(associationList[i])].nextReactionType = 3;
                    }
                    else if ( (nextHoppingTime < nextPhosphorylationTime) && (nextHoppingTime < nextDissociationTime) )
                    {
                        particleList[std::tr1::get<1>(associationList[i])].nextEventTime = nextHoppingTime;
                        particleList[std::tr1::get<1>(associationList[i])].nextEventTime += tSim;
                        particleList[std::tr1::get<1>(associationList[i])].nextEventType = 1;
                        particleList[std::tr1::get<1>(associationList[i])].nextReactionTime = nextHoppingTime;
                        particleList[std::tr1::get<1>(associationList[i])].nextReactionTime += tSim;
                        particleList[std::tr1::get<1>(associationList[i])].nextReactionType = 2;
                    }
                    ///add the substrate particle to the scheduler list
                    //schedulerList.push_back(std::tr1::get<1>(associationList[i]));
                }
                
                //change the nextEventTime of the enzyme
                particleList[std::tr1::get<0>(associationList[i])].nextEventTime = particleList[std::tr1::get<1>(associationList[i])].nextReactionTime;
                particleList[std::tr1::get<0>(associationList[i])].nextEventType = 1;
                
                
                std::cout << " The next event time and type and the next reaction time and type : "<< particleList[std::tr1::get<1>(associationList[i])].nextEventTime << " " << particleList[std::tr1::get<1>(associationList[i])].nextEventType << " " << particleList[std::tr1::get<1>(associationList[i])].nextReactionTime << " " << particleList[std::tr1::get<1>(associationList[i])].nextReactionType << std::endl;
                
            }
        }
        std::cout << "********************************************************************************************" <<std::endl;
        associationList.clear();
    }
    
}
