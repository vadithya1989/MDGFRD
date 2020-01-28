/*
 * File:   particle.cpp
 * Author: adithya
 *
 * Created on December 15, 2014, 7:23 PM
 */


#include "particle.h"
#include "eulerAngles.h"
#include "species.h"
#include "greensFunction.h"
#include "simulationParameters.h"


particle::particle()
{
    particleID = int(0);
    speciesType  = int(0);
    particleType = int(0);
    
    particlePosition.x = particlePosition.y = particlePosition.z = double(0);
    force.x = force.y = force.z = double(0);
    
    q.q0 = double(1);
    q.q1 = q.q2 = q.q3 = double(0);
    torque.q0 = torque.q1 = torque.q2 = torque.q3 = double(0);
    
    //patchPos.x = double(1);
    //patchPos.y = patchPos.z = double(0);
    
    std::vector<myVector> tempPP;
    patchPos = tempPP;
    
    nextEscapePosition.x = nextEscapePosition.y = nextEscapePosition.z = double(0);
    burstPosition.x = burstPosition.y = burstPosition.z = double(0);
    R = double(0);
    
    nextEventTime = double(0);
    nextEscapeTime = double(0);
    nextReactionTime = double(0);
    nextReactionType = int(0);
    nextEventType = int(0);
    timeOfConstructionOfDomain = double(0);
    domainRadius = double(0);
    
    //std::pair<int,int> tempPair(-1,-1);
    //std::vector<std::pair<int,int>> temp(2,tempPair);
    std::vector<std::pair<int,int>> temp;
    complexList = temp;
    
    std::vector<int> temp1;
    patchState = temp1;
}

particle::particle(const particle& orig)
{
    particleID = orig.particleID;
    speciesType = orig.speciesType;
    particleType = orig.particleType;
    
    particlePosition = orig.particlePosition;
    force = orig.force;
    
    q = orig.q;
    torque = orig.torque;
    
    patchPos = orig.patchPos;
    
    nextEscapePosition = orig.nextEscapePosition;
    burstPosition = orig.burstPosition;
    R = orig.R;
    
    nextEventTime = orig.nextEventTime;
    nextEscapeTime = orig.nextEscapeTime;
    nextReactionTime = orig.nextReactionTime;
    nextReactionType = orig.nextReactionType;
    nextEventType = orig.nextEventType;
    timeOfConstructionOfDomain = orig.timeOfConstructionOfDomain;
    domainRadius = orig.domainRadius;
    
    patchState = orig.patchState;
}

particle::~particle()
{
    
}

std::ostream& operator<< (std::ostream &out, particle &part)
{
    // Since operator<< is a friend of the Point class, we can access
    // Point's members directly.
    out <<          part.speciesType << " " <<
    part.particlePosition.x <<" " <<
    part.particlePosition.y <<" " <<
    part.particlePosition.z <<" " <<
    part.force.x    <<" " <<
    part.force.y    <<" " <<
    part.force.z    <<" " <<
    part.q.q0       <<" " <<
    part.q.q1       <<" " <<
    part.q.q2       <<" " <<
    part.q.q3       <<" " <<
    part.torque.q0  <<" " <<
    part.torque.q1  <<" " <<
    part.torque.q2  <<" " <<
    part.torque.q3  <<" " << std::endl;
    //part.patchPos.x <<" " <<
    //part.patchPos.y <<" " <<
    //part.patchPos.z <<" " << std::endl;
    return out;
}

void particle::drawNextEscapePosition()
{
    gsl_ran_dir_3d (rng, &nextEscapePosition.x, &nextEscapePosition.y, &nextEscapePosition.z);
    nextEscapePosition = nextEscapePosition * domainRadius;
    //update orientation
    eulerAngles ea0, ea;
    double timeInDomain = (nextEventTime-timeOfConstructionOfDomain);
    if(speciesList[speciesType].rotationalDiffusionConstant*timeInDomain < 0.05)
    {
        myVector axis;
        gsl_ran_dir_3d (rng, &axis.x, &axis.y, &axis.z);
        
        double SD = sqrt(2.0*speciesList[speciesType].rotationalDiffusionConstant*timeInDomain);
        myVector phi;
        phi.x = gsl_ran_gaussian ( rng, SD );
        phi.y = gsl_ran_gaussian ( rng, SD );
        phi.z = gsl_ran_gaussian ( rng, SD );
        
        double phiMag = phi.magnitude();
        
        q.rotate(phiMag, axis);
    }
    
    else
    {
        ea = drawAngles(ea0, speciesList[speciesType].rotationalDiffusionConstant, timeInDomain, rng);
        drawOrientation(*this, ea);
    }
}

void particle::drawBurstPosition(double time)
{
    greensFunction gf(speciesList[speciesType].diffusionConstant, domainRadius);
    double tempRNG =gsl_rng_uniform(rng);
    //std::cout <<"tempRNG " << time << std::endl;
     R = gf.drawR(tempRNG, time);
    gsl_ran_dir_3d (rng, &burstPosition.x, &burstPosition.y, &burstPosition.z);
    burstPosition = burstPosition * R;
    //update orientation
    eulerAngles ea0, ea;
    double timeInDomain = time;
    if(speciesList[speciesType].rotationalDiffusionConstant*timeInDomain < 0.05)
    {
        myVector axis;
        gsl_ran_dir_3d (rng, &axis.x, &axis.y, &axis.z);
        
        double SD = sqrt(2.0*speciesList[speciesType].rotationalDiffusionConstant*timeInDomain);
        myVector phi;
        phi.x = gsl_ran_gaussian ( rng, SD );
        phi.y = gsl_ran_gaussian ( rng, SD );
        phi.z = gsl_ran_gaussian ( rng, SD );
        
        double phiMag = phi.magnitude();
        
        q.rotate(phiMag, axis);
    }
    
    else
    {
        ea = drawAngles(ea0, speciesList[speciesType].rotationalDiffusionConstant, timeInDomain, rng);
        drawOrientation(*this, ea);
    }
}

void particle::drawNextEventTimeAndType()
{
    greensFunction gf(speciesList[speciesType].diffusionConstant, domainRadius);
    //draw time only returns the next escape time
    std::pair<double, int> temp = gf.drawTime(rng, speciesList, speciesType );
    nextEscapeTime = temp.first;
    nextEscapeTime += tSim;
}


void particle::buildADomain(double dist)
{
    //std::cout << "domain Built with ID: " << particleID << " and species type " << speciesType << std::endl;
    particleType = 1;
    domainRadius = dist;
    timeOfConstructionOfDomain = tSim;
    drawNextEventTimeAndType();
    //if escape happens first NET is escape
   
    if( nextEscapeTime < nextReactionTime)
    {
        nextEventTime = nextEscapeTime;
        //nextEventTime += tSim;
        nextEventType = 0;
    }
    //if dissociation happens first NET is escape
    if( nextReactionTime < nextEscapeTime)
    {
        nextEventTime = nextReactionTime;
        //nextEventTime += tSim;
        nextEventType = 1;
    }
    
}

void particle::burstTheDomain(std::string condition)
{
    //////////UPdate complex list position!!!!!!!!!!!!!!!!!!!!!!!!!!!
    particleType = 2;
    if(condition.compare("escape") == 0)
    {
        drawNextEscapePosition();
        particlePosition = particlePosition + nextEscapePosition;
        for(int i=0;i<complexList.size();++i)
        {
            //update positions of all particles bound to the substrate
            myVector patchDirection =speciesList[speciesType].bodyPatchVectors[complexList[i].second];
            particleList[complexList[i].first].particlePosition = q.bodyToSpace(patchDirection)* 2.0*speciesList[speciesType].particleRadius + particlePosition;
        }
    }
    if(condition.compare("burst") == 0)
    {
        //std::cout << particleID << " " <<tSim-timeOfConstructionOfDomain << std::endl;
        drawBurstPosition(std::abs(tSim-timeOfConstructionOfDomain));
        particlePosition = particlePosition + burstPosition;
        for(int i=0;i<complexList.size();++i)
        {
            myVector patchDirection =speciesList[speciesType].bodyPatchVectors[complexList[i].second];
            particleList[complexList[i].first].particlePosition = q.bodyToSpace(patchDirection)* 2.0*speciesList[speciesType].particleRadius + particlePosition;
        }
    }
    if(condition.compare("convert") == 0)
    {
        particlePosition = particlePosition;
        for(int i=0;i<complexList.size();++i)
        {
            myVector patchDirection =speciesList[speciesType].bodyPatchVectors[complexList[i].second];
            particleList[complexList[i].first].particlePosition = q.bodyToSpace(patchDirection)* 2.0*speciesList[speciesType].particleRadius + particlePosition;
        }
    }
    
    domainRadius = 0.0;/*the radius of the domain is the radius of the particle*/
    timeOfConstructionOfDomain = tSim;
    
    //the next event time and type is dissociation with the time drawn before
    nextEscapeTime = nextReactionTime;
    nextEventType = 1;/*for now since there is no domain around the particle, next event type is decay*/
    nextEventTime += tSim;/*to bring the next event time to the scale of the global simulation*/
    
}

void particle::dissociateParticleFromDomain()
{
    std::cout << "***************************************GFRD DISSOCIATION*********************************" << std::endl;
    ///first step is to burst the domain and get the new position and orientation
    burstTheDomain("burst");
    
    
    ///++++++++++++++++++++++++++++++++If the dissociating particle is a D particle++++++++++++++++++++++++++++++++++++++++++//
    if(complexList[0].first != -1 && complexList[1].first!=-1)
    {
        std::cout << "D particle GFRD dissociation" << std::endl;
        ///Draw a random number 0 or 1 to randomly dissociate a particle
        int n=1;
        u_long patchIndex = gsl_rng_uniform_int (rng, n);
        int otherPatchIndex=0;
        if(patchIndex == 0)
            otherPatchIndex = 1;
        
        ///get the ids of the enzyme that are dissociating
        int dissociatingEnzyme = (complexList[patchIndex]).first;
        
        //remove the enzyme from the idle particle list
        idleParticleList.erase( std::remove(idleParticleList.begin(), idleParticleList.end(), dissociatingEnzyme), idleParticleList.end() );
        
        //+++++++++++++++++++++++++If the reaction is dissociation+++++++++++++++++++++++++++++++//
        if (nextReactionType == 1)
        {
            ///For now the position of the dissociated particle is put just at a distance of 2Sigma
            myVector patchUnitVector = q.bodyToSpace(speciesList[speciesType].bodyPatchVectors[patchIndex])/ q.bodyToSpace(speciesList[speciesType].bodyPatchVectors[patchIndex]).magnitude();
            myVector newPosition = patchUnitVector*(2.0*prm.sigma);
            
            particleList[dissociatingEnzyme].particlePosition = newPosition;
            ///BD particle till domain is built
            particleList[dissociatingEnzyme].particleType = 2;
        }
        //+++++++++++++++++++++++++If the reaction is phosphorylation+++++++++++++++++++++++++++++++//
        else if(nextReactionType == 3)
        {
            //If the unbinding enzyme is a kinase then the substrate is phosphorylated
            if(particleList[dissociatingEnzyme].speciesType == 3)
            {
                //change the state of the patch to phosphorylated
                //patchState[particleList[dissociatingEnzyme].complexList[patchIndex].second] = 1;
                patchState[patchIndex] = 1;
                //check the other patch and then decide the species type of the substrate
                if( patchState[otherPatchIndex] == 1)
                    speciesType = 2;
                else if (patchState[otherPatchIndex] == 0)
                    speciesType = 1;
                
                //++++++++++++++++++++++++++++++++COUNTER++++++++++++++++++++++++++++++++++++++++++++//
                //if the other patch is also phosphorylated then turn ON D switch
                if(patchState[otherPatchIndex] == 1)
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
                else if(patchState[otherPatchIndex] == 0)
                {
                    ///counter for singly phosphorylated (C particle) C increases
                    tBoundEndC = tSim;
                    tBoundC = tBoundC + numberOfCParticles * (tBoundEndC - tBoundStartC);
                    probabilityOfBeingBoundC = tBoundC/(maxCparticles*tSim);
                    std::cout << " C " << numberOfCParticles << " " << (tBoundEndC - tBoundStartC) << " " << tSim << " " << probabilityOfBeingBoundC << std::endl;
                    tBoundStartC = tSim;
                    numberOfCParticles = numberOfCParticles + 1;
                }
                //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
            }
            //If the unbinding enzyme is a phosphotase then the substrate is dephosphorylated
            else if(particleList[dissociatingEnzyme].speciesType == 4)
            {
                //change the state of the patch to dephosphorylated
                patchState[patchIndex] = 0;
                
                //check the other patch and then decide the species type of the substrate
                if( patchState[otherPatchIndex] == 1)
                    speciesType = 1;
                else if (patchState[otherPatchIndex] == 0)
                    speciesType = 0;
                
                //++++++++++++++++++++++++++++++++COUNTER++++++++++++++++++++++++++++++++++++++++++++//
                //if the other patch is phosphorylated then turn ON D switch
                if(patchState[otherPatchIndex] == 1)
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
                else if(patchState[otherPatchIndex] == 0)
                {
                    ///counter for singly phosphorylated (C particle) C increases
                    tBoundEndC = tSim;
                    tBoundC = tBoundC + numberOfCParticles * (tBoundEndC - tBoundStartC);
                    probabilityOfBeingBoundC = tBoundC/(maxCparticles*tSim);
                    std::cout << " C "<< numberOfCParticles << " " << (tBoundEndC - tBoundStartC) << " " << tSim << " " << probabilityOfBeingBoundC << std::endl;
                    tBoundStartC = tSim;
                    numberOfCParticles = numberOfCParticles - 1;
                }
                //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
            }
        }
        
        ///draw next event times (dissociation times)
        double nextDissociationTime     = (1.0/dissociationRate) * (-log(gsl_rng_uniform(rng)));
        double nextPhosphorylationTime  = (1.0/phosphorylationRate) * (-log(gsl_rng_uniform(rng)));
        double nextHoppingTime          = (1.0/hoppingRate) * (-log(gsl_rng_uniform(rng)));
        
        ///the minimum of the above times is the next evebt time
        if( (nextDissociationTime < nextPhosphorylationTime) && (nextDissociationTime < nextHoppingTime) )
        {
            nextEventTime = nextDissociationTime;
            nextEventTime += tSim;
            nextReactionTime = nextDissociationTime;
            nextReactionTime += tSim;
            nextReactionType = 1;
            nextEventType = 1;
        }
        else if( (nextPhosphorylationTime < nextDissociationTime) && (nextPhosphorylationTime < nextHoppingTime) )
        {
            nextEventTime = nextPhosphorylationTime;
            nextEventTime += tSim;
            nextReactionTime = nextPhosphorylationTime;
            nextReactionTime += tSim;
            nextReactionType = 3;
            nextEventType = 1;
        }
        else if ( (nextHoppingTime < nextPhosphorylationTime) && (nextHoppingTime < nextDissociationTime) )
        {
            nextEventTime = nextHoppingTime;
            nextEventTime += tSim;
            nextReactionTime = nextHoppingTime;
            nextReactionTime += tSim;
            nextReactionType = 2;
            nextEventType = 1;
        }
        
        //change complex list of the substrate and enzyme
        std::pair<int, int>  tempPair(-1,-1);
        complexList[particleList[dissociatingEnzyme].complexList[0].second] = tempPair;
        particleList[dissociatingEnzyme].complexList[0] = tempPair;
        
        ///Remove the dissociated enzyme from the complex list of the substrate and viceversa
        //particleList[dissociatingEnzyme].complexList.pop_back();
        //complexList.erase(complexList.begin()+patchIndex);
    }
    
    
    
    ///++++++++++++++++++++++++++++++++If the dissociating particle is a C particle++++++++++++++++++++++++++++++++++++++++++//
    else if((complexList[0].first == -1 &&
             complexList[1].first != -1)
            ||
            (complexList[1].first == -1 &&
             complexList[0].first != -1)
            )
    {
        
        std::cout << "C particle GFRD dissociation "<< std::endl;
        //get the ids of the dissociating particles
        int dissociatingEnzyme=0;
        int patchIndex;
        int otherPatchIndex;
        if(complexList[0].first == -1 && complexList[1].first != -1)
        {
            dissociatingEnzyme = (complexList[1]).first;
        }
        else if(complexList[1].first == -1 && complexList[0].first != -1)
        {
            dissociatingEnzyme = (complexList[0]).first;
        }
        
        patchIndex = particleList[dissociatingEnzyme].complexList[0].second;
        
        otherPatchIndex=0;
        if(patchIndex == 0)
            otherPatchIndex = 1;
        
        
        
        
        //std::cout << " patch INdex " << patchIndex << " " << otherPatchIndex << std::endl;
        //remove the enzyme from the idle particle list
        idleParticleList.erase( std::remove(idleParticleList.begin(), idleParticleList.end(), dissociatingEnzyme), idleParticleList.end() );
        
        //+++++++++++++++++++++++++If the reaction is dissociation+++++++++++++++++++++++++++++++//
        if (nextReactionType == 1)
        {
            ///For now the position of the dissociated particle is put just at a distance of 2Sigma
            myVector patchUnitVector = q.bodyToSpace(speciesList[speciesType].bodyPatchVectors[patchIndex])/ q.bodyToSpace(speciesList[speciesType].bodyPatchVectors[patchIndex]).magnitude();
            myVector newPosition = patchUnitVector*(2.0*prm.sigma);
            
            particleList[dissociatingEnzyme].particlePosition = newPosition;
            ///BD particle till domain is built
            particleList[dissociatingEnzyme].particleType = 2;
        }
        //+++++++++++++++++++++++++If the reaction is hopping+++++++++++++++++++++++++++++++//
        else if(nextReactionType == 2)
        {
            ///put enzyme to other patch
            myVector patchUnitVector = q.bodyToSpace(speciesList[speciesType].bodyPatchVectors[otherPatchIndex])/ q.bodyToSpace(speciesList[speciesType].bodyPatchVectors[patchIndex]).magnitude();
            myVector newPosition = patchUnitVector*(2.0*prm.sigma);
            
            particleList[dissociatingEnzyme].particlePosition = newPosition;
            ///BD particle till domain is built
            particleList[dissociatingEnzyme].particleType = 2;
        }

        //+++++++++++++++++++++++++If the reaction is phosphorylation+++++++++++++++++++++++++++++++//
        else if(particleList[schedulerList[0]].nextReactionType == 3)
        {
            
            if(particleList[dissociatingEnzyme].speciesType == 3)
            {
                
                patchState[patchIndex] = 1;
                
                if( patchState[otherPatchIndex] == 1)
                    speciesType = 2;
                else if (patchState[otherPatchIndex] == 0)
                    speciesType = 1;
                //if the other patch is also phosphorylated then turn ON D switch
                if(patchState[otherPatchIndex] == 1)
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
                else if(patchState[otherPatchIndex] == 0)
                {
                    
                    ///counter for singly phosphorylated (C particle) C increases
                    tBoundEndC = tSim;
                    tBoundC = tBoundC + numberOfCParticles * (tBoundEndC - tBoundStartC);
                    probabilityOfBeingBoundC = tBoundC/(maxCparticles*tSim);
                    std::cout << "C " << numberOfCParticles << " " << (tBoundEndC - tBoundStartC) << " " << tSim << " " << probabilityOfBeingBoundC << std::endl;
                    tBoundStartC = tSim;
                    numberOfCParticles = numberOfCParticles + 1;
                }
                
            }
            else if(particleList[dissociatingEnzyme].speciesType == 4)
            {
                
                patchState[patchIndex] = 0;
                if( patchState[otherPatchIndex] == 1)
                    speciesType = 1;
                else if (patchState[otherPatchIndex] == 0)
                    speciesType = 0;
                
                //if the other patch is phosphorylated then turn ON D switch
                if(patchState[otherPatchIndex] == 1)
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
                else if(patchState[otherPatchIndex] == 0)
                {
                    
                    ///counter for singly phosphorylated (C particle) C increases
                    tBoundEndC = tSim;
                    tBoundC = tBoundC + numberOfCParticles * (tBoundEndC - tBoundStartC);
                    probabilityOfBeingBoundC = tBoundC/(maxCparticles*tSim);
                    std::cout<<"C " << numberOfCParticles << " " << (tBoundEndC - tBoundStartC) << " " << tSim << " " << probabilityOfBeingBoundC << std::endl;
                    tBoundStartC = tSim;
                    numberOfCParticles = numberOfCParticles - 1;
                }
            }
        }
        ///draw next event times (dissociation times)
        double nextDissociationTime     = (1.0/singleParticleDissociationRate) * (-log(gsl_rng_uniform(rng)));
        
        nextEventTime = nextDissociationTime;
        nextEventTime += tSim;
        nextReactionTime = nextDissociationTime;
        nextReactionTime += tSim;
        nextReactionType = 0;
        nextEventType = 1;
        
        //change complex list of the substrate and enzyme
        std::pair<int, int>  tempPair(-1,-1);
        complexList[particleList[dissociatingEnzyme].complexList[0].second] = tempPair;
        particleList[dissociatingEnzyme].complexList[0] = tempPair;
        
        
        ///Remove the dissociated enzyme from the complex list of the substrate and viceversa
        //complexList.pop_back();
        //particleList[dissociatingEnzyme].complexList.pop_back();
        
    }
    /*sort the particle list again based on the next event times*/
    sortSchedulerList();
    std::cout << "*********************************************************************************************" << std::endl;
}



