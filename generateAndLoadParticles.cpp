
//#include "simulationParameters.h"
//#include "species.h"
//#include "particle.h"
//#include <gsl/gsl_rng.h>
//#include <gsl/gsl_randist.h>
//#include <cmath>
//#include <iostream>
//#include <fstream>
//#include <sstream>


#include "main.h"

/**
 * This generates a random co-ordinate and accepts it only if there is no overlap
 * @param i species type
 * @param particleNumber is the particle id
 */
void getCoordinate( int i,
                    int particleNumber)
{
    /*GENERATE THE PARTICLE*/
    double rangeStart = 0.3*prm.xmin;
    double rangeEnd = 0.7*prm.xmax;
    double x = gsl_ran_flat (rng, rangeStart, rangeEnd);
    double y = gsl_ran_flat (rng, rangeStart, rangeEnd);
    double z = gsl_ran_flat (rng, rangeStart, rangeEnd);
    myVector coOrd;
    double dist = 0.0;
    coOrd.x = x; coOrd.y = y; coOrd.z = z;
    /*GET THE DISTANCE OF THE PARTICLE TO ITS NEIGHBORS*/
    std::vector<double> distances;
    for(int k=0; k < particleList.size(); ++k)
    {
        myVector particlePos = particleList[k].particlePosition;
        double distSq = (particlePos.x-x)*(particlePos.x-x) + (particlePos.y-y)*
            (particlePos.y-y) + (particlePos.z-z)*(particlePos.z-z);
        dist = sqrt(distSq);
        distances.push_back(dist);
    }
    /*CHECK FOR OVERLAPS*/
    /*if pass = 1, then the co-ordinate generated does not overlap with any neighboring particle, else it does and a new co-ordinate is generated*/
    int pass = 1; //int neighborNum = 0;
    for(std::vector<double>::iterator it=distances.begin(); it!=distances.end(); ++it)
    {
        if(*it < ( 1.2*pow(2,0.333333333333333)*prm.sigma))
        {
          //std::cout << "particles overlap" << std::endl;
          pass = 0 ;
          break;
        }
    }
    if(pass == 1)
    {
       
        particle temp;
        
        temp.particleID = particleNumber;
        temp.speciesType = i;
        temp.particleType = 2;
        temp.particlePosition = coOrd;
        temp.q.q0=1;
        temp.q.q1=0;
        temp.q.q2=0;
        temp.q.q3=0;
        
        temp.domainRadius = 0;/*the radius of the domain is the radius of the particle*/
        temp.nextEventTime = (1.0/singleParticleDissociationRate) * (-log(gsl_rng_uniform(rng)));/*Draw the next dissociation time of the particle based on the dissociation constant, which is got by an independent simulation and fitting the first escape times with an exponential function*/
        temp.nextEventTime += 0.0;/*to bring the next event time to the scale of the global simulation*/
        temp.nextReactionTime = temp.nextEventTime;/*reaction time is NET*/
        temp.nextEventType = 1;/*for now since there is no domain around the particle, next event type is decay*/
        temp.nextReactionType = 1;/*NRType is dissociation*/
        
        //temp.patchPos = temp.q.bodyToSpace() * speciesList[temp.speciesType].particleRadius + temp.particlePosition;
        
        ///if the particle is a substrate, then update patch type to dephosd. and if enzyme its active
        if(temp.speciesType == 0 || temp.speciesType == 1 || temp.speciesType == 2)
        {
            for(int i=0; i<speciesList[temp.speciesType].numberOfPatches;++i)
            {
                temp.patchState.push_back(0);
            }
            
        }
        else
        {
            temp.patchState.push_back(1);
        }
        
        for(int i=0;i<speciesList[temp.speciesType].numberOfPatches;++i)
        {
            //complexList
            std::pair<int,int> tempComplexList;
            tempComplexList = std::make_pair(-1, -1);
            temp.complexList.push_back(tempComplexList);
            //patchPos
            myVector patchDirection =speciesList[temp.speciesType].bodyPatchVectors[i];
            myVector tempPatchPos = temp.q.bodyToSpace(patchDirection) * speciesList[temp.speciesType].particleRadius + temp.particlePosition;
            temp.patchPos.push_back(tempPatchPos);
        }
        

        //add the particle into the particleList
        particleList.push_back(temp);
        //add the particle to the scheduler list
        schedulerList.push_back(particleNumber);
        

    }
    else
    {
        getCoordinate(i, particleNumber);
    }
}


void generateAndLoadParticles(int &readin)
{
    if(readin == 0)
    {
        int particleNumber=0;
        for (int i = 0; i < prm.totalNumberOfSpecies; ++i)
        {
            for (int j=0; j < speciesList[i].numberOfParticles; ++j)
            {
                getCoordinate(i, particleNumber);
                ++particleNumber;
            }
        }
        
        for(int k=0;k<particleList.size();++k)
        {
            for(int p=0; p<speciesList[particleList[k].speciesType].numberOfPatches;++p)
            {
                std::pair<int,int> temp1ComplexList(-1,-1);
                particleList[k].complexList.push_back(temp1ComplexList);
            }
        }
    }
    
    else if (readin == 1)
    {
        
        
        particle temp;
        
        
        temp.particleID = 0;
        temp.speciesType = 0;
        temp.particleType = 2;
        temp.particlePosition.x = prm.xmax/2;
        temp.particlePosition.y = prm.ymax/2;
        temp.particlePosition.z = prm.zmax/2;
        temp.q.q0=1;
        temp.q.q1=0;
        temp.q.q2=0;
        temp.q.q3=0;
        
        temp.domainRadius = 0;/*the radius of the domain is the radius of the particle*/
        temp.nextEventTime = (1.0/singleParticleDissociationRate) * (-log(gsl_rng_uniform(rng)));/*Draw the next dissociation time of the particle based on the dissociation constant, which is got by an independent simulation and fitting the first escape times with an exponential function*/
        temp.nextEventTime += 0.0;/*to bring the next event time to the scale of the global simulation*/
        temp.nextReactionTime = temp.nextEventTime;/*reaction time is NET*/
        temp.nextEventType = 1;/*for now since there is no domain around the particle, next event type is decay*/
        temp.nextReactionType = 1;/*NRType is dissociation*/
        
        //temp.patchPos = temp.q.bodyToSpace() * speciesList[temp.speciesType].particleRadius + temp.particlePosition;
        
        ///if the particle is a substrate, then update patch type to dephosd. and if enzyme its active
        if(temp.speciesType == 0 || temp.speciesType == 1 || temp.speciesType == 2)
        {
            
            for(int i=0; i<speciesList[temp.speciesType].numberOfPatches;++i)
            {
                temp.patchState.push_back(0);
                //temp.patchState.push_back(0);
            }
            
        }
        else
        {
            temp.patchState.push_back(1);
        }
        
        
        for(int i=0;i<speciesList[temp.speciesType].numberOfPatches;++i)
        {
            
            //complexList
            
            std::pair<int,int> tempComplexList;
            tempComplexList = std::make_pair(-1, -1);
            temp.complexList.push_back(tempComplexList);
            
            //patchPos
            myVector patchDirection =speciesList[temp.speciesType].bodyPatchVectors[i];
            myVector tempPatchPos = temp.q.bodyToSpace(patchDirection) * speciesList[temp.speciesType].particleRadius + temp.particlePosition;
            temp.patchPos.push_back(tempPatchPos);
        }
        
        //add the particle into the particleList
        particleList.push_back(temp);
        //add the particle to the scheduler list
        schedulerList.push_back(temp.particleID);
        
        
        particle temp1;
        
        temp1.particleID = 1;
        temp1.speciesType = 3;
        temp1.particleType = 2;
        temp1.particlePosition.x = prm.xmax/2+(1.001799126334586829134*5e-9);//(pow(2,0.3333333)*prm.sigma);
        temp1.particlePosition.y = prm.ymax/2;
        temp1.particlePosition.z = prm.zmax/2;
        temp1.q.q0=1;
        temp1.q.q1=0;
        temp1.q.q2=0;
        temp1.q.q3=0;
        
        temp1.domainRadius = 0;/*the radius of the domain is the radius of the particle*/
        temp1.nextEventTime = (1.0/singleParticleDissociationRate) * (-log(gsl_rng_uniform(rng)));/*Draw the next dissociation time of the particle based on the dissociation constant, which is got by an independent simulation and fitting the first escape times with an exponential function*/
        temp1.nextEventTime += 0.0;/*to bring the next event time to the scale of the global simulation*/
        temp1.nextReactionTime = temp1.nextEventTime;/*reaction time is NET*/
        temp1.nextEventType = 1;/*for now since there is no domain around the particle, next event type is decay*/
        temp1.nextReactionType = 1;/*NRType is dissociation*/
        
        //temp.patchPos = temp.q.bodyToSpace() * speciesList[temp.speciesType].particleRadius + temp.particlePosition;
        
        ///if the particle is a substrate, then update patch type to dephosd. and if enzyme its active
        if(temp1.speciesType == 0 || temp1.speciesType == 1 || temp1.speciesType == 2)
        {
            for(int i=0; i<speciesList[temp1.speciesType].numberOfPatches;++i)
            {
                temp1.patchState.push_back(0);
            }
            
        }
        else
        {
            temp1.patchState.push_back(1);
        }
        
        for(int i=0;i<speciesList[temp1.speciesType].numberOfPatches;++i)
        {
            //complexList
            std::pair<int,int> temp1ComplexList(-1,-1);
            temp1.complexList.push_back(temp1ComplexList);
            
            myVector patchDirection =speciesList[temp1.speciesType].bodyPatchVectors[i];
            myVector temp1PatchPos = temp1.q.bodyToSpace(patchDirection) * speciesList[temp1.speciesType].particleRadius + temp1.particlePosition;
            temp1.patchPos.push_back(temp1PatchPos);
        }
        
        
        //add the particle into the particleList
        particleList.push_back(temp1);
        //add the particle to the scheduler list
        schedulerList.push_back(temp1.particleID);
        
        
        particle temp2;
        
        temp2.particleID = 2;
        temp2.speciesType = 4;
        temp2.particleType = 2;
        temp2.particlePosition.x = prm.xmax/1.2;//+(1.001799126334586829134*5e-9);//(pow(2,0.3333333)*prm.sigma);
        temp2.particlePosition.y = prm.ymax/1.2;
        temp2.particlePosition.z = prm.zmax/2;
        temp2.q.q0=1;
        temp2.q.q1=0;
        temp2.q.q2=0;
        temp2.q.q3=0;
        
        temp2.domainRadius = 0;/*the radius of the domain is the radius of the particle*/
        temp2.nextEventTime = (1.0/singleParticleDissociationRate) * (-log(gsl_rng_uniform(rng)));/*Draw the next dissociation time of the particle based on the dissociation constant, which is got by an independent simulation and fitting the first escape times with an exponential function*/
        temp2.nextEventTime += 0.0;/*to bring the next event time to the scale of the global simulation*/
        temp2.nextReactionTime = temp2.nextEventTime;/*reaction time is NET*/
        temp2.nextEventType = 1;/*for now since there is no domain around the particle, next event type is decay*/
        temp2.nextReactionType = 1;/*NRType is dissociation*/
        
        //temp.patchPos = temp.q.bodyToSpace() * speciesList[temp.speciesType].particleRadius + temp.particlePosition;
        
        ///if the particle is a substrate, then update patch type to dephosd. and if enzyme its active
        if(temp2.speciesType == 0 || temp2.speciesType == 1 || temp2.speciesType == 2)
        {
            for(int i=0; i<speciesList[temp2.speciesType].numberOfPatches;++i)
            {
                temp2.patchState.push_back(0);
            }
            
        }
        else
        {
            temp2.patchState.push_back(1);
        }
        
        for(int i=0;i<speciesList[temp2.speciesType].numberOfPatches;++i)
        {
            //complexList
            std::pair<int,int> temp2ComplexList(-1,-1);
            temp2.complexList.push_back(temp2ComplexList);
            
            myVector patchDirection =speciesList[temp2.speciesType].bodyPatchVectors[i];
            myVector temp2PatchPos = temp2.q.bodyToSpace(patchDirection) * speciesList[temp2.speciesType].particleRadius + temp2.particlePosition;
            temp2.patchPos.push_back(temp2PatchPos);
        }
        
        
        //add the particle into the particleList
        particleList.push_back(temp2);
        //add the particle to the scheduler list
        schedulerList.push_back(temp2.particleID);
        
        std::pair<int,int> temp1ComplexList(-1,-1);
        particleList[0].complexList.push_back(temp1ComplexList);
        particleList[0].complexList.push_back(temp1ComplexList);
        particleList[1].complexList.push_back(temp1ComplexList);
        particleList[2].complexList.push_back(temp1ComplexList);
        
        
       
    }
    else if (readin == 2)
    {
        std::ifstream configInputFile(prm.configReadFile.c_str(), std::ifstream::in);
        if(configInputFile.good())
        {
            
            std::string line;
            configInputFile.clear();
            configInputFile.seekg(0, std::ios::beg);
            
            std::getline(configInputFile, line);
            std::stringstream s1(line);
            particle temp; 
            s1 >>   temp.speciesType        >> 
                    temp.particlePosition.x >>
                    temp.particlePosition.y >>
                    temp.particlePosition.z >>
                    temp.force.x            >>
                    temp.force.y            >>
                    temp.force.z            >>
                    temp.q.q0               >>
                    temp.q.q1               >>
                    temp.q.q2               >>
                    temp.q.q3               >>
                    temp.torque.q0          >>
                    temp.torque.q1          >>
                    temp.torque.q2          >>
                    temp.torque.q3          ;
                   // temp.patchPos.x         >>
                   // temp.patchPos.y         >>
                   // temp.patchPos.z ;
                  
                particleList.push_back(temp);
        }
    }
    
    else if (readin == 3)
    {
        std::ifstream configInputFile(prm.configReadFile.c_str(), std::ifstream::in);
        
        if(configInputFile.good())
        {
            
            std::string line;
            int numberOfLines =0;            
            
            while (getline (configInputFile, line) )
            {
                ++numberOfLines;
            }
           
            configInputFile.clear();
            configInputFile.seekg(0, std::ios::beg);
            
            int N = gsl_ran_flat (rng, 0, numberOfLines-1);    
           
            
            if(N%2 == 0)
            {
                for(int i = 0; i < N; ++i)
                    std::getline(configInputFile, line);
            }
            else
            {
                for(int i = 0; i < N-1; ++i)
                    std::getline(configInputFile, line);
            }
            
            

                std::getline(configInputFile, line);
                std::stringstream s1(line);
                particle temp; 
                
                s1 >>   temp.speciesType        >> 
                        temp.particlePosition.x >>
                        temp.particlePosition.y >>
                        temp.particlePosition.z >>
                        temp.force.x            >>
                        temp.force.y            >>
                        temp.force.z            >>
                        temp.q.q0               >>
                        temp.q.q1               >>
                        temp.q.q2               >>
                        temp.q.q3               >>
                        temp.torque.q0          >>
                        temp.torque.q1          >>
                        temp.torque.q2          >>
                        temp.torque.q3          ;
                        //temp.patchPos.x         >>
                        //temp.patchPos.y         >>
                        //temp.patchPos.z ;

                particleList.push_back(temp);
                
                std::getline(configInputFile, line);
                std::stringstream s2(line);

                particle temp1; 
                s2 >>   temp1.speciesType        >> 
                        temp1.particlePosition.x >>
                        temp1.particlePosition.y >>
                        temp1.particlePosition.z >>
                        temp1.force.x            >>
                        temp1.force.y            >>
                        temp1.force.z            >>
                        temp1.q.q0               >>
                        temp1.q.q1               >>
                        temp1.q.q2               >>
                        temp1.q.q3               >>
                        temp1.torque.q0          >>
                        temp1.torque.q1          >>
                        temp1.torque.q2          >>
                        temp1.torque.q3          ;
                        //temp1.patchPos.x         >>
                        //temp1.patchPos.y         >>
                        //temp1.patchPos.z ;
               
                particleList.push_back(temp1);
            

        }
    }
    
}
