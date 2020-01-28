
#include <vector>
#include <algorithm>
#include <cmath>
#include <iostream>
#include "simulationParameters.h"
#include "particle.h"
#include "externDefs.h"
#include "main.h"

std::pair <double,int> returnNearestDistance(
                                             int &cellID,
                                             int &particleId,
                                             std::string condition)
{
    
    double dist = double(0);/**stores the nearest particle distance*/
    std::vector<double> distances;/**list of all distances the minimum of which is the nearest particle distance*/
    int nnid = 0;/**the id of the nearest neighbor*/
    std::vector<int> nnidx;/**list of all ids of the neighbor, the id of the corresponding particle to the min distance is the nearest particle id*/
    
    int cIndex = cellID;
    int i = cells[cIndex];
    
    int cIndexX = cIndex % prm.xcells;
    cIndex /= prm.xcells;
    
    int cIndexY = cIndex % prm.ycells;
    cIndex /= prm.ycells;
    
    int cIndexZ = cIndex;
    
    int cIX = cIndexX;
    int cIY = cIndexY;
    int cIZ = cIndexZ;
    
    while ( i != -1 )
    {
        for (int dx = -1; dx <= 1; ++dx)
        {
            for (int dy = -1; dy <= 1; ++dy)
            {
                for (int dz = -1; dz <= 1; ++dz)
                {
                    cIndexX = cIX + dx;
                    cIndexY = cIY + dy;
                    cIndexZ = cIZ + dz;
                    
                    /*comment out the following three lines and uncomment the next and comment the next to not apply periodic boundary conditions*/
                    /*APPLYING PERIODIC BOUNDARY CONDITIONS*/
                    int cIndexXPeriodic = cIndexX - floor(double(cIndexX)/double(prm.xcells))*prm.xcells;
                    int cIndexYPeriodic = cIndexY - floor(double(cIndexY)/double(prm.ycells))*prm.ycells;
                    int cIndexZPeriodic = cIndexZ - floor(double(cIndexZ)/double(prm.zcells))*prm.zcells;
                    /*int cIndexCur = (cIndexZ * prm.xcells * prm.ycells) + (cIndexY * prm.xcells) + cIndexX;*/
                    int cIndexCur = (cIndexZPeriodic * prm.xcells * prm.ycells) + (cIndexYPeriodic * prm.xcells) + cIndexXPeriodic;/*mapping from co-ordinate to linear*/
                    
                    if (cIndexCur >= 0 && cIndexCur < (prm.xcells*prm.ycells*prm.zcells ) )
                    {
                        //std::cout << (particleList[0].particlePosition - particleList[1].particlePosition).magnitude() << std::endl;
                        int j = cells[cIndexCur];
                        while (j != -1)
                        {
                            if(i == particleId)
                            {
                                
                                /*vector r gives the distance between two domain centers*/
                                myVector rPb = particleList[i].particlePosition - particleList[j].particlePosition;
                                /*Applying periodic boundary condition*/
                                double delx = rPb.x;
                                double dely = rPb.y;
                                double delz = rPb.z;
                                myVector r;
                                r.x = rPb.x - ((floor((delx/prm.xmax)+0.5)) * prm.xmax );
                                r.y = rPb.y - ((floor((dely/prm.ymax)+0.5)) * prm.ymax );
                                r.z = rPb.z - ((floor((delz/prm.zmax)+0.5)) * prm.zmax );
                                if(condition.compare("particle") == 0)
                                {
                                    int pass = 1;
                                    for(std::vector<int>::iterator it=idleParticleList.begin();it!=idleParticleList.end();++it)
                                    {
                                        if(particleList[i].particleID == *it || particleList[j].particleID == *it)
                                        {
                                            pass =0;
                                        }
                                    }
                                    if(pass == 1)
                                    {
                                        if ((r.magnitude() > 0.0) && (particleList[i].particleType == 2) && (particleList[j].particleType == 2) ) /*measure distance only if both particles are LD type particles*/
                                        {
                                            distances.push_back(r.magnitude());
                                            nnidx.push_back(j);
                                        }
                                    }
                                    
                                }
                                if(condition.compare("domain") == 0)
                                {
                                    int pass =1 ;
                                    for(std::vector<int>::iterator it=idleParticleList.begin();it!=idleParticleList.end();++it)
                                    {
                                        if(particleList[i].particleID == *it || particleList[j].particleID == *it)
                                        {
                                            pass = 0;
                                        }
                                    }
                                    if(pass == 1)
                                    {
                                        if ((r.magnitude() > 0.0) && (particleList[i].particleType == 2 ) && (particleList[j].particleType == 1) ) /*measure distance only if given particle is a LD particle or zero single and the neighbor is a single domain*/
                                        {
                                            distances.push_back(r.magnitude() - particleList[j].domainRadius);
                                            nnidx.push_back(j);
                                        }
                                    }
                                    
                                }
                                if(condition.compare("neighbor") == 0)
                                {
                                    
                                    int pass=1;
                                    for(std::vector<int>::iterator it=idleParticleList.begin();it!=idleParticleList.end();++it)
                                    {
                                        
                                        if(particleList[i].particleID == *it || particleList[j].particleID == *it)
                                        {
                                            pass =0;
                                        }
                                    }
                                    if(pass == 1)
                                    {
                                        if (r.magnitude() > 0.0) /*measure distance only of the particles are different (i!=j)*/
                                        {
                                            if(particleList[j].particleType == 2 && particleList[i].particleType == 2 )/*if the neighboring particle is a zero single or a LD particle then */
                                            {
                                                //if the distance to a neighbouring particle is greater than twice the min size only then a domain can be built
                                                if( (r.magnitude() ) >= 2*prm.minDomainRadius+prm.rc)
                                                {
                                                    distances.push_back( (r.magnitude()-prm.rc) * double(0.5));
                                                    nnidx.push_back(j);
                                                }
                                                else/*if there is not enough space dist = particle radius*/
                                                {
                                                    distances.push_back(-1);
                                                    nnidx.push_back(-1);
                                                }
                                            }
                                            else if (particleList[j].particleType == 1 && particleList[i].particleType == 2 )/*if the neighbor is a single domain*/
                                            {
                                                /*if the distance available is greater than the min domain radius*/
                                                if( (r.magnitude() - particleList[j].domainRadius) >= prm.minDomainRadius+prm.rc)
                                                {
                                                    
                                                    /*then distance is till the domain surface*/
                                                    distances.push_back(r.magnitude() - particleList[j].domainRadius - prm.rc);
                                                    nnidx.push_back(j);
                                                }
                                                else
                                                {
                                                    /*if not then again make it a LD particle*/
                                                    distances.push_back(-1);
                                                    nnidx.push_back(-1);
                                                }
                                            }
                                        }
                                        
                                    }
                                    
                                }
                            }
                            j = particles[j];
                        }
                    }
                }
            }
        }
        i = particles[i];
    }
    
    if(!distances.empty() && distances.size() > 1)
    {
        //the least value in the distances vector gives the distance to the nearest neighbor
        dist = *std::min_element (distances.begin(), distances.end());
        nnid = nnidx[std::distance ( distances.begin(), std::min_element (distances.begin(), distances.end()) )];
    }
    else if (distances.size() == 1)
    {
        dist = distances[0];
        nnid = nnidx[0];
    }
    else
    {
        dist = -1.0;
        nnid = -1;
    }
    std::pair<double, int> temp;
    temp = std::make_pair (dist, nnid);
    return temp;
}

