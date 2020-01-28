
#include "simulationParameters.h"
#include "particle.h"
#include <cmath>
#include <iostream>

double gradPotential(double r,
                    double a,
                    double x0,
                    double dxstar,
                    double b,
                    double dxc,
                    const simulationParameters &prm)
{
        
    if( (x0-dxstar)<r && (x0+dxstar)>r )
    {
        return (-2.0*(a/prm.sigma)*(r/prm.sigma-x0/prm.sigma));
    }
    else if( (x0-dxstar)>r && (x0-dxc)<r )
    {
        return (2.0*(b/prm.sigma)*(r/prm.sigma-(x0/prm.sigma-dxc/prm.sigma)) );
    }
    else if( (x0+dxstar)<r && (x0+dxc)>r )
    {
        return (2.0*(b/prm.sigma)*(r/prm.sigma-(x0/prm.sigma+dxc/prm.sigma)) );
    }
    else
    {
        return 0.0;
    }

}

void drPatch_dq ( particle &part,
                        particle &partNei,
                        std::vector< std::pair<int, myVector> > &rPatchVector,
                        std::vector<species> &speciesList,  
                        std::vector<quaternion> &drPatchdq
                      )
{
    
    for (int m=0; m<rPatchVector.size(); ++m)
    {
        ///independent of the type of potential used
        ///l denotes the nth patch on the particle in question
        ///bodyPatchVectors[l] -> bodypatchVectors[m] where m =
        int l = int(m/speciesList[partNei.speciesType].numberOfPatches);
        double x = speciesList[part.speciesType].bodyPatchVectors[l].x;
        double y = speciesList[part.speciesType].bodyPatchVectors[l].y;
        double z = speciesList[part.speciesType].bodyPatchVectors[l].z;
        double q0 = part.q.q0;
        double q1 = part.q.q1;
        double q2 = part.q.q2;
        double q3 = part.q.q3;
        myVector delVec = rPatchVector[m].second;
        double del = rPatchVector[m].second.magnitude();
        quaternion temp;
        
        temp.q0 =  (-2.0*( q0*x-q3*y+q2*z)*delVec.x/del
                    -2.0*( q3*x+q0*y-q1*z)*delVec.y/del
                    -2.0*(-q2*x+q1*y+q0*z)*delVec.z/del)*speciesList[part.speciesType].particleRadius;
        
        temp.q1 =  (-2.0*( q1*x+q2*y+q3*z)*delVec.x/del
                    -2.0*( q2*x-q1*y-q0*z)*delVec.y/del
                    -2.0*( q3*x+q0*y-q1*z)*delVec.z/del)*speciesList[part.speciesType].particleRadius;
        
        temp.q2 =  (-2.0*(-q2*x+q1*y+q0*z)*delVec.x/del
                    -2.0*( q1*x+q2*y+q3*z)*delVec.y/del
                    -2.0*(-q0*x+q3*y-q2*z)*delVec.z/del)*speciesList[part.speciesType].particleRadius;
        
        temp.q3 =  (-2.0*(-q3*x-q0*y+q1*z)*delVec.x/del
                    -2.0*( q0*x-q3*y+q2*z)*delVec.y/del
                    -2.0*( q1*x+q2*y+q3*z)*delVec.z/del)*speciesList[part.speciesType].particleRadius;
        drPatchdq.push_back(temp);
    }
    
}

void dVatt_drPatch (    int &substrateID,
                        int &enzymeID,
                        std::vector< std::pair<int, myVector> > &rPatchVector,
                        std::vector<double> &dVattdrPatch,
                        myVector &rPb)
{
    ///for a 6-3 LJ potential of V=-epsilon*(sigmat/rPatch)^3 + V0 + V1 changes for different potentials
    for (int i=0; i<rPatchVector.size(); ++i)
    {
        double temp;
        double scale=0;
        if( (particleList[substrateID].patchState[rPatchVector[i].first]==0 && particleList[enzymeID].speciesType==3)
        ||  (particleList[substrateID].patchState[rPatchVector[i].first]==1 && particleList[enzymeID].speciesType==4))
        {
            
            scale = 20.0*prm.boltzmannConstant*prm.temperature;
        }
        
        double batt = 5.0;
        double dxcatt = 0.5*prm.sigma;        

        double aatt = 20;
        double x0att = 0;
        double dxstaratt = 0.1*prm.sigma;
    
        temp = -scale*gradPotential(rPatchVector[i].second.magnitude(), aatt, x0att, dxstaratt, batt, dxcatt, prm);
        
        dVattdrPatch.push_back(temp);
    }
}

void drPatch_dCoordinate(std::vector< std::pair<int, myVector> > &rPatchVector, std::vector<myVector> &drPatchdCoordinate)
{
    ///Independent of the type of potential used
    for (int i=0; i<rPatchVector.size(); ++i)
    {
        myVector temp;
        temp.x = -rPatchVector[i].second.x/rPatchVector[i].second.magnitude();
        temp.y = -rPatchVector[i].second.y/rPatchVector[i].second.magnitude();
        temp.z = -rPatchVector[i].second.z/rPatchVector[i].second.magnitude();
        drPatchdCoordinate.push_back(temp);
    }
}

myVector dVrep_dCoordiate( myVector &rPb)
{
    ///repulsion plus isotropic attraction
    double scale = 20.0*prm.boltzmannConstant*prm.temperature;
    double rfactor = 5;
    double brep = 2.6036;
    double dxcrep = 1.17647*prm.sigma;

    double bisoAtt = 2.6036;
    double dxcisoAtt = 1.17647*prm.sigma;
    double aisoAtt = 1.0;
    double x0isoAtt = 0.4*prm.sigma; 
    double dxstarisoAtt = 0.85*prm.sigma;
    
    double arep=1;
    double x0rep=0;
    double dxstarrep = 0.85*prm.sigma;

    ///dVrepdCoordiate = dVrep/dR*dR/dCoordinate
    myVector dVrepdCoordiate;
    
    double dVrepdR = scale*rfactor*gradPotential(rPb.magnitude(), arep, x0rep, dxstarrep, brep, dxcrep, prm)
                        -scale*0.5*gradPotential(rPb.magnitude(), aisoAtt, x0isoAtt, dxstarisoAtt, bisoAtt, dxcisoAtt, prm);
    myVector dRdCoordiate = rPb*-double(1)/rPb.magnitude();  
   
    
    dVrepdCoordiate = dRdCoordiate * dVrepdR;
    
    return dVrepdCoordiate;
    
}

myVector dVrepWCA_dr (const simulationParameters& prm, myVector &rPb)
{
    ///repulsion
    double scale = 20.0*prm.boltzmannConstant*prm.temperature;
    double rfactor = 5;
    double brep = 2.6036;
    double dxcrep = 1.17647*prm.sigma;
    
    double arep=1;
    double x0rep=0;
    double dxstarrep = 0.85*prm.sigma;

    ///dVrepdCoordiate = dVrep/dR*dR/dCoordinate
    myVector dVrepdCoordiate;
    
    double dVrepdR = scale*rfactor*gradPotential(rPb.magnitude(), arep, x0rep, dxstarrep, brep, dxcrep, prm);
                        
    myVector dRdCoordiate = rPb*-double(1)/rPb.magnitude();  
   
    
    dVrepdCoordiate = dRdCoordiate * dVrepdR;
    
    return dVrepdCoordiate;
}
