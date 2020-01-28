
#include "particle.h"
#include "simulationParameters.h"

#include <cmath>
#include <iostream>
#include <iomanip>

void updatePosition(
                    std::vector<myVector> &randomNoise,
                    std::vector<myVector> &randomAngularNoise
                    )
{
    std::vector<int> idleEnzymeList;
    for (int i=0; i<particleList.size(); ++i)
    {
        
        if(particleList[i].particleType == 2)
        {
            int pass=0;
            for(std::vector<int>::iterator it=idleParticleList.begin();it!=idleParticleList.end();++it)
            {
                if(*it == i)
                {
                    pass=1;
                }
            }
            
            if(pass == 0)
            {
                
                //UPDATE POSITION
                particleList[i].particlePosition = particleList[i].particlePosition +  particleList[i].force *  prm.deltaT/( (prm.gamma)*(speciesList[particleList[i].speciesType].particleMass) ) + randomNoise[i] ;
                
                //UPDATE ORIENTATION AND ANGULAR MOMENTUM
                long double M = (4.0/3.0)*(2.0/5.0)*speciesList[particleList[i].speciesType].particleMass * speciesList[particleList[i].speciesType].particleRadius * speciesList[particleList[i].speciesType].particleRadius;
                
                long double Y01 = (particleList[i].torque.q0*particleList[i].q.q1 - particleList[i].torque.q1*particleList[i].q.q0)*prm.deltaT/(prm.gammaRot*M) - randomAngularNoise[i].x;
                long double Y02 = (particleList[i].torque.q0*particleList[i].q.q2 - particleList[i].torque.q2*particleList[i].q.q0)*prm.deltaT/(prm.gammaRot*M) - randomAngularNoise[i].y;
                long double Y03 = (particleList[i].torque.q0*particleList[i].q.q3 - particleList[i].torque.q3*particleList[i].q.q0)*prm.deltaT/(prm.gammaRot*M) - randomAngularNoise[i].z;
                long double Y12 = (particleList[i].torque.q1*particleList[i].q.q2 - particleList[i].torque.q2*particleList[i].q.q1)*prm.deltaT/(prm.gammaRot*M) + randomAngularNoise[i].z;
                long double Y13 = (particleList[i].torque.q1*particleList[i].q.q3 - particleList[i].torque.q3*particleList[i].q.q1)*prm.deltaT/(prm.gammaRot*M) - randomAngularNoise[i].y;
                long double Y23 = (particleList[i].torque.q2*particleList[i].q.q3 - particleList[i].torque.q3*particleList[i].q.q2)*prm.deltaT/(prm.gammaRot*M) + randomAngularNoise[i].x;
                
                //CALCULATE THE EXPONENT OF THE Y-MATRIX
                //calculate a, b and mu
                long double a0 = (Y23*Y01 + Y03*Y12 - Y13*Y02) * (Y23*Y01 + Y03*Y12 - Y13*Y02);
                long double a2 = (Y23*Y23) + (Y01*Y01) + (Y03*Y03) +
                (Y12*Y12) + (Y13*Y13) + (Y02*Y02);
                long double del = sqrt( (a2*a2)- (double(4.0)*a0) );
                if((a2*a2)- (double(4.0)*a0)>=-1e-8 && (a2*a2)- (double(4.0)*a0) < 0)
                del = double(0);
                long double alpha = sqrt ( double(0.5)* (a2-del) );
                if( (double(0.5)*(a2-del))>=-1e-8 && (double(0.5)*(a2-del)) < 0)
                alpha = double(0);
                long double mu = sqrt ( double(0.5)* (a2+del) );
                long double a,b;
                
                if(del < 1e-8)
                {
                    long double gamma = sqrt(a2*double(0.5));
                    a = ( sin(gamma)/gamma - cos(gamma) )/a2;
                    b = sin(gamma)/(double(2)*gamma);
                }
                else if( a2 < 1e-8)
                {
                    a = double(1)/ double(6) -  double(1)/ double(120)*a2;
                    b = double(1)/ double(2) -  double(1)/ double(24)*a2;
                }
                else if( alpha < 1e-4)
                {
                    a = ( ((alpha*alpha)/20.0-1.0)*(alpha*alpha/6.0) + 1.0 - (sin(mu)/mu) ) / del;
                    b = ( cos(alpha)-cos(mu) )/del;
                }
                else
                {
                    a = ( (sin(alpha)/alpha) - (sin(mu)/mu) ) / del;
                    b = ( cos(alpha)-cos(mu) )/del;
                }
                //write expression for Exp(Ymatrix)
                long double expYmatrix[4][4];
                expYmatrix[0][0] = b* (mu*mu - Y01*Y01 - Y02*Y02 - Y03*Y03) + a* Y01*(-Y02*Y12 - Y03*Y13) + a* Y03*(Y01*Y13 + Y02*Y23) + a* Y02*(Y01*Y12 - Y03*Y23)   + cos(mu);
                expYmatrix[0][1] = b* (-Y02*Y12 - Y03*Y13) + a* Y01*(mu*mu - Y01*Y01 - Y12*Y12 - Y13*Y13) + a* Y03*(-Y01*Y03 + Y12*Y23) + a* Y02*(-Y01*Y02 - Y13*Y23) + (Y01*sin(mu)/mu);
                expYmatrix[0][2] = a* Y03*(-Y02*Y03 - Y12*Y13) + b* (Y01*Y12 - Y03*Y23) +  a* Y01*(-Y01*Y02 - Y13*Y23) + a* Y02*(mu*mu - Y02*Y02 - Y12*Y12 - Y23*Y23) + (Y02*sin(mu)/mu);
                expYmatrix[0][3] = a* Y02*(-Y02*Y03 - Y12*Y13) + b* (Y01*Y13 + Y02*Y23) + a* Y01*(-Y01*Y03 + Y12*Y23) + a* Y03*(mu*mu - Y03*Y03 - Y13*Y13 - Y23*Y23)  + (Y03*sin(mu)/mu);
                expYmatrix[1][0] = -a* Y01*(mu*mu - Y01*Y01 - Y02*Y02 - Y03*Y03) + b* (-Y02*Y12 - Y03*Y13) + a* Y13*(Y01*Y13 + Y02*Y23) + a* Y12*(Y01*Y12 - Y03*Y23)  - (Y01*sin(mu)/mu);;
                expYmatrix[1][1] = -a*Y01*(-Y02*Y12 - Y03*Y13) + b* (mu*mu - Y01*Y01 - Y12*Y12 - Y13*Y13) + a*Y13*(-Y01*Y03 + Y12*Y23) + a*Y12*(-Y01*Y02 - Y13*Y23)   + cos(mu);
                expYmatrix[1][2] = a*Y13*(-Y02*Y03 - Y12*Y13) - a*Y01*(Y01*Y12 - Y03*Y23) + b* (-Y01*Y02 - Y13*Y23) + a*Y12*(mu*mu - Y02*Y02 - Y12*Y12 - Y23*Y23)     + (Y12*sin(mu)/mu);;
                expYmatrix[1][3] = a*Y12*(-Y02*Y03 - Y12*Y13) - a*Y01*(Y01*Y13 + Y02*Y23) + b* (-Y01*Y03 + Y12*Y23) + a*Y13*(mu*mu - Y03*Y03 - Y13*Y13 - Y23*Y23)     + (Y13*sin(mu)/mu);;
                expYmatrix[2][0] = -a*Y02*(mu*mu - Y01*Y01 - Y02*Y02 - Y03*Y03) - a*Y12*(-Y02*Y12 - Y03*Y13) +  a*Y23*(Y01*Y13 + Y02*Y23) + b* (Y01*Y12 - Y03*Y23)    - (Y02*sin(mu)/mu);;
                expYmatrix[2][1] = -a*Y02*(-Y02*Y12 - Y03*Y13) - a*Y12*(mu*mu - Y01*Y01 - Y12*Y12 - Y13*Y13) + a*Y23*(-Y01*Y03 + Y12*Y23) + b* (-Y01*Y02 - Y13*Y23)   - (Y12*sin(mu)/mu);;
                expYmatrix[2][2] = a*(-Y02*Y03 - Y12*Y13)*Y23 - a*Y02*(Y01*Y12 - Y03*Y23) - a*Y12*(-Y01*Y02 - Y13*Y23) + b* (mu*mu - Y02*Y02 - Y12*Y12 - Y23*Y23)     +  cos(mu);
                expYmatrix[2][3] = b* (-Y02*Y03 - Y12*Y13) - a*Y02*(Y01*Y13 + Y02*Y23) - a*Y12*(-Y01*Y03 + Y12*Y23) + a*Y23*(mu*mu - Y03*Y03 - Y13*Y13 - Y23*Y23)     + (Y23*sin(mu)/mu);;
                expYmatrix[3][0] = -a*Y03*(mu*mu - Y01*Y01 - Y02*Y02 - Y03*Y03) - a*Y13*(-Y02*Y12 - Y03*Y13) + b* (Y01*Y13 + Y02*Y23) - a*Y23*(Y01*Y12 - Y03*Y23)     - (Y03*sin(mu)/mu);;
                expYmatrix[3][1] = -a*Y03*(-Y02*Y12 - Y03*Y13) - a*Y13*(mu*mu - Y01*Y01 - Y12*Y12 - Y13*Y13) +  b* (-Y01*Y03 + Y12*Y23) - a*Y23*(-Y01*Y02 - Y13*Y23)  - (Y13*sin(mu)/mu);;
                expYmatrix[3][2] = b* (-Y02*Y03 - Y12*Y13) - a*Y03*(Y01*Y12 - Y03*Y23) - a*Y13*(-Y01*Y02 - Y13*Y23) - a*Y23*(mu*mu - Y02*Y02 - Y12*Y12 - Y23*Y23)     - (Y23*sin(mu)/mu);;
                expYmatrix[3][3] =-a*(-Y02*Y03 - Y12*Y13)*Y23 -  a*Y03*(Y01*Y13 + Y02*Y23) - a*Y13*(-Y01*Y03 + Y12*Y23) + b* (mu*mu - Y03*Y03 - Y13*Y13 - Y23*Y23)    + cos(mu);
                
                /*std::cout << expYmatrix[0][0] << " " << expYmatrix[0][1] << " " << expYmatrix[0][2] << " " << expYmatrix[0][3] << std::endl << " " << expYmatrix[1][0] << " "
                 << expYmatrix[1][1] << " " << expYmatrix[1][2] << " " << expYmatrix[1][3] << std::endl<< " " << expYmatrix[2][0] << " " << expYmatrix[2][1] << " "
                 << expYmatrix[2][2] << " " << expYmatrix[2][3] << std::endl<< " " << expYmatrix[3][0] << " " << expYmatrix[3][1] << " " << expYmatrix[3][2] << " "
                 << expYmatrix[3][3] << " " << std::endl << std::endl;*/
                
                quaternion temp = particleList[i].q;
                
                particleList[i].q.q0 = expYmatrix[0][0]*temp.q0 + expYmatrix[0][1]*temp.q1 + expYmatrix[0][2]*temp.q2 + expYmatrix[0][3]*temp.q3;
                particleList[i].q.q1 = expYmatrix[1][0]*temp.q0 + expYmatrix[1][1]*temp.q1 + expYmatrix[1][2]*temp.q2 + expYmatrix[1][3]*temp.q3;
                particleList[i].q.q2 = expYmatrix[2][0]*temp.q0 + expYmatrix[2][1]*temp.q1 + expYmatrix[2][2]*temp.q2 + expYmatrix[2][3]*temp.q3;
                particleList[i].q.q3 = expYmatrix[3][0]*temp.q0 + expYmatrix[3][1]*temp.q1 + expYmatrix[3][2]*temp.q2 + expYmatrix[3][3]*temp.q3;
                
                
            }
            else if(pass == 1)
            {
                
                idleEnzymeList.push_back(i);
            }
        }
    }
    
    //now update the positions of the particles in the idleEnzymeList based on their bound substrate
   
    if(idleEnzymeList.size()>0)
    {
        for(std::vector<int>::iterator it=idleEnzymeList.begin();it!=idleEnzymeList.end();++it)
        {
            int correspondingSubstrate = particleList[*it].complexList[0].first;
            myVector patchDirection =speciesList[particleList[correspondingSubstrate].speciesType].bodyPatchVectors[particleList[*it].complexList[0].second];
            
            particleList[*it].particlePosition = particleList[correspondingSubstrate].q.bodyToSpace(patchDirection) * 2.0*speciesList[particleList[correspondingSubstrate].speciesType].particleRadius + particleList[correspondingSubstrate].particlePosition;
        }
    }
    //UPDATE PARTICLE PATCH POSITION
    for(int p=0; p<particleList.size();++p)
    {
        
        for(int i=0;i<speciesList[particleList[p].speciesType].numberOfPatches;++i)
        {
            
            myVector patchDirection = speciesList[particleList[p].speciesType].bodyPatchVectors[i];
            myVector tempPatchPos = particleList[p].q.bodyToSpace(patchDirection) * speciesList[particleList[p].speciesType].particleRadius + particleList[p].particlePosition;
            particleList[p].patchPos[i] = tempPatchPos;
            
        }
        
    }
    
}
