
//#include <gsl/gsl_rng.h>
//#include <gsl/gsl_randist.h>
//
//#include <cmath>
//#include <iostream>
//
//#include "simulationParameters.h"
//#include "particle.h"
#include "main.h"
#include "externDefs.h"

/**
 Updates the position and the orientation of the langevin particle
 @param randomNoise randomNoise
 @param randomAngularNoise randomAngularNoise
 */
void updatePosition (
                     std::vector<myVector> &randomNoise,
                     std::vector<myVector> & randomAngularNoise);


/**
 Calculates the derivative of the potential in order to calculate the force
 */
void calcForce  ();


/**
 This removes the component of torque that is responsible for the elongation of the
 torque and makes it length differ from unity
 */
void removeTorqueComponent ();


void integrator ()
{
    std::vector<myVector> randomNoise;/**a vector that holds the random noise part for each particle*/
    std::vector<myVector> randomAngularNoise;/**vector that holds the random noise part for the orientations for each particle*/
    
    for ( int i=0; i<particleList.size(); ++i)
    {
        long double Mass = (4.0/3.0)*(2.0/5.0)*speciesList[particleList[i].speciesType].particleMass * speciesList[particleList[i].speciesType].particleRadius * speciesList[particleList[i].speciesType].particleRadius;
        long double standardDeviation = sqrt ( prm.deltaT * double(2)*prm.boltzmannConstant*prm.temperature/(prm.gamma*speciesList[particleList[i].speciesType].particleMass) );/**standard deviation of the random noise*/
        long double standardDeviationRotation = sqrt (prm.deltaT * double(2)*prm.boltzmannConstant*prm.temperature/( prm.gammaRot*Mass) );/**standard deviation of the random angular noise*/
        

        
        myVector temp;
        temp.x =  gsl_ran_gaussian ( rng, standardDeviation );
        temp.y =  gsl_ran_gaussian ( rng, standardDeviation );
        temp.z =  gsl_ran_gaussian ( rng, standardDeviation );
        randomNoise.push_back(temp);
        
        myVector tempRot;
        tempRot.x = gsl_ran_gaussian ( rng, standardDeviationRotation );
        tempRot.y = gsl_ran_gaussian ( rng, standardDeviationRotation );
        tempRot.z = gsl_ran_gaussian ( rng, standardDeviationRotation );
        randomAngularNoise.push_back(tempRot);
        
        
    }
    
    
    calcForce       ( );
    //std::cout << "calc Force" << std::endl;
    removeTorqueComponent( );
    //std::cout << "rem tor comp" << std::endl;
    updatePosition  ( randomNoise, randomAngularNoise);
    //std::cout << "update position" << std::endl;
}
