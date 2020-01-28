
//#include "particle.h"
//#include "simulationParameters.h"
#include "main.h"

//#include <cmath>
//#include <iostream>

/**
 This calculates the derivative of the patch vector wrt the quaternions. This
 quaternion depends only on the orientation and disance of the particles and not
 the potential
 
 @param pi particle i
 @param pj particle j
 @param rPatchVector inter patch vector
 @param speciesList list of species
 @param drPatchDq gradient of the patch vector w.r.t the quaternions
 */
void drPatch_dq (particle &pi,
                       particle &pj,
                       std::vector< std::pair<int, myVector> > &rPatchVector,
                       std::vector<species> &speciesList,
                       std::vector<quaternion> &drPatchDq );

/**
 Calculates the derivative of the potential wrt to the dist bet patches.
 This value changes with the potential
 @param substrateID substrate particle ID
 @param rPatchVector inter patch vector
 @param dvAttDrPatch derivative of the potential wrt to the distance bet patches
 @param r inter particle vector
 */
void dVatt_drPatch (int &substrateID,
                    int &enzymeID,
                    std::vector< std::pair<int, myVector> >  &rPatchVector,
                    std::vector<double> &dvAttDrPatch,
                    myVector &r);


/**
 This function calculates the derivative of the inter particle patch distance
 wrt to the co-ordinate of the given particle. Independent of the potential used

 @param rPatchVector all inter patch vectors
 @param drPatchDCoordinate derivative of the inter particle patch distance
 wrt to the co-ordinate
 */
void drPatch_dCoordinate(std::vector< std::pair<int, myVector> > &rPatchVector,
                         std::vector<myVector> &drPatchDCoordinate);

/**
 This function calculates the derivative of the repulsive potential wrt to coordiate of the current particle
 Depends on the potential used
 @param r inter particle vector
 @return dVrepdCoord which is the derivative of the repulsive potential wrt to coordiate of the current particle
 */
myVector dVrep_dCoordiate(myVector &r);


/**
 Calculates the derivative of the WCA repulsive potential wrt r for like particles to repel each other

 @param prm simualtion parameters
 @param r inter particle distance vector
 @return derivative of the WCA repulsive potential wrt r
 */
myVector dVrepWCA_dr (const simulationParameters &prm,
                      myVector &r);

void calcForce ()
{
    for ( int i=0; i<particleList.size(); ++i)
    {
        /*RECALCULATING FORCES*/
        particleList[i].force  = double(0);
        particleList[i].torque = double(0);
        
        for( int j=0; j<particleList.size(); ++j)
        {
            if ( i!=j )
            {
                /*calculate the distance between the particles*/
                myVector r;               
                measureInterParticleDistace(r, particleList[i].particlePosition, particleList[j].particlePosition);
                
                /*calculate the patch vector array*/
                std::vector< std::pair<int, myVector> > rPatchVector;
                measureInterPatchDistance(r, rPatchVector, particleList[i], particleList[j] );
                
                
                /*FORCE CALCULATION*/
                /*Force = -(dV/dCoord) = -(dV/drPatch*drPatch/dCoord)*/
                myVector dVtotdCoord ;/**Gives the gradient of the total potential wrt the coordiate of the particle*/
                myVector dVattdCoord ;
                myVector dVrepdCoord ;/**Gives the gradient of the repulsive potential wrt the coordiate of the particle*/
                std::vector<myVector> drPatchDCoordinate;
                std::vector<double> dVattDrPatch;
                myVector Force;               
               ///the particles have arrtaction only if one particle is a substrate and the other is an enzyme()
                if(( (particleList[i].speciesType == 0 || particleList[i].speciesType == 1 || particleList[i].speciesType == 2) && (particleList[j].speciesType == 3 || particleList[j].speciesType == 4)
                   ) ||
                   ( (particleList[j].speciesType == 0 || particleList[j].speciesType == 1 || particleList[j].speciesType == 2) && (particleList[i].speciesType == 3 || particleList[i].speciesType == 4)
                   )
                  )
                {
                    int substrate;
                    int enzyme;
                    if(( (particleList[i].speciesType == 0 || particleList[i].speciesType == 1 || particleList[i].speciesType == 2) && (particleList[j].speciesType == 3 || particleList[j].speciesType == 4)
                        )
                      )
                    {
                        substrate = i;
                        enzyme = j;
                    }
                    else
                    {
                        substrate = j;
                        enzyme = i;
                    }
                    dVrepdCoord = dVrep_dCoordiate(r);
                    drPatch_dCoordinate(rPatchVector, drPatchDCoordinate);
                    dVatt_drPatch(substrate, enzyme, rPatchVector, dVattDrPatch, r);

                    for(int p=0;p<rPatchVector.size();++p)
                    {                    
                        dVattdCoord.x += drPatchDCoordinate[p].x*dVattDrPatch[p];
                        dVattdCoord.y += drPatchDCoordinate[p].y*dVattDrPatch[p];
                        dVattdCoord.z += drPatchDCoordinate[p].z*dVattDrPatch[p];
                    }

                    dVtotdCoord = dVrepdCoord + dVattdCoord;                
                    Force = dVtotdCoord*-1;

                    particleList[i].force = particleList[i].force + Force;
                    
                    /*Torque calculation*/
                    /*Torque = -dV/dq = - (dV/drPatchVector*drPatchVector/dq)*/

                    quaternion dVtotdq;
                    quaternion dVrepdq;
                    quaternion dVattdq;
                    quaternion Torque;   
                    std::vector<quaternion> drPatchDq;
                    dVrepdq = double(0);

                    drPatch_dq(particleList[i], particleList[j], rPatchVector, speciesList, drPatchDq) ;

                    for(int p=0;p<rPatchVector.size();++p)
                    {                    
                        dVattdq.q0 += drPatchDq[p].q0*dVattDrPatch[p];
                        dVattdq.q1 += drPatchDq[p].q1*dVattDrPatch[p];
                        dVattdq.q2 += drPatchDq[p].q2*dVattDrPatch[p];
                        dVattdq.q3 += drPatchDq[p].q3*dVattDrPatch[p];

                        quaternion temp;
                        temp.q0 = drPatchDq[p].q0;//*dVattDrPatch[p];
                        temp.q1 = drPatchDq[p].q1;//*dVattDrPatch[p];
                        temp.q2 = drPatchDq[p].q2;//*dVattDrPatch[p];
                        temp.q3 = drPatchDq[p].q3;//*dVattDrPatch[p];                    
                    }




                    dVtotdq = dVrepdq + dVattdq;
                    Torque = dVtotdq*-1.0;              


                    particleList[i].torque = particleList[i].torque + Torque;  

                    //particleList[i].torque = particleList[i].torque + dUdq*-1;  
                }
                else
                {
                    //calculate the derivative of the repulsive potential wrt r
                    myVector dVrepdr = dVrepWCA_dr( prm, r);
                    Force = dVrepdr*-1;
                    /*the gradient of the potential is obtained by chain rule. dvdx = dvdr*drdx*/
                    particleList[i].force = particleList[i].force + Force;
                }
            }
        }
    }
        
}
