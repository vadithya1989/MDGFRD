//
//  externDefs.h
//  mapK
//
//  Created by Adithya on 30/04/2017.
//  Copyright © 2017 Adithya. All rights reserved.
//

#ifndef externDefs_h
#define externDefs_h

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <vector>

//#include "particle.h"
//#include "species.h"
//#include "simulationParameters.h"


class particle;
class species;
class simulationParameters;

extern double tfactorial[171];


extern gsl_rng * rng;

extern double tBoundC;
extern double tBoundStartC;
extern double tBoundEndC;
extern double numberOfCParticles;
extern double maxCparticles;
extern double probabilityOfBeingBoundC;

extern double tBoundD;
extern double tBoundStartD;
extern double tBoundEndD;
extern double numberOfDParticles;
extern double maxDparticles;
extern double probabilityOfBeingBoundD;


extern double tSim;

extern simulationParameters prm;/**Object of simulationParameters class that contains all the parameters*/

extern std::vector<particle> particleList;/**list of particles*/
extern std::vector<species> speciesList;/**list of various species used in the simulation */
extern std::vector<int> schedulerList;/**list of particles from the particle list arranged in the order of
                                       their next event times*/
extern std::vector<int> idleParticleList;/**list of all kinases/phosphotases in a complex*/

extern std::vector<int> cells;/**cell list for the neighbor build*/
extern std::vector<int> particles;/**particleS list for the neighbor build*/

extern double dissociationRate;
extern double blockedDissociationRate;
extern double singleParticleDissociationRate;
extern double hoppingRate;
extern double phosphorylationRate;

#endif /* externDefs_h */
