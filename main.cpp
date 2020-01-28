/* 
 * File:   main.cpp
 * Author: adithya
 *
 * Created on December 13, 2014, 11:21 AM
 */

#include "main.h"


/*GLOBAL RANDOM NUMBER GENERATOR*/
 gsl_rng * rng;


double tfactorial[171]; 

/*GLOBAL COUBTER VARIABLES*/
double tBoundC;
double tBoundStartC;
double tBoundEndC;
double numberOfCParticles;
double maxCparticles=1;
double probabilityOfBeingBoundC;

double tBoundD;
double tBoundStartD;
double tBoundEndD;
double numberOfDParticles;
double maxDparticles=1;
double probabilityOfBeingBoundD;
double tSim = 0.0;

/**LISTS*/
std::vector<particle> particleList;/**list of particles*/
std::vector<species> speciesList;/**list of various species used in the simulation */
std::vector<int> schedulerList;/**list of particles from the particle list arranged in the order of
                                their next event times*/
std::vector<int> idleParticleList;/**list of all kinases/phosphotases in a complex*/
simulationParameters prm;/**Object of simulationParameters class that contains all the parameters*/
/**CELL AND PARTICLE LIST FOR THE NEIGHBOR LIST*/
std::vector<int> cells;/**cell list for the neighbor build*/
std::vector<int> particles;/**particleS list for the neighbor build*/

/**DISSOCIATION, HOPPING, PHOSPHORYLATION RATES*/
double dissociationRate=0.73;
double blockedDissociationRate=1.35;
double singleParticleDissociationRate=0;
double hoppingRate=0.706;
double phosphorylationRate=1500;

int main(int argc, char** argv)
{
    //++++++++++++++Routine that reads in the factorial values into the array+++++++++++++++++//
    std::ifstream factorialFile("factorial.dat", std::ifstream::in);
    if(factorialFile.good())
    {
        std::string line;
        int c=0;
        while (getline (factorialFile, line) )
        {
            std::stringstream s1(line);
            s1 >> tfactorial[c];
            ++c;                    
        }
    }
    
    //+++++++++++Check if the number of command line arguments is correct++++++++//wr
    if (argc != 4)
    {
        std::cout << "./particleList  <input_file> rng_seed readin_type" << std::endl;
    }
  
    /*LOADING THE SIMULATION PARAMETERS FROM INPUT FILE*/
    prm.loadParameters(argv[1]);
    //prm.printInfo();
    
    /* Create a list of species and load from species file*/
    loadSpeciesInfo ();

    
    /*CALCULATE THE DERIVED PARAMETERS*/
    prm.calculateDerivedParameters();
    /*PRINT THE PARAMETERS*/
   
    
    /*RANDOM NUMBER GENERATOR*/
    /*type of the random number generator*/
    const gsl_rng_type * T;
    T = gsl_rng_mt19937;
    /*initializing instance of rng*/
    rng = gsl_rng_alloc(T);
    /*setting the seed to the rng*/
    gsl_rng_set(rng, atoi(argv[2]) );
    std::cout << "rng : " << gsl_rng_uniform(rng) << std::endl;
     int readin = atoi(argv[3]);/**Gives how particles should be generated...0 if new particles
                                 * placed randomly..1 if particles placed with positions
                                 * entered manually..
                                 * 2 if n particles should be loaded from some configuration file
                                 * 3 for two particle FFS type loading from configuration file/
                                 */
    
    //GENERATE AND LOAD PARTICLES
    generateAndLoadParticles( readin );
     
    //check if domains can be built on the generated particles
    checkIfDomainsCanBeBuiltOnLangevinParticles();
    //after building domains sort the scheduler list
    sortSchedulerList();
    //for visualisation
    int vtkCount=0; int step=0;
    //for keeping track
    int ldSteps=0;
    int gfrdSteps=0;
    double ldTime = 0;
    double gfrdTime = 0;
    
    //+++++++++++++++++++++++++++++++++ROTATING THE PARTICLES+++++++++++++++++++++++++++++++++++++/
    // particleList[0].patchPos = particleList[0].q.bodyToSpace() * speciesList[particleList[0].speciesType].particleRadius + particleList[0].particlePosition;
    myVector axis;
    axis.x=0.0;
    axis.y=0.0;
    axis.z=1.0;
    //particleList[0].q.rotate(-M_PI/2.0, axis);
    myVector axis1;
    axis1.x=0.0;
    axis1.y=1.0;
    axis1.z=0.0;
    particleList[1].q.rotate(-M_PI, axis1);
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++/
    
    
    //++++++++++++++++++++++++++++++MAIN TIME LOOP+++++++++++++++++++++++++++++++++++//
    while ( tSim < prm.tEnd)
    {
        ++step;
        //std:: cout << "tSim = " << tSim << std::endl;
       
        //Calculate the number of BD and GFRD particles in the system
        std::vector<particle>::iterator iter = particleList.begin();
        u_long numberOfLDParticles=0;/**Number of LD particles at any given time of the simulation*/
        u_long numberOfGfrdParticles=0;/**Number of particles in the eGFRD regime*/
        
        while(iter != particleList.end())
        {
            if(iter->particleType == 2)
            {
                ++numberOfLDParticles;
            }
            else if(iter->particleType == 1)
            {
                ++numberOfGfrdParticles;
            }
            ++iter;
        }
        //std::cout << tSim << " time Of Constrn " << particleList[4].timeOfConstructionOfDomain << std::endl;
        numberOfLDParticles = numberOfLDParticles-idleParticleList.size();
        //std::cout << numberOfLDParticles << " " << numberOfGfrdParticles << std::endl;
        if(numberOfLDParticles > 0)
        {
            ++ldSteps;
            ldTime = ldTime + prm.deltaT;
            //std::cout << " LD" << std::endl;
             //keeps track of the particles that associate energy, particleI, particleJ, patchID
            std::vector<std::tr1::tuple<int, int,int> > associationList;
            dissociationBD();
            //std::cout << " dissociation" << std::endl;
            sortSchedulerList();
            //std::cout << " sort" << std::endl;
            integrator ();
            //std::cout << " inte" << std::endl;
            tSim += prm.deltaT;
            updateAssociationList(associationList);
            //std::cout << " up" << std::endl;
            associationReaction(associationList);
            //std::cout << " ass" << std::endl;
            sortSchedulerList();
            //std::cout << " sort" << std::endl;
            createDomains();
            //std::cout << " create" << std::endl;
            
        }
        //std::cout << "rng : " << gsl_rng_uniform(rng) << std::endl;
        numberOfLDParticles=0;
        numberOfGfrdParticles=0;
        iter = particleList.begin();
        while(iter != particleList.end())
        {
            if(iter->particleType == 2)
            {
                ++numberOfLDParticles;
            }
            else if(iter->particleType == 1)
            {
                ++numberOfGfrdParticles;
            }
            ++iter;
        }
        numberOfLDParticles = numberOfLDParticles-idleParticleList.size();
        
        //
        
//        std::cout << numberOfLDParticles << " " << numberOfGfrdParticles << std::endl;
//        std::cout << particleList[schedulerList[0]].nextEventTime << " " <<
//        particleList[schedulerList[1]].nextEventTime << " " << particleList[schedulerList[2]].nextEventTime << " " << tSim<< std::endl;
        if(step%100000 == 0)
        {
            std:: cout << "tSim = " << tSim << " ; ldSteps = " << ldSteps << " ; gfrdSteps = " << gfrdSteps << " ; ldTime: " << ldTime <<" ; gfrdTime " << gfrdTime << std::endl;
            //writeOutput(vtkCount);
            ++vtkCount;
            ldSteps = 0;
            gfrdSteps = 0;
            ldTime = 0;
            gfrdTime = 0;
        }
        //std::cout << step << std::endl;
        if(numberOfGfrdParticles > 0)
        {      
            
            if(
                ((tSim >= particleList[schedulerList[0]].nextEventTime) &&
                (particleList[schedulerList[0]].particleType == 1))||
               (idleParticleList.size() == 2 && particleList.size()==3) ||
               (idleParticleList.size() == 1 && particleList.size()==2)||
               numberOfLDParticles == 0
              )
            {
                //std::cout << particleList[schedulerList[0]].nextEventTime << " " <<
                //particleList[schedulerList[1]].nextEventTime << " " << tSim<< std::endl;
                //particleList[schedulerList[2]].nextEventTime << std::endl;
                 //std::cout << " GFRD" << std::endl;
                if(particleList[schedulerList[0]].nextEventType == 0)
                {
                    ++gfrdSteps;
                    gfrdTime = gfrdTime + std::abs(particleList[schedulerList[0]].nextEventTime-tSim);
                    //std::cout << " escape " << particleList[schedulerList[0]].particleID << std::endl;
                    tSim = particleList[schedulerList[0]].nextEventTime;
                    particleList[schedulerList[0]].burstTheDomain("escape");
                    createDomains();
                }
                else if(particleList[schedulerList[0]].nextEventType == 1)
                {
                    //std::cout << particleList[0].patchState[0] << " patch states " << particleList[0].patchState[1] << std::endl;
                    //std::cout << tSim << " " << particleList[schedulerList[0]].speciesType << " " << particleList[schedulerList[0]].particleID<< " GFRD dissociation" << std::endl;
                    tSim = particleList[schedulerList[0]].nextEventTime;
                    particleList[schedulerList[0]].dissociateParticleFromDomain();
                    //std::cout << particleList[0].patchState[0] << " patch states " << particleList[0].patchState[1] << std::endl;
                    //std::cout << "ipl size : " << idleParticleList.size() << std::endl;
                }
            }
        }
    }
    return 0;
}


