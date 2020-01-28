/* 
 * File:   main.h
 * Author: adithya
 *
 * Created on December 13, 2014, 12:40 PM
 */

#ifndef MAIN_H
#define	MAIN_H

#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <limits>
#include <tr1/tuple>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include <sys/time.h>
#include <sys/types.h>

#include "externDefs.h"
#include "myVector.h"
#include "simulationParameters.h"
#include "species.h"
#include "particle.h"
#include "eulerAngles.h"


/**WRITE INTO ERROR LOG FILE*/
 inline void writeToErrorLogFile( const std::string &text )
{
    std::stringstream s2;
    std::string filename2;
    s2 << "err.log" ;
    s2 >> filename2;
     std::ofstream errorFile(filename2.c_str(),  std::fstream::app);
    if ( errorFile.good() )
    {
        errorFile << text << std::endl;
    }
}

/**
 * @function
 Function that loads the species.dat file into the species list
 */
void loadSpeciesInfo ();

/**
 * Generates a random set of particles in the simulation world and adds it to the particle list
 @param readin Gives how particles should be generated :
 0 if new particles placed randomly.
 1 if particles placed with positions entered manually.
 2 if n particles should be loaded from some configuration file.
 3 for two particle FFS type loading from configuration file.
 */
void generateAndLoadParticles(int &readin);

/**
 * The main routine of the Langevin dynamics integrator
 */
void integrator ();


/**
 Writes out the VTK file for visualization
 @param vtkCount vtk file count
 */
void writeOutput (
                  int &vtkCount );


/**
 * This function calculates the interparticle distance between two particles 
 * also taking into account periodic boundaries
 * @param r inter particle distance
 * @param posA position of i particle
 * @param posB position of j particle
 */
void measureInterParticleDistace(myVector &r,
                                 myVector &posA, 
                                 myVector &posB);


/**
 * This function measures the distance between two patches of two particles.
 * @param r inter particle vector
 * @param rPatchVector inter patch vector
 * @param pi particle i
 * @param pj particle j
 */
void measureInterPatchDistance( myVector &r, 
                                std::vector< std::pair<int, myVector> > &rPatchVector,
                                particle &pi,
                                particle &pj
                                );

/**
 * This calculates the energy of the system i.e the potential energy
 * @param r inter particle distance
 * @param rPatchVector intter patch distance
 * @param E energy
 * @param pa particle A
 * @param pb particle B
 */
void calcEnergy (   myVector &r, 
                    std::vector< std::pair<int, myVector> > &rPatchVector,
                    double &E,
                    particle& pa,
                    particle& pb);



/**
 * This function adds all the pairs of particles which have energy lesser than a threshold
 * value => pairs that will associate
 * @param associationList list of particles that needs to associate

 */
void updateAssociationList(std::vector< std::tr1::tuple<int, int, int> > &associationList);


/**
 * This function replaces the two associating particles by a single particle and does all the necessary changes
 * @param associationList list of particles that needs to associate
 */
void associationReaction(std::vector< std::tr1::tuple<int, int, int> > &associationList
                         );

/**
 * This function checks if a C, BD particle dissociates and if it does it places the particles with the pre generated configuration
 */
void dissociationBD();

/**
 * This function sorts the particle list and then reassigns the particleID to the particles based on their position in the particle list
 */
void sortSchedulerList();


/**
 * This function initializes all values in the cells list to -1
 */
void blankCells ();

/**
 * This function finds the cell number to which the particle belongs to
 * @param pi the particle whose cell has to be found
 * @return the cell number(ID) of the particle
 */
int findCellOfTheParticle (particle &pi);
/**
 * This function sorts particles into their respective cells
 */
void sortParticles ();

/**
 * This function is called after the langevin particles are propogated. It checks if these particles bursts any domains.
 */
void checkIfDomainsCanBeBuiltOnLangevinParticles();


/**
 * This function creates domain for LD particles after bursting all possible domains locally
 */
void createDomains();

/**
 * This function return the nearest distance of a langevin particle to another langevin particle
 * @param cellID neighbor list
 * @param particleID neighbor list
 * @param instruction giving instructions
 * @return pair of nearest distance and the id of the nearest particle
 */
std::pair<double, int> returnNearestDistance(
                                             int &cellID,
                                             int &particleID,
                                             std::string instruction);





/**
 * Calculates the Wigner small d function
 * @param j to calculate wigner small d
 * @param m to calculate wigner small d
 * @param k to calculate wigner small d
 * @param ea euler angles
 */
double wignerSmalld(double &j,
                    double &m,
                    double &k,
                    eulerAngles &ea);


/**
 * This function return the Greens function for a given in and out orientation
 * @param ea0 intitial euler angles
 * @param ea new euler angles
 * @param D_r is the rotational diffusion const
 * @param t time since building of the domain
 * @return rotational greens function
 */
double rotationalGreensFunction(eulerAngles &ea0,
                                eulerAngles &ea,
                                double &D_r,
                                double &t);


/**
 * Draws angles and checks it with the acceptance pbty and passes it
 * @param ea0 initial euler angles
 * @param D_r rotational diffusion constant
 * @param t difference between the global simulation time and the time of construction of the domain
 * @param rng random number generator
 * @return final orientation
 */
eulerAngles drawAngles(eulerAngles &ea0,
                       double &D_r,
                       double &t,
                        gsl_rng *rng);

/**
 * Rotates the particle based on alpha beta and gamma
 * @param part particle to be rotated
 * @param ea euler angles of the particle
 */
void drawOrientation(particle &part,
                     eulerAngles &ea);

#endif	/* MAIN_H */

