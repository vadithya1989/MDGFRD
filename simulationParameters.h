/* 
 * File:   simulationParameters.h
 * Author: adithya
 *
 * Created on December 15, 2014, 1:34 PM
 */

#ifndef SIMULATIONPARAMETERS_H
#define	SIMULATIONPARAMETERS_H

#include <string>
#include <vector>

#include <sstream>
#include <fstream>
#include <iostream>

class simulationParameters {
public:
    simulationParameters();
    simulationParameters(const simulationParameters& orig);
    virtual ~simulationParameters();
    
    /**
     * Function reads the input and loads it to the simulation parameters object
     * @param filename name of the input file
     */
    void loadParameters( const std::string filename);
    /**
     * Function simply prints out the values of the simulation parameters
     */
    void printInfo();
    /**
     * Function calculates the derived parameters of the simulation parameters object
     */
    void calculateDerivedParameters();
public:
                             /*THE INDEPENDENT PARAMETERS*/
    /*The parameters of the world*/
    double xmin;            /**The min value of the size simulation world in x direction in m*/
    double xmax;            /**The max value of the size simulation world in x direction in m*/
    double ymin;            /**The min value of the size simulation world in y direction in m*/
    double ymax;            /**The max value of the size simulation world in y direction in m*/
    double zmin;            /**The min value of the size simulation world in z direction in m*/
    double zmax;            /**The max value of the size simulation world in z direction in m*/
    int xcells;             /**The number of cells for neighbor build in the x direction*/
    int ycells;             /**The number of cells for neighbor build in the y direction*/
    int zcells;             /**The number of cells for neighbor build in the z direction*/
    double temperature;     /**The ambient temperature of the simulation world in K*/
    double boltzmannConstant;/**The value of Boltzmann constant in kg.m^2.s^-2.k^-1*/
    
    /*The parameters of the particles*/
    std::string patchInputFile;    /**The name of the file that contains the body frame directions of all the patches in a particular species*/
    std::string speciesInputFile;   /**The name of the file that contains information about the species' of particles used in the simulation*/
    std::string configWriteFile;    /**File to which configuration of particles are written to after/whenever necessary*/
    std::string configReadFile;     /**File from which the initial configuration is read from */
    int totalNumberOfSpecies;          /**The total number of different species in the simulation*/
    
    /*The parameters of Lennard-Jones potential*/
    double deltaT;          /**The time step length of the Langevin dynamics integrator in s*/
    double tEnd;            /**The time at which the simulation ends in s*/
    double depthOfThePotentialWell;/**The pre-factor to n*k_B*T which is used to calculate epsilon*/
    double sigmaPreFactor;/**The value "n" that goes into the calculation of sigma, sigma = n*particleRadius*/
    double reactionDistancePreFactor;/**The value "n" that goes into the calculation of reactionDistance, reactionDistance=n*sigma;*/
    
                               /*THE DERIVED PARAMETERS*/
    /*The parameters of the world*/
    double h;               /**The size of the cell of the linked cell neighbor list in m*/
    int totalNumberOfCells; /**The total number of linked cells in all the three directions put together*/
    
    /*The parameters of Lennard-Jones potential*/
    double sigma;           /**The distance between the particles at which the potential is zero in m*/
    double epsilon;         /**The potential depth of the well in Joules*/
    double gamma;           /**The damping coefficient in s^-1*/
    double gammaRot;        /**The rotational damping constant in s^-1*/
    double reactionDistance;/**The peak in the free energy curve, the particle position which decides if the particle is bound or not*/
    double minDomainRadius;  /**The min allowed size of domain....length scale of sigma/2*/
    double rc;                 /**The distance where LJ potential is cut off*/
    
};/**Class which holds all the parameters that govern the simulation*/

#endif	/* SIMULATIONPARAMETERS_H */

