/* 
 * File:   species.h
 * Author: adithya
 *
 * Created on December 13, 2014, 11:28 AM
 */


#ifndef SPECIES_H
#define	SPECIES_H

#include<vector>
#include "myVector.h"

 
/**
 * @class
 * The species class contains the physical properties of the particle
 * the number of particles of the particular species
 * mass in kg
 * diffusion constant in m^2/s
 * radius in m
 * dissociation constant in s^-1
 */
 /**0 is a substrate particle with two patches, both unphosphorylated
 1 is a substrate particle with two patches, one unphosphorylated and one phosphorylated
 2 is a substrate particle with two patches, both phosphorylated
 3 is a kinase    particle with one patch,   which is active
 4 is a phosphotase particle with one patch, which is active
 5 is a kinase    particle with one patch,   which is inactive
 6 is a phosphotase particle with one patch, which is inactive*/

class species 
{
public:
    species();
    species(const species& orig);
    virtual ~species();
public:
    /**
     * @memberVariableOfSpeciesClass is the number of particles of a species type
     */
    int     numberOfParticles;
    
    /**
     * @memberVariableOfSpeciesClass is the number of patches on the particular species
     */
    int numberOfPatches;
    
    /**
     * @memberVariableOfSpeciesClass is a vector that contains the 
     * direction of all the patches in the body frame 
     */
    std::vector<myVector> bodyPatchVectors;
    
    
    /**
     * @memberVariableOfSpeciesClass is the mass of the particle in kg
     */
    double  particleMass;
    /**
     * @memberVariableOfSpeciesClass is the radius of the particle in m
     */
    double  particleRadius;
    /**
     * @memberVariableOfSpeciesClass is the diffusion constant of the particle in m^2/s
     */
    double  diffusionConstant;
    /**
     * @memberVariableOfSpeciesClass is the rotational diffusion constant of the particle in rad^2/s
     */
    double  rotationalDiffusionConstant;
    /**
     * @memberVariableOfSpeciesClass is the dissociation constant of the particle in s^-1
     */
    double  dissociationConstant;
    double  massMomOfInertia;    /**The mass moment of inertia of the particle*/

};

#endif	/* SPECIES_H */

