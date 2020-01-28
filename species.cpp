/* 
 * File:   species.cpp
 * Author: adithya
 * 
 * Created on December 13, 2014, 11:28 AM
 */

#include "species.h"

species::species() /*constructor*/
{
    numberOfParticles   = 0;
    numberOfPatches     = 0;
    particleMass        = double(0);
    particleRadius      = double(0);
    diffusionConstant   = double(0);
    rotationalDiffusionConstant = double(0);
    dissociationConstant= double(0);
    massMomOfInertia    = double(0);

}

species::species(const species& orig) /*copy constructor*/
{
    numberOfParticles   = orig.numberOfParticles;
    numberOfPatches     = orig.numberOfPatches;
    particleMass        = orig.particleMass;
    particleRadius      = orig.particleRadius;
    diffusionConstant   = orig.diffusionConstant;
    rotationalDiffusionConstant = orig.rotationalDiffusionConstant;
    dissociationConstant= orig.dissociationConstant;
    massMomOfInertia    = orig.massMomOfInertia;
    bodyPatchVectors    = orig.bodyPatchVectors;
}

species::~species() /*destructor*/
{
}

