/* 
 * File:   particle.h
 * Author: adithya
 *
 * Created on December 15, 2014, 7:23 PM
 */

#ifndef PARTICLE_H
#define	PARTICLE_H

//#include "myVector.h"
//#include "quaternion.h"
//#include "species.h"
//#include "simulationParameters.h"
//
//#include <fstream>
//#include <sstream>
//#include <vector>
//
//#include <gsl/gsl_rng.h>
//#include <gsl/gsl_randist.h>

#include "myVector.h"
#include "quaternion.h"
#include "main.h"

class particle 
{
public:
    particle();
    particle(const particle& orig);
    virtual ~particle();
public:
    int particleID;/**This ID helps identifying the particle in question*/
    int speciesType;/**Is the species type of the particle, which in turn holds information about the particles' physical properties.*/
    int particleType;/**Describes the type of the particle, 1 if a particle with a domain(eGFRD), 2 if LD*/
    //std::vector < std::pair <int, int> > complexList;
    
    std::vector < std::pair <int, int> > complexList;/**Gives the list of all other particles (if any)
                                                      connectected to this particle in a complex
                                                      and the patches of the particle that are occupied*/
    /**
     * @memberVariableOfSpeciesClass is a vector that contains the
     * the state of all patches
     * for a substrate
     * 0 is unphosphorylated
     * 1 is phosphorylated
     * for an enzyme its always 1
     */
    std::vector<int> patchState;
    
    myVector particlePosition;  /**3D position of the centre of the particle*/
    myVector force;             /**force on the particle*/
    
    quaternion q;               /**The orientation of the particle in 3 dimensions*/
    quaternion torque;          /**The torque of the particle*/
    
    std::vector<myVector> patchPos;          /** 3D position of the centre of the patch of the particle which is used to
                                 visualize the rotation*/
    
    /*GREENS FUNCTIONS RELATED VARIABLES*/
    myVector nextEscapePosition;/**The center of the particle when the particle escapes*/
    myVector burstPosition;/**center of the particle when the domain bursts*/
    double R;/**This is the radius of displacement of the particle from the original position, after bursting*/
    
     /*DOMAIN RELATED MEMBER VARIABLES*/
    int nextReactionType;/**if 0 then nothing (just Brownian step)
                            if 1 then monomolecular dissociation
                            if 2 hopping
                            if 3 phosphorylation/dephosphorylation*/
    int nextEventType;/**If 0 its escape(from the domain)
                         if 1 its reaction
                         */
    double nextEscapeTime;/**The time at which the next escape occurs*/
    double nextReactionTime;/**The time at which the next reaction occurs
                             minimum of dissociation time, hopping time and phosphorylation time*/
    double nextEventTime;/**The time at which the next event occurs, the min of the above two values*/
    double timeOfConstructionOfDomain;/**Global time at which the domain is constructed*/
    double domainRadius;/**radius of the domain around the particle*/
    
    bool operator < (const particle temp) const
    {
        return (this->nextEventTime < temp.nextEventTime);
    }
    
     /**Function prints out the particle class
     */
    friend std::ostream& operator<< (std::ostream &out, particle &part);
    
    
    /**
     This function draws a random position at the sphere of domain radius if the particle escapes
     */
    void drawNextEscapePosition();


    /**
     This function draws a new position for the particle if the domain is prematurely burst
     */
    void drawBurstPosition (double time);
    

    /**
     This function gets a next event time and type using greens functions for a single domain
     */
    void drawNextEventTimeAndType ();
    

    
    /**
     This function builds a domain on the particle which calls the function

     @param nearestNeighborDistance radius of the domain to be built
     */
    void buildADomain(double nearestNeighborDistance);
    

    /**
     This function bursts the given domain it can be used for bursting a domain, moving a particle or converting a particle to a LD particle

     @param instruction escape, burst or convert
     */
    void burstTheDomain(std::string instruction);
    

    /**
     This dissociates a particle in a domain
     */
    void dissociateParticleFromDomain();
    
};/**Particle holds all information related to the particle in the simulation*/

#endif	/* PARTICLE_H */

