/* 
 * File:   eulerAngles.cpp
 * Author: touldrid
 * 
 * Created on 20 November 2015, 16:27
 */

#include "eulerAngles.h"

eulerAngles::eulerAngles() 
{
    alpha = double(0);
    beta  = double(M_PI/2.0);
    gamma = double(0);
    
}

eulerAngles::eulerAngles(const eulerAngles& orig) 
{
    alpha = orig.alpha;
    beta  = orig.beta;
    gamma = orig.gamma;
}

void eulerAngles::operator= (const eulerAngles &orig)
{
   alpha = orig.alpha;
   beta  = orig.beta;
   gamma = orig.gamma;
}

std::ostream& operator<< (std::ostream &out, eulerAngles &ea)
{
    // Since operator<< is a friend of the Point class, we can access
    // Point's members directly.
    out <<  ea.alpha << " " <<
            ea.beta << " " <<
            ea.gamma << std::endl;
    return out;
}



eulerAngles::~eulerAngles() 
{
}

