/* 
 * File:   greensFunction.h
 * Author: vijaykumar
 *
 * Created on February 7, 2015, 4:09 PM
 */

#ifndef GREENSFUNCTION_H
#define	GREENSFUNCTION_H

#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_trig.h>
#include <gsl/gsl_sum.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_sf_expint.h>
#include <gsl/gsl_sf_elljac.h>
#include <gsl/gsl_roots.h>
#include "species.h"


class greensFunction
{
    public:
        greensFunction ( double DD, double aa ) : D(DD), a(aa) {}
        virtual ~greensFunction () {}

        const double D;   //diffusion constant
        const double a;   //domain size

        static constexpr double CUTOFF   = 1e-10;
        static constexpr double CUTOFF_H = 6.0;

    public:

        //*****************************FUNCTIONS TO DRAW TIME*******************************************************************//
        std::pair<double, int> drawTime    ( gsl_rng *, std::vector<species> &, int ) const; //samples the next event time and type that is the minimum between dissociation and escape
        double p_survival                   ( double ) const; //caculates the survival probability which is the integral of greens function over the whole domain
        static double ellipticTheta4Zero   ( double ); //helper function to find the survival probability
        //***********************************************************************************************************************//

        //*****************************FUNCTIONS TO DRAW r***********************************************************************//
         double drawR           ( double, double ) const; //samples the position of the particle
         double p_int_r         (double r, double t) const;
         double p_int_r_free    (double r, double t) const;
         double p_r_fourier     (double r, double t) const;
         //**********************************************************************************************************************//


};



#endif	/* GREENSFUNCTION_H */

/* 
 * File:   greensFunction.h
 * Author: vijaykumar
 *
 * Created on December 24, 2015, 1:42 PM
 */

#ifndef GREENSFUNCTION_H
#define	GREENSFUNCTION_H

class greensFunction {
public:
    greensFunction();
    greensFunction(const greensFunction& orig);
    virtual ~greensFunction();
private:

};

#endif	/* GREENSFUNCTION_H */
