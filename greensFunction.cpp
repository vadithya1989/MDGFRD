/* 
 * File:   greensFunction.cpp
 * Author: vijaykumar
 * 
 * Created on February 7, 2015, 4:09 PM
 */

#include "greensFunction.h"
#include "findRoot.h"
#include <assert.h>
#include "externDefs.h"

double greensFunction::ellipticTheta4Zero ( double q )
{
    if ( std::abs(q) > 1.0 )
    {
        std::cout << "Invalid q argument into ellipticTheta4Zero" << std::endl;
    }

    const long int N   (1000);
    double value        (1.0);
    double q_n          (q);
    double q_2n         (1.0);

    for (long int n(1); n <=N; ++n)
    {
        const double term2(1.0 - q_2n * q);
        q_2n = q_n * q_n;
        const double term1(1.0 - q_2n);
        const double term(term1 * term2 * term2);
        const double value_prev(value);
        value *= term;

        if (fabs(value - value_prev) < 1e-18)
        {
            return value;   //normal exit
        }

        q_n *= q;
    }

    std::cout << "Elliptictheta4zero did not converge" << std::endl;
    return value;   //return value even if there is no convergence
}

double greensFunction::p_survival(double t) const
{
    const double asq(a * a);
    const double PIsq(M_PI * M_PI);
    const double q(- D * PIsq * t / asq);

    return 1.0 - ellipticTheta4Zero(exp(q));
}

struct p_survival_params
{
    const greensFunction* const gf;
    const double rnd;
};

static double p_survival_F( double t, p_survival_params const* params )
{
    return params->rnd - params->gf->p_survival(t);
}


std::pair<double, int> greensFunction::drawTime(gsl_rng* r, std::vector<species> &speciesList, int speciesType) const
{
    //event type 1 is escape and 2 is dissociation
    //*******************************************GET AN ESCAPE TIME*************************************************//
    double tEscape = 0.0;
    double rnd = gsl_rng_uniform(r);
    if (rnd >= 1.0 || rnd < 0.0)
    {
        std::cout <<  "Invalid randon number into drawTime" << std::endl;
    }

    /*if (D == 0.0 || a == INFINITY)
    {
        tEscape = INFINITY;
        return tEscape; // If the diffusion constant is zero or the domain size is infinite
    }

    if (a == 0.0)
    {
        tEscape = 0.0;
        return tEscape; // If the domain size is zero
    }*/

    //parameters for the function that goes into the root finder
    p_survival_params params = { this, rnd };
    //Function that goes into the root finder brought to the right form
    gsl_function F =
    {
        reinterpret_cast<typeof(F.function)>(&p_survival_F),
        &params
    };

    //the upper and lower bounds for the root finder
    const double t_guess(a * a / (6. * D));
    double low (t_guess);
    double high (t_guess);

    const double value ( GSL_FN_EVAL(&F, t_guess) );
    if (value < 0.0)
    {
        high *= 10.0;

        for ( ; ; )
        {
            const double high_value ( GSL_FN_EVAL(&F, high) );

            if (high_value >= 0.0)
            {
                break;
            }

            if (fabs(high) >= t_guess * 1e6)
            {
                std::cout << "Could not adjust high value in drawTime for the root solver";
            }
            high *= 10.0;
        }
    }
    else
    {
        double low_value_prev(value);
        low *= .1;

        for (;;)
        {
            const double low_value(GSL_FN_EVAL(&F, low));

            if (low_value <= 0.0)
            {
                break;
            }

            if ( fabs(low) <= t_guess * 1e-6 || fabs(low_value - low_value_prev) < CUTOFF )
            {
                std::cout << "could not adjust hi value .. so returning lo value" << std::endl;
                tEscape = low;
                 std::pair<double, int> temp;
                 temp = std::make_pair(tEscape, 0);
                return temp;
            }
            low_value_prev = low_value;
            low *= .1;
        }
    }
    const gsl_root_fsolver_type* solverType(gsl_root_fsolver_brent);
    gsl_root_fsolver* solver(gsl_root_fsolver_alloc(solverType));

    tEscape = (findRoot(F, solver, low, high, 1e-18, 1e-12,
                            "greensFunction::drawTime"));
    std::pair<double, int> escapePair;
    escapePair = std::make_pair(tEscape, 0);
    gsl_root_fsolver_free(solver);
    //*****************************************************************************************************************//

    //**********************************GET A DISSOCIATION TIME********************************************************//
    /*double tDissociation = 0.0;
    tDissociation = (1.0/speciesList[speciesType].dissociationConstant) * (-log(gsl_rng_uniform(r)));
    std::pair<double, int> dissociationPair;
    dissociationPair = std::make_pair(tDissociation, 1);
    
    if( (speciesList[speciesType].dissociationConstant != 0.0) )
    {
        if(escapePair.first <= dissociationPair.first)
            return escapePair;
        else if(dissociationPair.first < escapePair.first)
            return dissociationPair;
        else
        {
            std::pair<double, int> temp;
            temp = std::make_pair(-1, -1);
            return temp;
        }
    }
    else*/
        return escapePair;

    //*******************************************************************************//
}

//*************************************************************************************************************************************************************************************************************//

double greensFunction::p_int_r_free(double r, double t) const
{
    const double Dt(D * t);
    const double sqrtDt(sqrt(Dt));
    const double sqrtPI(sqrt(M_PI));

    return erf(r / (sqrtDt + sqrtDt)) - r * exp(- r * r / (4.0 * Dt)) / (sqrtPI * sqrtDt);
}

double greensFunction::p_int_r(double r, double t) const
{
    double value(0.0);

    const double p_free(this->p_int_r_free(r, t));

    // p_int_r is always smaller than p_free.
    if (fabs(p_free) < CUTOFF)
    {
        return 0.0;
    }

    const double asq(a * a);
    const double PIsq(M_PI * M_PI);

    const double PIr(M_PI * r);
    const double PIr_a(PIr / a);
    const double DtPIsq_asq(D * t * PIsq / asq);

    const double factor(2.0 / (a * M_PI));

    const double maxn((a / M_PI) * sqrt(log(exp(DtPIsq_asq) / CUTOFF) /
                                          (D * t)));

    const long int N_MAX(10000);

    const long int N(std::min(static_cast<long int>(ceil(maxn) + 1),
                               N_MAX));
    if (N == N_MAX)
    {
        std::cout << "p_int_r dint converge" << std::endl;
    }


    for (long int n(1); n <= N; ++n)
    {
        const double term1(exp(- n * n * DtPIsq_asq));

        const double angle_n(n * PIr_a);
        double sin_n;
        double cos_n;
        __sincos(angle_n, &sin_n, &cos_n);
        const double term2(a * sin_n);
        const double term3(n * PIr * cos_n);

        const double term(term1 * (term2 - term3) / n);
        value += term;
    }

    return value * factor;
}

struct p_r_params
{
    const greensFunction* const gf;
    const double t;
    const double target;
};

static double p_r_free_F (double r, p_r_params const* params)
{
    return params->gf->p_int_r_free(r, params->t) - params->target;
}


static double p_r_F (double r, p_r_params const* params)
{
    return params->gf->p_int_r(r, params->t) - params->target;
}


double greensFunction::drawR(double rnd, double t) const
{
    //the uniform random number must be [0,1)
    if (rnd >= 1.0 || rnd < 0.0)
    {
        std::cout << "Invalid random number into drawR function" << std::endl;
    }
    //time sampled must be positive
    if (t < 0.0)
    {        
        std::cout << "invalid time into the drawR function" << std::endl;
    }
    //domain must be greater than zero, time must be > 0 and diffusion constant > 0
    if (a == 0.0 || t == 0.0 || D == 0.0)
    {
        return 0.0;
    }
    const double thresholdDistance(CUTOFF_H * sqrt(6.0 * D * t));

    gsl_function F;
    double psurv;

    if (a <= thresholdDistance)
    {
        //psurv = p_survival(t);  // this causes a problem when p_survival is very small.
        psurv = p_int_r(a, t);
        if (psurv == 0.0)
        {
            std::cout << "warning:: high next event time" << std::endl;
            return a;
        }

        assert(psurv >= 0.0);

        F.function = reinterpret_cast<typeof(F.function)>(&p_r_F);
    }
    else
    {
        // p_int_r < p_int_r_free
        if (p_int_r_free(a, t) < rnd)
        {
            std::cout << "p_int_r_free(a, t) < rnd, returning a" << std::endl;
            return a;
        }

        psurv = 1.0;
        F.function = reinterpret_cast<typeof(F.function)>(&p_r_free_F);
    }
    //std::cout << thresholdDistance/a << " " << a*a / (6*D) << " " << t << std::endl;
    const double target(psurv * rnd);
    p_r_params params = { this, t, target };

    F.params = &params;

    const double low(0.0);
    const double high(a);
    //const double high(std::min(thresholdDistance, a));

    const gsl_root_fsolver_type* solverType(gsl_root_fsolver_brent);
    gsl_root_fsolver* solver(gsl_root_fsolver_alloc(solverType));

    const double r(findRoot(F, solver, low, high, 1e-18, 1e-12,
                            "GreensFunction3DAbsSym::drawR"));

    gsl_root_fsolver_free(solver);

    return r;



}
