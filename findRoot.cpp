#include <gsl/gsl_roots.h>
#include <gsl/gsl_errno.h>
#include <iostream>
#include "findRoot.h"

double findRoot ( gsl_function const& F,
                   gsl_root_fsolver* solver,
                   double low,
                   double high,
                   double tol_abs,
                   double tol_rel,
                   char const* funcName )

{
    // low and high should constitute an interval that straddles a root.
    double l(low);
    double h(high);

    gsl_root_fsolver_set(solver, const_cast<gsl_function*>(&F), l, h);

    const unsigned int maxIter(100);

    unsigned int i(0);

    for ( ; ; )
    {
        // iterate
        gsl_root_fsolver_iterate(solver);

        // get found bracketing interval
        l = gsl_root_fsolver_x_lower(solver);
        h = gsl_root_fsolver_x_upper(solver);

        // see if this is acceptable
        const int status(gsl_root_test_interval(l, h, tol_abs,
                                                  tol_rel));

        // stop finder if convergence or too much iterations
        if (status == GSL_CONTINUE)
        {
            if (i >= maxIter)
            {
                gsl_root_fsolver_free(solver);
                std::cout << "root finder failed to converge" << funcName << std::endl;
            }
        }
        else
        {
            break;
        }

        ++i;
    }
    const double root(gsl_root_fsolver_root(solver));

    return root;
}


