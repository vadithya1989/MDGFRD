/* 
 * File:   findRoot.h
 * Author: vijaykumar
 *
 * Created on February 7, 2015, 4:12 PM
 */

#ifndef FINDROOT_H
#define	FINDROOT_H

#include <gsl/gsl_roots.h>


double findRoot(gsl_function const& F, gsl_root_fsolver* solver, double low,
              double high, double tol_abs, double tol_rel, char const* funcName);


#endif	/* FINDROOT_H */

