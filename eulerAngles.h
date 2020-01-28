/* 
 * File:   eulerAngles.h
 * Author: touldrid
 *
 * Created on 20 November 2015, 16:27
 */

#ifndef EULERANGLES_H
#define	EULERANGLES_H

#include <cmath>
#include <iostream>

class eulerAngles {
public:
    eulerAngles();
    eulerAngles(const eulerAngles& orig);
    virtual ~eulerAngles();
    
     void operator = (const eulerAngles &);
     
     friend std::ostream& operator<< (std::ostream &out, eulerAngles &ea);

    
    double alpha;
    double beta;
    double gamma;
private:

};

#endif	/* EULERANGLES_H */

