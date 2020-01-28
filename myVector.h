/* 
 * File:   myVector.h
 * Author: adithya
 *
 * Created on December 15, 2014, 6:07 PM
 */

#ifndef MYVECTOR_H
#define	MYVECTOR_H
#include <iostream>

///@class myVector is a user defined class to store 3D positions, forces etc
class myVector {
    
public:
    myVector();
    ~myVector();
    
    /**
     Overloaded operator for vector addition
     
     @param vec vector to be added to the current vector
     @return sum of the two vector
     */
    myVector operator + (myVector vec);
    
    /**
     Overloaded operator for vector subtraction
     
     @param vec vector to be subtracted from the current vector
     @return difference of the vector
     */
    myVector operator - (myVector vec);
    
    /**
     Overloaded operator for scalar multiplication
     
     @param scalar number that scales the vector
     @return scaled vector
     */
    myVector operator * (long double scalar);
    
    /**
     Overloaded operator for scalar division
     
     @param scalar number that divises the vector
     @return scaled vector
     */
    myVector operator / (long double scalar);
    
    
    /**
     Calculates cross product
     
     @param vec vector with which cross product is calculated
     @return cross product
     */
    myVector cross (myVector vec);
    
    
    /**
     Equates all the components of vector to the scalar number
     
     @param scalar value to be equated to vector
     */
    void operator = (long double scalar);
    
    /**
     Equates two vectors
     
     @param vec vector to be equated with
     */
    void operator = (myVector vec);
    
    /**
     Calculates the dot product
     
     @param vec vector
     @return scalar dot product
     */
    long double dot(myVector vec);
    
    /**
     Magnitude of the vector
     
     @return scalar magnitude
     */
    long double magnitude();
    
    
    /**
     Prints out all the components of the vector
     
     @param out output stream
     @param myVec vector
     */
    friend std::ostream& operator<< (std::ostream &out, myVector &myVec);
    
public:
    double x;
    double y;
    double z;
    
};

/** myVector class defines a user defined 3D vector*/

#endif	/* MYVECTOR_H */

