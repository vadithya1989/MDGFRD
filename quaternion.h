/* 
 * File:   quaternion.h
 * Author: adithya
 *
 * Created on December 15, 2014, 6:37 PM
 */

#ifndef QUATERNION_H
#define	QUATERNION_H

#include "myVector.h"

///@class Class that stores a 4D vector that can be used to store orientation or torques
class quaternion {
public:
    quaternion();
    quaternion(double, double, double, double);
    ~quaternion();
    
    double q0;
    double q1;
    double q2;
    double q3;
    
    /**
     Overloaded operator for quaternion-quaternion multiplication
     
     @param q qusternion to be multiplied with
     @return qusternion
     */
    quaternion operator * (quaternion q);
    
    /**
     Overloaded operator for quaternion-scalar multiplication
     
     @param scalar number
     @return quaternion
     */
    quaternion operator * (double scalar);
    
    /**
     Overloaded operator for quaternion-scalar division
     
     @param scalar number
     @return quaternion
     */
    quaternion operator / (double scalar);
    
    /**
     Overloaded operator for quaternion addition
     
     @param q quaternion
     @return quaternion
     */
    quaternion operator + (quaternion q);
    
    /**
     Overloaded operator for quaternion subtraction
     
     @param q quaternion
     @return quaternion
     */
    quaternion operator - (quaternion q);
    
    /**
     copies the input quaternion to the current one
     
     @param q quaternion
     */
    void operator = (quaternion q);
    
    /**
     copies the input double value to the all components of the current quaternion
     
     @param scalar number
     */
    void operator = (double scalar);
    
    /**
     calculates the dot produce of two quaternion
     
     @param q quaternion
     @return double
     */
    double dot(quaternion q);
    
    /**
     rotates a quaternion by a given angle about an axis
     
     @param angle angle in radians about which the particle is to be rotated
     @param axis axis about which the quaternion is to be rotated
     */
    void rotate(double angle, myVector axis);
    
    /**
     converting the body vector to the space frame the body vector being (1,0,0)
     
     @return vector in the space frame
     */
    myVector bodyToSpace();
    
    /**
     converting the body vector to the space frame the body vector being (x,y,z)
     
     @param bodyVector body vector
     @return vector in the space frame
     */
    myVector bodyToSpace(myVector &bodyVector);
    
    /**
     Prints out all the quaternion components
     */
    friend std::ostream& operator<< (std::ostream &out, quaternion &quat);
    
    /**
     Calculates the magnitude of the quaternion
     
     @return magnitude
     */
    long double magnitude();
};/**Quaternion is a user defined 4D vector useful to simulate rotations*/

#endif	/* QUATERNION_H */

