/* 
 * File:   quaternion.cpp
 * Author: adithya
 * 
 * Created on December 15, 2014, 6:37 PM
 */

#include "quaternion.h"
#include <cmath>


quaternion::quaternion()
{
    q0 = double(0);
    q1 = double(0);
    q2 = double(0);
    q3 = double(0);
}

quaternion::quaternion(double a0, double a1, double a2, double a3)
{
    q0 = a0;
    q1 = a1;
    q2 = a2;
    q3 = a3;
}

quaternion::~quaternion()
{
    //dtor
}


quaternion quaternion::operator* (quaternion q)
{
    quaternion temp;
    temp.q0 = (q0*q.q0) - (q1*q.q1) - (q2*q.q2) - (q3*q.q3);
    temp.q1 = (q0*q.q1) + (q1*q.q0) + (q2*q.q3) - (q3*q.q2);
    temp.q2 = (q0*q.q2) - (q1*q.q3) + (q2*q.q0) + (q3*q.q1);
    temp.q3 = (q0*q.q3) + (q1*q.q2) - (q2*q.q1) + (q3*q.q0);
    return temp;
}

quaternion quaternion::operator* (double k)
{
    quaternion temp;
    temp.q0 = q0*k;
    temp.q1 = q1*k;
    temp.q2 = q2*k;
    temp.q3 = q3*k;
    return temp;
}

quaternion quaternion::operator/ (double k)
{
    quaternion temp;
    temp.q0 = q0/k;
    temp.q1 = q1/k;
    temp.q2 = q2/k;
    temp.q3 = q3/k;
    return temp;
}

quaternion quaternion::operator+ (quaternion in)
{
    quaternion temp;
    temp.q0 = q0 + in.q0;
    temp.q1 = q1 + in.q1;
    temp.q2 = q2 + in.q2;
    temp.q3 = q3 + in.q3;
    return temp;
}

quaternion quaternion::operator- (quaternion in)
{
    quaternion temp;
    temp.q0 = q0 - in.q0;
    temp.q1 = q1 - in.q1;
    temp.q2 = q2 - in.q2;
    temp.q3 = q3 - in.q3;
    return temp;
}

void quaternion::operator= (quaternion in)
{
   q0 = in.q0;
   q1 = in.q1;
   q2 = in.q2;
   q3 = in.q3;
}

void quaternion::operator= (double in)
{
   q0 = in;
   q1 = in;
   q2 = in;
   q3 = in;
}

double quaternion::dot (quaternion b)
{
   double temp;
   temp = (q0*b.q0) + (q1*b.q1) + (q2*b.q2) + (q3*b.q3);
   return temp;
}

void quaternion::rotate(double angle, myVector axis) //for rotating the quaternion about an axis along an angle
{
    axis = axis/axis.magnitude();
    
    double cosThetaByTwo = cos(angle /double(2));
    double sinThetaByTwo = sin(angle /double(2));

    quaternion temp;
    temp.q0 = cosThetaByTwo;
    temp.q1 = axis.x * sinThetaByTwo;
    temp.q2 = axis.y * sinThetaByTwo;
    temp.q3 = axis.z * sinThetaByTwo;

    quaternion temp1;
    temp1.q0 = q0*temp.q0 - q1*temp.q1 - q2*temp.q2 - q3*temp.q3;
    temp1.q1 = q0*temp.q1 + q1*temp.q0 + q2*temp.q3 - q3*temp.q2;
    temp1.q2 = q0*temp.q2 + q2*temp.q0 + q3*temp.q1 - q1*temp.q3;
    temp1.q3 = q0*temp.q3 + q3*temp.q0 + q1*temp.q2 - q2*temp.q1;

    q0 = temp1.q0;
    q1 = temp1.q1;
    q2 = temp1.q2;
    q3 = temp1.q3;
}

myVector quaternion::bodyToSpace()          // to convert from the body frame to the space frame with the body frame as (1 0 0)
{
    myVector temp;
    temp.x = q0*q0 + q1*q1 - q2*q2 - q3*q3;
    temp.y = double(2) * (q1*q2 + q0*q3);
    temp.z = double(2) * (q1*q3 - q0*q2);

    double mag = temp.magnitude();
    temp.x = temp.x/mag;
    temp.y = temp.y/mag;
    temp.z = temp.z/mag;

    //My name is code, Dr code!
    return temp;
}

myVector quaternion::bodyToSpace(myVector &input)          // to convert from the body frame to the space frame with the body frame as (1 0 0)
{
    myVector temp;
    temp.x = (q0*q0 + q1*q1 - q2*q2 - q3*q3) * input.x + double(2) * (q1*q2 - q0*q3)     * input.y + double(2) * (q1*q3 + q0*q2)     * input.z;
    temp.y = double(2) * (q1*q2 + q0*q3)     * input.x + (q0*q0 - q1*q1 + q2*q2 - q3*q3) * input.y + double(2) * (q2*q3 - q0*q1)     * input.z;
    temp.z = double(2) * (q1*q3 - q0*q2)     * input.x + double(2) * (q2*q3 + q0*q1)     * input.y + (q0*q0 - q1*q1 - q2*q2 + q3*q3) * input.z;

    double mag = temp.magnitude();
    temp.x = temp.x/mag;
    temp.y = temp.y/mag;
    temp.z = temp.z/mag;
    
    //My name is code, Dr code!
    return temp;
}

std::ostream& operator<< (std::ostream &out, quaternion &quat)
{
    // Since operator<< is a friend of the Point class, we can access
    // Point's members directly.
    std::cout << "(" << quat.q0 << ", " <<
        quat.q1 << ", " <<
        quat.q2 << ", " << quat.q3 << ")";
    return out;
}

long double quaternion::magnitude ()
{
   long double temp;
   temp = sqrt(q0*q0 + q1*q1 + q2*q2 + q3*q3);
   return temp;
}

