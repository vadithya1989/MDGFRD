//
//  myVector.cpp
//  associationMultiplePatches
//
//  Created by Adithya on 23/02/2017.
//  Copyright © 2017 Adithya. All rights reserved.
//

#include "myVector.h"

#include <cmath>

myVector::myVector()
{
    x = double(0);
    y = double(0);
    z = double(0);
}

myVector::~myVector()
{
    
}

myVector myVector::operator+ (myVector in)
{
    myVector temp;
    temp.x = x + in.x;
    temp.y = y + in.y;
    temp.z = z + in.z;
    return temp;
}

myVector myVector::operator- (myVector in)
{
    myVector temp;
    temp.x = x - in.x;
    temp.y = y - in.y;
    temp.z = z - in.z;
    return temp;
}

myVector myVector::operator* (long double k)
{
    myVector temp;
    temp.x = x*k;
    temp.y = y*k;
    temp.z = z*k;
    return temp;
}

myVector myVector::operator/ (long double k)
{
    myVector temp;
    temp.x = x/k;
    temp.y = y/k;
    temp.z = z/k;
    return temp;
}

void myVector::operator= (long double k)
{
    x = k;
    y = k;
    z = k;
}

void myVector::operator= (myVector in)
{
    x = in.x;
    y = in.y;
    z = in.z;
}

long double myVector::dot (myVector b)
{
    long double temp;
    temp = (x*b.x) + (y*b.y) + (z*b.z);
    return temp;
}

myVector myVector::cross (myVector b)
{
    myVector temp;
    temp.x = y*b.z - z*b.y;
    temp.y = z*b.x - x*b.z;
    temp.z = x*b.y - y*b.x;
    return temp;
}

long double myVector::magnitude ()
{
    long double temp;
    temp = sqrt(x*x + y*y + z*z);
    return temp;
}

std::ostream& operator<< (std::ostream &out, myVector &myVec)
{
    // Since operator<< is a friend of the Point class, we can access
    // Point's members directly.
    std::cout << "(" << myVec.x << ", " <<
    myVec.y << ", " <<
    myVec.z << ")";
    return out;
}
