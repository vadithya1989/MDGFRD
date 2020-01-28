
#include "eulerAngles.h"
#include "particle.h"

void drawOrientation(particle &part, eulerAngles &ea)
{  
    myVector axis;
    
    //rotate around y axis by -90
    
    axis.x=0.0;
    axis.y=1.0;
    axis.z=0.0;
    axis = part.q.bodyToSpace(axis);
    part.q.rotate(-M_PI/2.0, axis);
    
    
    myVector yaxis, temp;
    yaxis.x=0.0;
    yaxis.y=1.0;
    yaxis.z=0.0;
    yaxis = part.q.bodyToSpace(yaxis);
   
    myVector zaxis;
    zaxis.x=0.0;
    zaxis.y=0.0;
    zaxis.z=1.0;
    zaxis = part.q.bodyToSpace(zaxis);
     
    part.q.rotate(ea.gamma, zaxis);
    part.q.rotate(ea.beta,  yaxis);
    part.q.rotate(ea.alpha, zaxis);  
    
}

