#include "myVector.h"
#include "simulationParameters.h"
#include "externDefs.h"
#include <cmath>

void measureInterParticleDistace(myVector &r,
                                 myVector &posA, 
                                 myVector &posB)
{
    r = posB - posA;

    //applying periodic boundary condition
    double delx = r.x;
    double dely = r.y;
    double delz = r.z;
    myVector rPb;///distance after periodic boundary condition is applied
    rPb.x = r.x - ((floor((delx/prm.xmax)+0.5)) * prm.xmax );
    rPb.y = r.y - ((floor((dely/prm.xmax)+0.5)) * prm.xmax );
    rPb.z = r.z - ((floor((delz/prm.xmax)+0.5)) * prm.xmax );
    
    r = rPb;
}
