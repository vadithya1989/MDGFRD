
#include <vector>

#include "particle.h"

void removeTorqueComponent ()
{
    for(int i=0;i<particleList.size();++i)
    {   
        quaternion torque3d;
        
        double tdotq = (particleList[i].torque.dot(particleList[i].q));
        particleList[i].torque = particleList[i].torque - (particleList[i].q * tdotq);
        tdotq = (particleList[i].torque.dot(particleList[i].q));
        
        torque3d.q0 =  particleList[i].torque.q0*particleList[i].q.q0 + particleList[i].torque.q1*particleList[i].q.q1 + particleList[i].torque.q2*particleList[i].q.q2 + particleList[i].torque.q3*particleList[i].q.q3;
        torque3d.q1 = -particleList[i].torque.q0*particleList[i].q.q1 + particleList[i].torque.q1*particleList[i].q.q0 + particleList[i].torque.q2*particleList[i].q.q3 - particleList[i].torque.q3*particleList[i].q.q2;
        torque3d.q2 = -particleList[i].torque.q0*particleList[i].q.q2 - particleList[i].torque.q1*particleList[i].q.q3 + particleList[i].torque.q2*particleList[i].q.q0 + particleList[i].torque.q3*particleList[i].q.q1;
        torque3d.q3 = -particleList[i].torque.q0*particleList[i].q.q3 + particleList[i].torque.q1*particleList[i].q.q2 - particleList[i].torque.q2*particleList[i].q.q1 + particleList[i].torque.q3*particleList[i].q.q0;
        
        torque3d = torque3d*double(0.5);
        
        myVector torqueSpace;
        torqueSpace.x = torque3d.q1 ;
        torqueSpace.y = torque3d.q2 ;
        torqueSpace.z = torque3d.q3 ;
        
        
        torqueSpace = particleList[i].q.bodyToSpace(torqueSpace);
    }
    
    /*test by printing out torque in cartesian(3d) coordinates*/
                 
    
}
