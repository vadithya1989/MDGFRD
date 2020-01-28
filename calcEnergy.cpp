
#include "simulationParameters.h"
#include "myVector.h"
//#include "particle.h"
//
#include<limits>
#include <cmath>
#include<vector>

#include "externDefs.h"

//CHANGE rPATCH FOR MULTIPLE PATCHES!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
double potential(   double r,
                 double a,
                 double x0,
                 double dxstar,
                 double b,
                 double dxc,
                 const simulationParameters &prm)
{
    
    if( (x0-dxstar)<r && (x0+dxstar)>r )
    {
        return (1-a*(r/prm.sigma-x0/prm.sigma)*(r/prm.sigma-x0/prm.sigma));
    }
    else if( (x0/prm.sigma-dxstar)>r && (x0/prm.sigma-dxc)<r )
    {
        return (b*(r/prm.sigma-(x0/prm.sigma-dxc/prm.sigma))*(r/prm.sigma-(x0/prm.sigma-dxc/prm.sigma)) );
    }
    else if( (x0+dxstar)<r && (x0+dxc)>r )
    {
        return (b*(r/prm.sigma-(x0/prm.sigma+dxc/prm.sigma))*(r/prm.sigma-(x0/prm.sigma+dxc/prm.sigma)) );
    }
    else
    {
        return 0.0;
    }
    
    
}


void calcEnergy (myVector &r,
                 std::vector< std::pair<int, myVector> > &rPatchVector,
                 double &E,
                 particle &pa,
                 particle &pb
                 )
{
    double Vrep=0.0;
    double Vatt=0.0;
    double Vtot=0.0;
    double VisoAtt=0.0;
    
    double scale = 20.0*prm.boltzmannConstant*prm.temperature;
    double rfactor = 5;
    double brep = 2.6036;
    double dxcrep = 1.17647*prm.sigma;
    double batt = 5.0;
    double dxcatt = 0.5*prm.sigma;
    
    double arep=1;
    double x0rep=0;
    double dxstarrep = 0.85*prm.sigma;
    
    double aatt = 20;
    double x0att = 0;
    double dxstaratt = 0.1*prm.sigma;
    
    double bisoAtt = 2.6036;
    double dxcisoAtt = 1.17647*prm.sigma;
    double aisoAtt = 1.0;
    double x0isoAtt = 0.4*prm.sigma;
    double dxstarisoAtt = 0.85*prm.sigma;
    
    
    Vrep = scale*rfactor*potential(r.magnitude(), arep, x0rep, dxstarrep, brep, dxcrep, prm);
    VisoAtt = -scale*0.5*potential(r.magnitude(), aisoAtt, x0isoAtt, dxstarisoAtt, bisoAtt, dxcisoAtt, prm);
    for(int i=0; i<rPatchVector.size(); ++i)
    {
        Vatt += -scale*potential(rPatchVector[i].second.magnitude(), aatt, x0att, dxstaratt, batt, dxcatt, prm);
    }
    
    Vtot = Vatt + Vrep + VisoAtt;
    E = Vtot;
}

