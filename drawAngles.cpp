
#include "eulerAngles.h"
#include "main.h"

eulerAngles drawAngles(eulerAngles &ea0, double &D_r, double &t, gsl_rng *rng)
{
    double p0=rotationalGreensFunction(ea0, ea0, D_r, t);
    while(1)
    {
        eulerAngles ea;
        ea.alpha = gsl_ran_flat (rng, 0.0, 2.0*M_PI);
        ea.beta  = gsl_ran_flat (rng, 0.0,     M_PI);
        ea.gamma = gsl_ran_flat (rng, 0.0, 2.0*M_PI);
        double earlyRejectionFactor = 1.0;
        

        
         /*if(D_r*t < 0.01)
        {
            if((ea.gamma>0.4 && ea.gamma<2.0*M_PI-0.4) || (ea.alpha>0.4 && ea.alpha<2.0*M_PI-0.4)
                    || (abs(ea.beta-M_PI/2.0)>0.4))
            {
                earlyRejectionFactor = 0.1;
            }
        }
        
        else if(D_r*t < 0.05)
        {
            if((ea.gamma>0.7 && ea.gamma<2.0*M_PI-0.7) || (ea.alpha>0.7 && ea.alpha<2.0*M_PI-0.7)
                    || (abs(ea.beta-M_PI/2.0)>0.7))
            {
                earlyRejectionFactor = 0.1;
            }
        }*/
        if(D_r*t <= 0.06)
        {
            if((ea.gamma>1.1 && ea.gamma<2.0*M_PI-1.1) || (ea.alpha>1.1 && ea.alpha<2.0*M_PI-1.1)
               || (std::abs(ea.beta-M_PI/2.0)>1))
            {
                earlyRejectionFactor = 0.01;
            }
            else if((ea.gamma>0.75 && ea.gamma<2.0*M_PI-0.75) || (ea.alpha>0.75 && ea.alpha<2.0*M_PI-0.75)
                    || (std::abs(ea.beta-M_PI/2.0)>0.75))
            {
                earlyRejectionFactor = 0.1;
            }
        }
        else if(D_r*t <= 0.075)
        {
            if((ea.gamma>1.2 && ea.gamma<2.0*M_PI-1.2) || (ea.alpha>1.2 && ea.alpha<2.0*M_PI-1.2)
               || (std::abs(ea.beta-M_PI/2.0)>1.2))
            {
                earlyRejectionFactor = 0.01;
            }
            else if((ea.gamma>0.8 && ea.gamma<2.0*M_PI-0.8) || (ea.alpha>0.8 && ea.alpha<2.0*M_PI-0.8)
                    || (std::abs(ea.beta-M_PI/2.0)>0.8))
            {
                earlyRejectionFactor = 0.1;
            }
        }
        
        else if(D_r*t <= 0.1)
        {
            if((ea.gamma>1.4 && ea.gamma<2.0*M_PI-1.4) || (ea.alpha>1.4 && ea.alpha<2.0*M_PI-1.4)
               || (std::abs(ea.beta-M_PI/2.0)>1.25))
            {
                earlyRejectionFactor = 0.01;
            }
            else if((ea.gamma>1.0 && ea.gamma<2.0*M_PI-1.0) || (ea.alpha>1.0 && ea.alpha<2.0*M_PI-1.0)
                    || (std::abs(ea.beta-M_PI/2.0)>0.9))
            {
                earlyRejectionFactor = 0.1;
            }
        }
        
        else if(D_r*t <= 0.2)
        {
            if((ea.gamma>1.5 && ea.gamma<2.0*M_PI-1.5) || (ea.alpha>1.5 && ea.alpha<2.0*M_PI-1.5)
                    || (std::abs(ea.beta-M_PI/2.0)>1.2))
            {
                earlyRejectionFactor = 0.1;
            }
        }
        else if(D_r*t<=0.4){
            if((ea.gamma>2 && ea.gamma<2.0*M_PI-2) || (ea.alpha>2 && ea.alpha<2.0*M_PI-2))
            {
                earlyRejectionFactor = 0.1;
            }
        }
        
        double Ner =  gsl_ran_flat (rng, 0.0, 1.0);
        
        if(earlyRejectionFactor > Ner)
        {
            double p_acc = sin(ea.beta)*rotationalGreensFunction(ea0, ea, D_r, t)/p0;
            double N =  gsl_ran_flat (rng, 0.0, 1.0);
            if(p_acc/earlyRejectionFactor>1.0)
            {
                writeToErrorLogFile("GFRD rotations rejecting more");
            }
            if (N<p_acc/earlyRejectionFactor)
            {
                return ea;
            }
        }
    }
}

