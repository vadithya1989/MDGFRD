
#include "eulerAngles.h"
#include "main.h"

double rotationalGreensFunction(eulerAngles &ea0, eulerAngles &ea, double &D_r, double &t)
{
    double j=0.0, residual = 1.0;
    double D=0.0;
    while (residual > 1e-8)
    {
        //std::cout << j << " residual: " << residual << std::endl;
        //+++++++++++++++++++For m=0 and k=0++++++++++++++++++++//
        double m_00=0.0, k_00=0.0;
        double D_00 = wignerSmalld(j, m_00, k_00, ea0)*wignerSmalld(j, m_00, k_00, ea);
        //+++++++++++++++++++++++++++++++++++++++++++++++++++++//
        
        //+++++++++++++++++++For m!=0 and k=0++++++++++++++++++++//
        double m_m0=0.0, k_m0=0.0, D_m0=0.0;
        for(m_m0=1.0; m_m0<=j; ++m_m0)
        {
            D_m0 += wignerSmalld(j, m_m0, k_m0, ea0)*wignerSmalld(j, m_m0, k_m0, ea)
                    *double(2)*cos(m_m0*(ea.alpha-ea0.alpha));
        }
        //+++++++++++++++++++++++++++++++++++++++++++++++++++++//
        
        //+++++++++++++++++++For m=0 and k!=0++++++++++++++++++++//
        double m_0k=0.0, k_0k=0.0, D_0k=0.0;
        for(k_0k=1.0; k_0k<=j; ++k_0k)
        {
            D_0k += wignerSmalld(j, m_0k, k_0k, ea0)*wignerSmalld(j, m_0k, k_0k, ea)
                    *double(2)*cos(k_0k*(ea.gamma-ea0.gamma));
        }
        //+++++++++++++++++++++++++++++++++++++++++++++++++++++//
        
        //+++++++++++++++++++For m>0 and k!=0++++++++++++++++++++//
        double m_pmk=0.0, k_pmk=0.0, D_pmk=0.0;
        for(m_pmk=1.0; m_pmk<=j; ++m_pmk)
        {
            for(k_pmk=1.0; k_pmk<=j; ++k_pmk)
            {               
                D_pmk += wignerSmalld(j, m_pmk, k_pmk, ea0)*wignerSmalld(j, m_pmk, k_pmk, ea)
                            * double(2)*cos( m_pmk*(ea.alpha-ea0.alpha) + k_pmk*(ea.gamma-ea0.gamma));
            }
        }
        //+++++++++++++++++++++++++++++++++++++++++++++++++++++//
        
         //+++++++++++++++++++For m<0 and k!=0++++++++++++++++++++//
        double m_nmk=0.0, k_nmk=0.0, D_nmk=0.0;
        for(m_nmk=-1.0; m_nmk>=-j; --m_nmk)
        {
            for(k_nmk=1.0; k_nmk<=j; ++k_nmk)
            {
            
                D_nmk += wignerSmalld(j, m_nmk, k_nmk, ea0)*wignerSmalld(j, m_nmk, k_nmk, ea)
                        * double(2)*cos( k_nmk*(ea.gamma-ea0.gamma)+m_nmk*(ea.alpha-ea0.alpha) );
            }
        }
        //+++++++++++++++++++++++++++++++++++++++++++++++++++++//
        
        double D_prev = D;
        D += (D_00 + D_m0 + D_0k + D_pmk + D_nmk) * (2.0*j+1.0)/(8.0*M_PI*M_PI) * 
                    exp(-D_r*j*(j+1)*t);
        residual = std::abs((D_prev-D)/D) ;
        ++j;
        
    }
    return D;

}
