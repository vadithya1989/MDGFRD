/* 
 * File:   simulationParameters.cpp
 * Author: adithya
 * 
 * Created on December 15, 2014, 1:34 PM
 */

#include "simulationParameters.h"


simulationParameters::simulationParameters() {
    xmin    = double(0);
    xmax    = double(0);
    ymin    = double(0);
    ymax    = double(0);
    zmin    = double(0);
    zmax    = double(0);
    xcells  = 0;
    ycells  = 0;
    zcells  = 0;
    temperature = double(0);
    boltzmannConstant = double(0);
    
    speciesInputFile = std::string();
    patchInputFile = std::string();
    configWriteFile = std::string();
    configReadFile = std::string();
    totalNumberOfSpecies = 0;
    
    deltaT = double(0);
    tEnd = double(0);
    depthOfThePotentialWell = double(0);
    sigmaPreFactor = double(0);
    reactionDistancePreFactor = double(0);
    
    h = double(0);
    totalNumberOfCells = 0;
    
    sigma = double(0);
    epsilon = double(0);
    gamma = double(0);
    gammaRot = double(0);
    reactionDistance = double(0);
    minDomainRadius = double(0);
    rc = double(0);
}

simulationParameters::simulationParameters(const simulationParameters& orig) {
}

simulationParameters::~simulationParameters() {
}

/**
 * Additional function to read in from the input file
 * @param line line that has to be compared
 * @return the type of the parameter read in from the file
 */
template<typename T>
T getToken( const std::string& line )
{
    T value;
    std::stringstream ss (line);
    ss >> value;
    return value;
}

void simulationParameters::loadParameters( const std::string filename )
{
    std::ifstream infile( filename.c_str(), std::ifstream::in );
    if( infile.good() )
    {
        std::string line;
        while (getline (infile, line) )
        {
            std::string tmp, value;
            std::stringstream ss (line);
            ss >> tmp >> value;
            
            if (tmp.compare("xmin") == 0)
            {
                xmin = getToken<double>(value);
            }
            else if (tmp.compare("xmax") == 0)
            {
                xmax = getToken<double>(value);
            }
            else if (tmp.compare("ymin") == 0)
            {
                ymin = getToken<double>(value);
            }
            else if (tmp.compare("ymax") == 0)
            {
                ymax = getToken<double>(value);
            }
            else if (tmp.compare("zmin") == 0)
            {
                zmin = getToken<double>(value);
            }
            else if (tmp.compare("zmax") == 0)
            {
                zmax = getToken<double>(value);
            }
            else if (tmp.compare("xcells") == 0)
            {
                xcells = getToken<int>(value);
            }
            else if (tmp.compare("ycells") == 0)
            {
                ycells = getToken<int>(value);
            }
            else if (tmp.compare("zcells") == 0)
            {
                zcells = getToken<int>(value);
            }
            else if (tmp.compare("temperature") == 0)
            {
                temperature = getToken<double>(value);
            }
            else if (tmp.compare("boltzmannConstant") == 0)
            {
                boltzmannConstant = getToken<double>(value);
            }
            else if (tmp.compare("patchInputFile") == 0)
            {
                patchInputFile = getToken<std::string>(value);
            }
            else if (tmp.compare("configWriteFile") == 0)
            {
                configWriteFile = getToken<std::string>(value);
            }
            else if (tmp.compare("configReadFile") == 0)
            {
                configReadFile = getToken<std::string>(value);
            }            
             else if (tmp.compare("speciesInputFile") == 0)
            {
                speciesInputFile = getToken<std::string>(value);
            }
             else if (tmp.compare("totalNumberOfSpecies") == 0)
            {
                totalNumberOfSpecies = getToken<int>(value);
            }
            else if (tmp.compare("deltaT") == 0)
            {
                deltaT = getToken<double>(value);
            }
            else if (tmp.compare("tEnd") == 0)
            {
                tEnd = getToken<double>(value);
            }
            else if (tmp.compare("depthOfThePotentialWell") == 0)
            {
                depthOfThePotentialWell = getToken<double>(value);
            }
            else if (tmp.compare("sigmaPreFactor") == 0)
            {
                sigmaPreFactor = getToken<double>(value);
            }
            else if (tmp.compare("reactionDistancePreFactor") == 0)
            {
                reactionDistancePreFactor = getToken<double>(value);
            }
            else if (tmp.compare("gamma") == 0)
            {
                gamma = getToken<double>(value);
            }
            else if (tmp.compare("rc") == 0)
            {
                rc = getToken<double>(value);
            }
            else if (tmp.compare("minDomainRadius") == 0)
            {
                minDomainRadius = getToken<double>(value);
            }

        }
    }
}

void simulationParameters::calculateDerivedParameters()
{
    h = xmax/xcells;
    totalNumberOfCells = xcells*ycells*zcells;
    sigma = 5e-9;
    //sigma = sigmaPreFactor * speciesList[0].particleRadius;
    epsilon = depthOfThePotentialWell * boltzmannConstant * temperature;
   // epsilon = 4e-20*5.0;
    //gamma = boltzmannConstant * temperature / ( speciesList[0].particleMass * speciesList[0].diffusionConstant);
    gamma = 1e14;
    gammaRot = 1e15;
    reactionDistance = reactionDistancePreFactor * sigma;
    rc = 1.6*sigma;
    minDomainRadius = sigma/2.0;
}

void simulationParameters::printInfo()
{
    std::cout << std::endl;
    std::cout << "/++++++++++++++++++++THE INDEPENDENT PARAMETERS++++++++++++++++++++/"     << std::endl;
    std::cout << "/++++++++++++++++++++Parameters of the world++++++++++++++++++++/"     << std::endl;
    std::cout << "xmin " << xmin << std::endl;
    std::cout << "xmax " << xmax << std::endl;
    std::cout << "ymin " << ymin << std::endl;
    std::cout << "ymax " << ymax << std::endl;
    std::cout << "zmin " << zmin << std::endl;
    std::cout << "zmax " << zmax << std::endl;
    std::cout << "xcells " << xcells << std::endl;
    std::cout << "ycells " << ycells << std::endl;
    std::cout << "zcells " << zcells << std::endl;
    std::cout << "temperature " << temperature << std::endl;
    std::cout << "Boltzmann Constant " << boltzmannConstant << std::endl;
    std::cout << "/++++++++++++++++++++Parameters of particles/species++++++++++++++++++++/"     << std::endl;
    std::cout << "Species Input File: "  << speciesInputFile  << std::endl;
    std::cout << "Patch Input File: "  << patchInputFile  << std::endl;
    std::cout << "Total number of species: "   << totalNumberOfSpecies   << std::endl;
    std::cout << "/++++++++++++++++++++Parameters of Lennard-Jones potential ++++++++++++++++++++/"     << std::endl;
    std::cout << "deltaT " << deltaT << std::endl;
    std::cout << "tEnd " << tEnd << std::endl;
    std::cout << "depthOfThePotentialWell " << depthOfThePotentialWell << std::endl;
    std::cout << "sigmaPreFactor " << sigmaPreFactor << std::endl;
    std::cout << "reactionDistancePreFactor " << reactionDistancePreFactor << std::endl;
    std::cout << std::endl;
    std::cout << "/++++++++++++++++++++THE DERIVED PARAMETERS++++++++++++++++++++/"     << std::endl;
    std::cout << "/++++++++++++++++++++Parameters of the world++++++++++++++++++++/"     << std::endl;
    std::cout << "h " << h << std::endl;
    std::cout << "totalNumberOfCells "    << totalNumberOfCells << std::endl;
    std::cout << "/++++++++++++++++++++Parameters of Lennard-Jones potential ++++++++++++++++++++/"     << std::endl;
     std::cout << "epsilon " << epsilon << std::endl;
    std::cout << "sigma " << sigma << std::endl;
    std::cout << "gamma " << gamma << std::endl;
    std::cout << "gammaRot " << gammaRot << std::endl;
    std::cout << "reactionDistance " << reactionDistance << std::endl;
    std::cout << "/++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++/" << std::endl;
    std::cout << std::endl;   
}
