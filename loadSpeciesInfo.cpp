//#include <fstream>
//#include <sstream>
//#include <iostream>
//#include <cstdlib>
//
//#include "simulationParameters.h"
#include "main.h"
void loadSpeciesInfo ()
{
    /**************************load the species.dat file*/
    std::ifstream speciesInputFile(prm.speciesInputFile.c_str(), std::ifstream::in);
    if(speciesInputFile.good())
    {
        std::string line;
        if(speciesInputFile)
        {
            for(int i=0; i<prm.totalNumberOfSpecies; ++i)
            {
                getline(speciesInputFile, line);
                std::stringstream s1(line);
                species temp;
                s1  >> temp.numberOfParticles 
                    >> temp.numberOfPatches
                    >> temp.diffusionConstant 
                    >> temp.rotationalDiffusionConstant
                    >> temp.particleRadius 
                    >> temp.particleMass >> temp.dissociationConstant 
                    >> temp.massMomOfInertia;
                speciesList.push_back(temp);
            }
        }
        else
        {
            std::cerr << "Error in species file" << std::endl;
        }
    }
    else
    {
        std::cerr << "Error in loading species.dat file" << std::endl;
        exit(-1);
    }
    
    /************************************load the patch.dat file**********/
    std::ifstream patchInputFile(prm.patchInputFile.c_str(), std::ifstream::in);
    
    if(patchInputFile.good())
    {
        std::string line1;
        for(int l=0; l<prm.totalNumberOfSpecies; ++l)
        {
            //speciesList[l].bodyPatchVectors.resize(speciesList[l].numberOfPatches);
            for(int m=0; m<speciesList[l].numberOfPatches; ++m)
            {
                //std::cout << "hi" ;
                getline(patchInputFile, line1);
                myVector t;
                int patchstate;
                std::stringstream s2(line1);
                s2  >> t.x
                    >> t.y
                    >> t.z
                    >> patchstate;
                speciesList[l].bodyPatchVectors.push_back(t);
            }
            //std::cout << std::endl;
        }
       
        
    }
    else
    {
        std::cerr << "Error in loading patch.dat file" << std::endl;
        exit(-1);
    }
}

