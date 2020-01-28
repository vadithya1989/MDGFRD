
#include <iosfwd>
#include <iostream>
#include <cmath>
#include <sstream>
#include <fstream>

#include "particle.h"
#include "simulationParameters.h"

void writeOutput ( 
                    int &vtkCount)
{
    std::stringstream s;
    std::string file_VTK;

    s << "particle" << vtkCount << ".vtk";
    s >> file_VTK;
    std::cout << tSim << "  .....saving VTK number " << vtkCount << std::endl;
    
    std::ofstream vtkFile(file_VTK.c_str(), std::ifstream::out);
    if ( vtkFile.good() )
    {
        vtkFile.precision(15);

        /*Header*/
        vtkFile << "# vtk DataFile Version 3.0" << std::endl;
        vtkFile << "MD particles" << std::endl;
        vtkFile << "ASCII" << std::endl;
        vtkFile << "DATASET UNSTRUCTURED_GRID" << std::endl;
        
        /*write the grid from particlePositionitions of particles*/
        vtkFile << "POINTS " << size_t(particleList.size()) << " float" << std::endl;
        for (int i = 0; i < particleList.size(); ++i )
        {
            double delx = particleList[i].particlePosition.x;
            double dely = particleList[i].particlePosition.y;
            double delz = particleList[i].particlePosition.z;
            myVector rPb;
            rPb.x = particleList[i].particlePosition.x - ( floor( double(delx)/double(prm.xmax) ) * prm.xmax );
            rPb.y = particleList[i].particlePosition.y - ( floor( double(dely)/double(prm.ymax) ) * prm.ymax );
            rPb.z = particleList[i].particlePosition.z - ( floor( double(delz)/double(prm.zmax) ) * prm.zmax );

            vtkFile << std::scientific << "  " <<rPb.x << "  " << rPb.y << "  " << rPb.z << std::endl;
        }
        /*header for the remaining scalars and vectors*/
        vtkFile << "POINT_DATA " << size_t(particleList.size()) << std::endl;

        /*write radius*/
        vtkFile << "SCALARS radius float 1"  << std::endl;
        vtkFile << "LOOKUP_TABLE default"  << std::endl;
        for (int i = 0; i < particleList.size(); ++i )
        {
            vtkFile << "  " << speciesList[particleList[i].speciesType].particleRadius << std::endl;
        }
        vtkFile << "SCALARS speciesType float 1"  << std::endl;
        vtkFile << "LOOKUP_TABLE default"  << std::endl;
        for (int i = 0; i < particleList.size(); ++i )
        {
            vtkFile << "  " << particleList[i].speciesType << std::endl;
        }
        
        /*write velocity*/
       /* vtkFile << "VECTORS velocity float" << std::endl;
        for (int i = 0; i < particleList.size(); ++i )
        {
            vtkFile << std::scientific << " " << particleList[i].momentum.x << "  " << particleList[i].momentum.y << "  " << particleList[i].momentum.z << std::endl;
        }*/
        
        /*write force*/
        vtkFile << "VECTORS force float" << std::endl;
        for (int i = 0; i < particleList.size(); ++i )
        {
            vtkFile << std::scientific << "  " << particleList[i].force.x << "  " << particleList[i].force.y << "  " << particleList[i].force.z << std::endl;
        }
    }
    
    std::stringstream s1;
    std::string file_patch_VTK;

    s1 << "patch_" << vtkCount << ".vtk";
    s1 >> file_patch_VTK;
   
    std::ofstream vtkPatchFile(file_patch_VTK.c_str(), std::ifstream::out);
      if(vtkPatchFile.good())
      {
        vtkPatchFile.precision(15);

        //Header
        vtkPatchFile << "# vtk DataFile Version 3.0" << std::endl;
        vtkPatchFile << "MD patch particles" << std::endl;
        vtkPatchFile << "ASCII" << std::endl;
        vtkPatchFile << "DATASET UNSTRUCTURED_GRID" << std::endl;

        //write the grid from positions of particles
        vtkPatchFile << "POINTS " << size_t( particleList.size()) << " float" << std::endl;
        for (int i = 0; i <  particleList.size(); ++i )
        {
            for (int j=0; j<speciesList[particleList[i].speciesType].numberOfPatches; ++j)
            {
                
                vtkPatchFile << std::fixed << "  " << particleList[i].patchPos[j].x << "  " << particleList[i].patchPos[j].y << "  " << particleList[i].patchPos[j].z << std::endl;
            }
        }
        //header for the remaining scalars and vectors
        vtkPatchFile << "POINT_DATA " << size_t( particleList.size()) << std::endl;

        // write radius
        vtkPatchFile << "SCALARS radius float 1"  << std::endl;
        vtkPatchFile << "LOOKUP_TABLE default"  << std::endl;
        for (int i = 0; i <  particleList.size(); ++i )
        {
            for (int j=0; j<speciesList[particleList[i].speciesType].numberOfPatches; j++)
            {
                vtkPatchFile << "  " << 0.3*speciesList[particleList[i].speciesType].particleRadius << std::endl;
            }
        }

      }
      
      //NOW FOR THE DOMAIN
    std::stringstream s2;
    std::string fileD_VTK;

    s2 << "domain" << vtkCount << ".vtk";
    s2 >> fileD_VTK;

    std::ofstream vtkDFile(fileD_VTK.c_str(), std::ifstream::out);
    if ( vtkDFile.good() )
    {
        vtkDFile.precision(15);

        //Header
        vtkDFile << "# vtk DataFile Version 3.0" << std::endl;
        vtkDFile << "eGFRD domains" << std::endl;
        vtkDFile << "ASCII" << std::endl;
        vtkDFile << "DATASET UNSTRUCTURED_GRID" << std::endl;

        //write the grid from positions of particles
        vtkDFile << "POINTS " << size_t(particleList.size()) << " float" << std::endl;
        for (int i = 0; i < particleList.size(); ++i )
        {
            double delx = particleList[i].particlePosition.x;
            double dely = particleList[i].particlePosition.y;
            double delz = particleList[i].particlePosition.z;
            myVector rPb;
            rPb.x = particleList[i].particlePosition.x - ( floor( double(delx)/double(prm.xmax) ) * prm.xmax );
            rPb.y = particleList[i].particlePosition.y - ( floor( double(dely)/double(prm.ymax) ) * prm.ymax );
            rPb.z = particleList[i].particlePosition.z - ( floor( double(delz)/double(prm.zmax) ) * prm.zmax );

            vtkDFile << std::scientific << "  " <<rPb.x << "  " << rPb.y << "  " << rPb.z << std::endl;
        }
        //header for the remaining scalars and vectors
        vtkDFile << "POINT_DATA " << size_t(particleList.size()) << std::endl;

        //write radius of the domain
        vtkDFile << "SCALARS radius float 1"  << std::endl;
        vtkDFile << "LOOKUP_TABLE default"  << std::endl;
        for (int i = 0; i < particleList.size(); ++i )
        {
            if(particleList[i].particleType == 1)
                vtkDFile << "  " << particleList[i].domainRadius << std::endl;
            else
                vtkDFile << "  " << 1e-10 << std::endl;
        }
    }

      
      
}

