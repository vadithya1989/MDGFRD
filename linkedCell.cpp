
#include "simulationParameters.h"
#include "externDefs.h"
#include "main.h"
#include <vector>

/**
 Function resets all the cells for a new linked cell sorting
 */
void blankCells ()
{
    cells.clear();
    for(int i=0; i<prm.totalNumberOfCells; ++i)
    {
        cells.push_back(-1);
    }
}

/**
 Function determines the cell to which the particle belongs in

 @param Particle particle whose cell is to be determined
 @return cell number
 */
int findCellOfTheParticle(
                          particle &Particle)
{
    myVector rPb;
    double delx = Particle.particlePosition.x;
    double dely = Particle.particlePosition.y;
    double delz = Particle.particlePosition.z;
    
    rPb.x = Particle.particlePosition.x - (floor(double(delx)/double(prm.xmax)) * prm.xmax);
    rPb.y = Particle.particlePosition.y - (floor(double(dely)/double(prm.ymax)) * prm.ymax);
    rPb.z = Particle.particlePosition.z - (floor(double(delz)/double(prm.zmax)) * prm.zmax);
    
    int coord = rPb.z/prm.h;
    int cellID = coord * prm.ycells;

    coord = rPb.y/prm.h;
    cellID += coord;
    cellID *= prm.xcells;

    coord = rPb.x/prm.h;
    cellID += coord;
    return cellID;
}

/**
 Sorts the particles into cells

 */
void sortParticles ()
{
    particles.clear();
    for (uint i=0; i<particleList.size(); ++i)
    {
        int coord = findCellOfTheParticle ( particleList[i]);
        particles.push_back(cells[coord]);
        cells[coord] = i;
    }
}
