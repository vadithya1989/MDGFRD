#include "externDefs.h"
#include "myVector.h"
#include "particle.h"
#include "species.h"

void measureInterPatchDistance( myVector &r, 
                                std::vector< std::pair<int, myVector> > &rPatchVector,
                                particle &pi,
                                particle &pj
                              )
{
    for (int l=0; l<speciesList[pi.speciesType].numberOfPatches;++l)
    {
        for (int m=0; m<speciesList[pj.speciesType].numberOfPatches;++m)
        {
            //vectors from the particle centre to the patch centre
            myVector riPatch;   ///vector from the centre of particle i to the patch centre
            myVector rjPatch;   ///vector from the centre of particle j to the patch centre
            riPatch = pi.q.bodyToSpace(speciesList[pi.speciesType].bodyPatchVectors[l]) * speciesList[pi.speciesType].particleRadius;
            rjPatch = pj.q.bodyToSpace(speciesList[pj.speciesType].bodyPatchVectors[m]) * speciesList[pj.speciesType].particleRadius;
            int patchID=0;
            if ( (pi.speciesType == 0) || (pi.speciesType == 1) || (pi.speciesType == 2) )
            {
                patchID = l;
            }
            else if( (pj.speciesType == 0) || (pj.speciesType == 1) || (pj.speciesType == 2) )
            {
                patchID = m;
            }
            
            std::pair<int, myVector> patchVectorPair(patchID, r + (rjPatch - riPatch));
            rPatchVector.push_back(patchVectorPair);
        }
    }
}
