//
//  phosphorylationReaction.cpp
//  mapK
//
//  Created by Adithya on 08/06/2017.
//  Copyright © 2017 Adithya. All rights reserved.
//

#include "externDefs.h"
#include "simulationParameters.h"
#include "particle.h"

void phosphorylationReaction()
{
    ///if simulation time is larger than first next event time in scheduler list
    if (particleList[schedulerList[0]].particleType == 2     &&
        tSim >= particleList[schedulerList[0]].nextEventTime &&
        particleList[schedulerList[0]].nextEventType == 1)
    {
        
    }
}
