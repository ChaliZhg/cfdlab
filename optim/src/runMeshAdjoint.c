#include <stdio.h>
#include <stdlib.h>
#include "opt.h"

void runMeshAdjoint()
{
   printf("Executing %s\n", adjmesh);
   system(adjmesh);
}
