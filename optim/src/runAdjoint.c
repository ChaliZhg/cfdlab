#include <stdio.h>
#include <stdlib.h>
#include "opt.h"

void runAdjoint()
{
   printf("Executing %s\n", adjsol);
   system(adjsol);
}
