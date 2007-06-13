#include<stdio.h>
#include<stdlib.h>

void runSolver()
{
   extern char flosol[];
   printf("Executing %s\n", flosol);
   system(flosol);
}
