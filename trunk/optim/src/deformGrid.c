#include<stdio.h>
#include<stdlib.h>

void deformGrid()
{
   extern char deform[];
   printf("Executing %s\n", deform);
   system(deform);
}
