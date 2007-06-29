#ifndef _OPT_H

#define _OPT_H 1

double cl, cd, clref, cdref, cost0;
double CP0[1000], CP[1000];
int costfun;
char rundir[100];
char deform[100];
char flosol[100];
char adjsol[100];
char adjmesh[100];

double objFun(int, double *);
void objGrad(int, double*, double*);

#endif
