#ifndef __CONSTANTS_H__
#define __CONSTANTS_H__

// Coefficients for 3-stage RK scheme of Shu-Osher
static const double a_rk[] = {0.0, 3.0/4.0, 1.0/3.0};
static const double b_rk[] = {1.0, 1.0/4.0, 2.0/3.0};

static const double KKK = 1.0/3.0;

enum GridType {gmsh, bamg, delaundo};
enum CellType {median, voronoi};
enum Dimension {two, axi};

#endif
