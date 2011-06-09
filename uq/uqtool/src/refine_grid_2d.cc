/*
 *  refine_grid_2d.cc
 *  uqtool
 *
 *  Created by Praveen Chandrashekar on 07/06/11.
 *  Copyright 2011 TIFR-CAM, Bangalore. All rights reserved.
 *
 */


#include <iostream>
#include <cassert>
#include "uq.h"

using namespace std;

// Refine the stochastic grid in 2-D
// * = new samples
// Edge 0-1 is the largest edge and it gets divided
//
// Linear element:
//                     2                              2
//                    / \                            /|\
//                   /   \                          / | \
// old element      /     \       new element      /  |  \
//                 /       \                      /   |   \
//                /         \                    /    |    \
//               0-----------1                  0-----*-----1
//
//
// Quadratic element:
//                     2
//                    / \
//                   /   \
//                  5     4
//                 /       \
//                /         \
// old element:  0-----3-----1
//
//                     2
//                    /|\
//                   / | \
//                  5  *  4
//                 /   |   \
//                /    |    \
// old element:  0--*--3--*--1
//
template <>
void UQProblem<2>::refine_grid (bool eno_mode)
{
   unsigned int n_element = grid.element.size();
   
   for(unsigned int i=0; i<n_element; ++i)
      if(grid.element[i].refine_flag)
      {         
         // Refine this element into two
         cout << "Refining element = " << i;
         
         // If in eno mode, we need quadratic element
         if(eno_mode)
            assert (grid.element[i].order == 2);
         
         unsigned int i0 = grid.element[i].idof[0];
         unsigned int i1 = grid.element[i].idof[1];
         unsigned int i2 = grid.element[i].idof[2];
         
         valarray<double> x0 = sample[i0].x;
         valarray<double> x1 = sample[i1].x;
         valarray<double> x2 = sample[i2].x;
         
         vector<unsigned int> e = grid.element[i].largest_face (x0, x1, x2);
         
         if(grid.element[i].order == 1 && order == 1)
         {
            // Divide linear element into two linear elements
            cout << ": dividing P1 element into two P1 elements\n";
            
            // sample at middle of largest edge
            valarray<double> x = (sample[e[0]].x + sample[e[1]].x) / 2.0;
            
            // check if this sample already exists
            int s = find_sample (x);
            
            if (s == -1) // sample does not exist
            {
               // One new sample
               Sample<2> new_sample1 (n_var, n_cell, n_moment, sample.size());
               new_sample1.x = x;
               sample.push_back (new_sample1);
               s = sample.size() - 1;
            }
            
            // Two new linear elements
            Element<2> new_element1 (1, n_moment, n_cell, grid.element.size());
            new_element1.idof[0] = e[0];
            new_element1.idof[1] = s;
            new_element1.idof[2] = e[2];
            grid.element.push_back (new_element1);
            
            Element<2> new_element2 (1, n_moment, n_cell, grid.element.size());
            new_element2.idof[0] = s;
            new_element2.idof[1] = e[1];
            new_element2.idof[2] = e[2];
            grid.element.push_back (new_element2);
         }
         else if(grid.element[i].order == 1 && order == 2)
         {
            // Convert linear element into quadratic element
            cout << ": converting P1 element into P2 element\n";
            abort ();

         }
         else if(grid.element[i].order == 2 && eno_mode == false)
         {
            // Divide quadratic element into two quadratic elements
            cout << ": dividing P2 element into two P2 elements\n";
            
//            cout << e[0] << "  " << sample[e[0]].x[0] << "  " <<
//            sample[e[0]].x[1] << endl;
//            cout << e[1] << "  " << sample[e[1]].x[0] << "  " <<
//            sample[e[1]].x[1] << endl;
            // One new sample at middle of median
            Sample<2> new_sample0 (n_var, n_cell, n_moment, sample.size());
            new_sample0.x = (sample[e[3]].x + sample[e[2]].x) / 2.0;
            sample.push_back (new_sample0);
            int s0 = sample.size() - 1;
            
            valarray<double> x1 = (sample[e[0]].x + sample[e[3]].x) / 2.0;
            int s1 = find_sample (x1);

            if (s1 == -1) // sample does not exist
            {
               // One new sample
               Sample<2> new_sample1 (n_var, n_cell, n_moment, sample.size());
               new_sample1.x = x1;
               sample.push_back (new_sample1);
               s1 = sample.size() - 1;
            }
            
            valarray<double> x2 = (sample[e[1]].x + sample[e[3]].x) / 2.0;
            int s2 = find_sample (x2);
       
            if (s2 == -1) // sample does not exist
            {
               // One new sample
               Sample<2> new_sample2 (n_var, n_cell, n_moment, sample.size());
               new_sample2.x = x2;
               sample.push_back (new_sample2);
               s2 = sample.size() - 1;
            }
            
            // Two new quadratic elements
            Element<2> new_element1 (2, n_moment, n_cell, grid.element.size());
            new_element1.idof[0] = e[0];
            new_element1.idof[1] = e[3];
            new_element1.idof[2] = e[2];
            new_element1.idof[3] = s1;
            new_element1.idof[4] = s0;
            new_element1.idof[5] = e[5];
            grid.element.push_back (new_element1);
            
            Element<2> new_element2 (2, n_moment, n_cell, grid.element.size());
            new_element2.idof[0] = e[3];
            new_element2.idof[1] = e[1];
            new_element2.idof[2] = e[2];
            new_element2.idof[3] = s2;
            new_element2.idof[4] = e[4];
            new_element2.idof[5] = s0;
            grid.element.push_back (new_element2);
         }
         else if(grid.element[i].order == 2 && eno_mode == true)
         {
            // Divide quadratic element into two linear elements
            cout << ": dividing P2 element into two P1 elements\n";
            abort ();
         }
         else
         {
            cout << "Cannot refine grid for order =" 
                 << grid.element[i].order
                 << endl;
            abort ();
         }
         
         grid.element[i].active      = false;
         grid.element[i].refine_flag = false;
      }   
}
