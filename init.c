#include <iostream>
#include <cmath>
#include <vector>
#include <stdio.h>
#include <functional>
#include "grid.h"


static const int        i_min = 1, i_max = 2;
static const int        j_min = 1, j_max = 2;
static const int        k_min = 1, k_max = 2;
static const char*		space = "linear";
static const char*    bound_11 = "outflow";
static const char*    bound_12 = "outflow";

double d_i = ((float)i_max - (float)i_min)/(float)idim ;//  dX of numerical grid

double x[2*idim+1+2*dim_b]; // +2 due to ghost zones , +1 due to i/2
// double x[2*idim+1]; // trying without dim_b cuz its not needed in Rieman Solver. If needed, than will revert back
double z[0]; 
double y[0];


int startup() {

  if (space == "linear"){
      if (ndim == 1){
        for (int i = 0; i < 2*idim+4*dim_b+1; ++i) {// + 2* cells to check from ghost zones
        // for (int i = 0; i < 2*idim+1; ++i) { // same as L16
        x[i] = (i_min - d_i*dim_b)+i*d_i/2;
        //x[i] = i_min +i*d_i/2; // same as L16
        y[0] = j_min;
        z[0] = k_min;
        }
      }
  } 

  int p;
  // Scanning the grid for Initial Conditions. main_grid the ndvector from "grid.h"
  for (int i = dim_b; i < idim+dim_b; ++i) {
    for (int j = 0; j < main_grid.j; ++j) {
      for (int k = 0; k < main_grid.k; ++k) {
        p = i*main_grid.j*main_grid.k+j*main_grid.k+k; // fix the pointer to the "array"
        
        // x,y,z from grid.h used as reference to the real numerical distance rather than cell ID
        //x[2*p+1] an theloume sinartisi tis thesis
        main_grid.array[p].p_th = sin(2*M_PI*(x[2*p+1]-1)); 
        main_grid.array[p].rho  = 1;
        //printf("%d %f \n", p,main_grid.array[p].p_th);

        for (int ii = 0; ii < n_comp; ++ii) {
          main_grid.array[p].u[ii]  = 0.5; 
        }
        // builed conservables
        main_grid.conservables_[p].build(main_grid.array[p]);
      }
    }
  }

  return 0; }

int boundaries() {

int p,p1;

  if (bound_11 == "outflow"){

  for (int i = 0; i < dim_b; ++i) {
    for (int j = 0; j < main_grid.j; ++j) {
      for (int k = 0; k < main_grid.k; ++k) {
        p = i*main_grid.j*main_grid.k+j*main_grid.k+k; // fix the pointer to the "array"
        p1= dim_b*main_grid.j*main_grid.k+j*main_grid.k+k; // fix the pointer to the first non-boundary cell

        // x,y,z from grid.h used as reference to the real numerical distance rather than cell ID
        //x[2*p+1] an theloume sinartisi tis thesis
        main_grid.array[p].p_th = main_grid.array[p1].p_th  ; 
        main_grid.array[p].rho  = main_grid.array[p1].rho;

        for (int ii = 0; ii < n_comp; ++ii) {
          main_grid.array[p].u[ii]  = main_grid.array[p1].u[ii]; 
        }
        // builed conservables
        // main_grid.conservables_[p].build(main_grid.array[p]);
      }
    }
  }
}

  if (bound_12 == "outflow"){

  for (int i = idim+dim_b; i < idim+2*dim_b; ++i) {
    for (int j = 0; j < main_grid.j; ++j) {
      for (int k = 0; k < main_grid.k; ++k) {
        p = i*main_grid.j*main_grid.k+j*main_grid.k+k; // fix the pointer to the "array"
        p1= (idim+dim_b-1)*main_grid.j*main_grid.k+j*main_grid.k+k; // fix the pointer to the first non-boundary cell

        // x,y,z from grid.h used as reference to the real numerical distance rather than cell ID
        //x[2*p+1] an theloume sinartisi tis thesis

        main_grid.array[p].p_th = main_grid.array[p1].p_th  ; 
        main_grid.array[p].rho  = main_grid.array[p1].rho;


        for (int ii = 0; ii < n_comp; ++ii) {
          main_grid.array[p].u[ii]  = main_grid.array[p1].u[ii]; 
        }
        // builed conservables
        main_grid.conservables_[p].build(main_grid.array[p]);
      }
    }
  }
}





  return 0;
}
