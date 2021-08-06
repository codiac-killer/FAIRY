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


int main(int argc, char const *argv[])
{


double x[2*idim+1+2*dim_b]; // +2 due to ghost zones , +1 due to i/2
double z[0]; 
double y[0];
if (space == "linear"){
		if (ndim == 1){


			double d_i = (i_max - i_min)/idim ;//  dX of numerical grid
			for (int i = 0; i < 2*idim+dim_b+1; ++i) // + 2 cells to check from ghost zones
			{
				x[i] = (i_min - d_i/2)+i*d_i/2;
			}
			y[0] = j_min;
			z[0] = k_min;
			


		}

} 


int p;
// Scanning the grid for Initial Conditions. main_grid the ndvector from "grid.h"
for (int i = 2; i < main_grid.i; ++i)
{
	for (int j = 0; j < main_grid.j; ++j)
	{
		for (int k = 0; k < main_grid.k; ++k)
		{
			//p = i*main_grid.j*main_grid.k+j*main_grid.k+k; // fix the pointer to the "array"

			// x,y,z from grid.h used as reference to the real numerical distance rather than cell ID
			//main_grid.array[p].p_th = x[i]; 
			//main_grid.array[p].rho  = y[j];

			for (int ii = 0; ii < ndim; ++ii)
			{
				//main_grid.array[p].u[ii]  = x[i]/(x[i]+z[k]);
			}
		}
	}
}
	return 0;
}