#include <iostream>
#include <cmath>
#include <vector>
#include <stdio.h>
#include "init_params.h"
#include "basic_io.h"
#include <functional>

// point value to numerical distance
switch(space){
	case "linear":
		switch(ndim){
			case 1:
			double real_position[2*idim+1+dim_b]; // +2 due to ghost zones , +1 due to i/2
			double d_i = (i_max - i_min)/idim ;//  dX of numerical grid
			for (int i = 0; i < 2*idim+4; ++i) // + 2 cells to check from ghost zones
			{
				real_position[i] = (i_min - d_i/2)+i*d_i/2;
			}
			break;


		}
	break;
}


class primitive {
public:
// orizoume to u pinaka me ta primitves 
	float p_th, rho;
	float u[n_comp];
	primitive(){
		// default times se periptosi pou dn dosei o user 
		p_th = 1;
		rho  = 10;
		for (int i = 0; i < n_comp; ++i)
		{
			u[i]=0;
		}
	}
};



