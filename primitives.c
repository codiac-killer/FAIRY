#include <iostream>
#include <cmath>
#include <vector>
#include <stdio.h>
#include <functional>
#include "init.c"


double F(double Lor, double q2, double e, double d){
	return q2/(1.-1./Lor/Lor)-pow((e*gamma_par*Lor*Lor-d*Lor)/(gamma_par*Lor*Lor-1.0),2);
}


double Fd(double Lor, double q2, double e, double d){
	return -2*q2/pow((1.-1./Lor/Lor),2)*pow(Lor,-3)-2*((e*gamma_par*Lor*Lor-d*Lor)/(gamma_par*Lor*Lor-1.0))*((2*e*gamma_par*Lor-d)*(gamma_par*Lor*Lor-1.0)-(e*gamma_par*Lor*Lor-d*Lor)*(2*gamma_par*Lor))/pow((gamma_par*Lor*Lor-1.0),2);
}

// Quantities for Solution (8,9,10 delzana)
double W,E,D;
double Q2;
//error for the lorentz factor 
double error = 0.01;
// how many itteration for N-R numerical solution
int n = 1; 
int MAXITER = 100;
double x_lor;
int p;

// conversion from conservatives to primitves needed for each step
int conv_to_prims(ndvector main_grid){
for (int i = 2; i < main_grid.i; ++i){
	for (int j = 0; j < main_grid.j; ++j)
	{
		for (int k = 0; k < main_grid.k; ++k)
		{


			p = i*main_grid.j*main_grid.k+j*main_grid.k+k; // fix the pointer to the "array"







			// Build W,D,Q,E from (8,9,10) Del_Zana
			E = main_grid.conservables_[p].e_;
			D = main_grid.conservables_[p].rho_;
			Q2=0;
						for (int ii = 0; ii < n_comp; ++ii)
			{
				Q2 = Q2+pow(main_grid.conservables_[p].p_[ii],2);
			}




			
			n = 1;
			x_lor = 1.2;

			while( ( fabs(F(x_lor, Q2, E, D)) > error ) && ( n <= MAXITER ) )
			{
			    x_lor = x_lor - ( F(x_lor, Q2, E, D)/ Fd(x_lor, Q2, E, D) );

			    n++;
			}
			std::cout << x_lor << ',' <<'\n';

			
			





		}
	}
}





	return 0;
}

int main()
{

	startup();
	conv_to_prims(main_grid);
	for (int i = 0; i < idim; ++i)
	{
	//std::cout << main_grid.array[i] << ' '  << '\n';
	}
	
	return 0;

}


/*

			// x,y,z from grid.h used as reference to the real numerical distance rather than cell ID
			main_grid.array[p].p_th = x[i]; 
			main_grid.array[p].rho  = y[j];

			for (int ii = 0; ii < ndim; ++ii)
			{
				main_grid.array[p].u[ii]  = x[i]/(x[i]+z[k]);
			}
			// builed conservable
			main_grid.conservables_[p].build(main_grid.array[p]);
*/