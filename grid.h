#include <iostream>
#include <cmath>
#include <math.h>  
#include <vector>
#include <stdio.h>
#include "init_params.h"
#include <functional>

static const int        i_min = 1, i_max = 2;
static const int        j_min = 1, j_max = 2;
static const int        k_min = 1, k_max = 2;
static const char*		space = "linear";
static const char*    bound_11 = "outflow";
static const char*    bound_12 = "outflow";

double d_i = ((float)i_max - (float)i_min)/(float)idim ;//  dX of numerical grid

double x[2*idim+1+4*dim_b];
double z[0]; 
double y[0];


class primitive {
public:
// orizoume to u pinaka me ta primitves 
	double p_th, rho;
	double u[n_comp];
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


class conservable {
public:
// orizoume to u pinaka me ta primitves 
	float rho_, e_;
	float p_[n_comp];
	
	conservable(){
		// default times se periptosi pou dn dosei o user 
		rho_ = 1;
		e_  = 10;
		for (int i = 0; i < n_comp; ++i)
		{
			p_[i]=0;
		}
	}
	int build(primitive input_){
		double w_;
		double lor =0;
		//First add the u**2 part of the lorentz factor.
		for (int i = 0; i < n_comp; ++i)
		{
			lor = lor + input_.u[i]*input_.u[i];
		}
		//now we make the correct lorentz factor
		lor = pow(1/(1-lor),0.5);
		
		//definition of w
		w_= input_.rho+gamma_par*input_.p_th;


		// Now values of conserved ( eq. (5) of delzana_HD )
		rho_ = input_.rho*lor;
		e_	 = w_*lor*lor-input_.p_th;


		for (int i = 0; i < n_comp; ++i)
		{
			p_[i]= w_*lor*lor*input_.u[i];
		}


		return 0;
	}
};



class flux {
public:
// orizoume to u pinaka me ta primitves 
	float rhof_[n_comp], ef_[n_comp];
	float pf_[n_comp][n_comp];
	
	flux(){
		// default times se periptosi pou dn dosei o user 

		for (int i = 0; i < n_comp; ++i)
		{
			rhof_[i] = 1;
			ef_[i]  = 1;		
			pf_[i][i]=1;
		}
	}
	int build(primitive input_){
		double w_;
		double lor =0;
		//First add the u**2 part of the lorentz factor.
		for (int i = 0; i < n_comp; ++i)
		{
			lor = lor + input_.u[i]*input_.u[i];
		}
		//now we make the correct lorentz factor
		lor = pow(1/(1-lor),0.5);
		
		//definition of w
		w_= input_.rho+gamma_par*input_.p_th;





		for (int i = 0; i < n_comp; ++i)
			{
			// Now values of fluxes ( eq. (6) of delzana_HD )
			rhof_[i] = input_.rho*lor*input_.u[i];
			ef_[i]	 = w_*lor*lor*input_.u[i];

			for (int j = 0; j < n_comp; ++j)
			
			
				{

					pf_[i][j]= w_*lor*lor*input_.u[i]*input_.u[j]+input_.p_th*(i==j);
				}
			}

		return 0;
	}
};



std::ostream& operator<<(std::ostream& os, const primitive& pr){
    os << pr.p_th << ',' << pr.rho << ',' << pr.u[0] << ',' << pr.u[1] << ',' << pr.u[2] << ','<< (1/(1-pow(pr.u[0],2)-pow(pr.u[1],2)-pow(pr.u[2],2)));
    return os;
}


// kanoume klassi gia to grid
class ndvector {
public: 
	// ta int einai oi diastasis
	int i,j,k;
	std::vector<primitive> array; // TODO : array se center
	std::vector<primitive> interfaces;
	std::vector<conservable> intf_cons;
	std::vector<conservable> conservables_;
	std::vector<flux> fluxes_;
	std::vector<flux> intf_flx;
	std::vector<flux> solv_flx;
	std::vector<double> eigens;
	// apaititai mono i 1i diastasi, ta alla mpainoun 1
	ndvector(int ii, int jj=1, int kk=1){
		i = ii;
		j = jj;
		k = kk;
		array.resize(2*(i+(i>1)*(2*dim_b+1))*(j+(j>1)*(2*dim_b+1))*(k+(k>1)*(2*dim_b+1)));  // 2 elements (center & left interface) for each cell plus one extra if for each non ekfulismeni dimension
		//array.resize(i*i*j*k);
		interfaces.resize(2*(i+(i>1)*2*dim_b)*j*k);
		intf_cons.resize(2*(i+(i>1)*2*dim_b)*j*k);
		intf_flx.resize(2*(i+(i>1)*2*dim_b)*j*k);
		solv_flx.resize((i+(i>1)*2*dim_b)*j*k);
		conservables_.resize(2*(i+(i>1)*2*dim_b)*j*k);
		fluxes_.resize(2*(i+(i>1)*2*dim_b)*j*k);
		eigens.resize(2*(i+(i>1)*2*dim_b)*j*k);
	}
	// me tin print tiponoume ta incules tou primitive
	void print_cells(){
		for (int Di=0; Di<i; Di++){
			for (int Dj=0; Dj<j; Dj++){
				for (int Dk=0; Dk<k; Dk++){
					std::cout<<array[Di*j*k+Dj*k+Dk]<<" ";
				}
				std::cout<<"\n";
			}
			std::cout<<"\n";
		}
	}
	void print_grid(){
		for (int Di=0; Di<2*idim+1+4*dim_b; Di++){
			printf("%d %f \n", Di, x[Di]);
		}
		printf("\n");
	}

	void print_interfaces(){
    printf("|");
		for (int Di=0; Di<i+2*dim_b; Di++){
      	    printf("%f %f \n", interfaces[2*Di].p_th, interfaces[2*Di+1].p_th);
    }
    printf("\n");
	}

	void test_interfaces(int error_margin){
		int Di;
		for (Di=0; Di<i+1; Di++){
      	    printf("%f %f %f\n", interfaces[2*(Di+dim_b)-1].p_th, interfaces[2*(Di+dim_b)].p_th, x[2*(Di+dim_b)]);  // sin(2*M_PI*(x[2*Di+1+dim_b]-1)));
    }
    printf("\n");
	}
};


// dim_b 
ndvector main_grid(idim, (ndim>1) ? jdim : 1, (ndim>2) ? kdim : 1);
