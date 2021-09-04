#include <iostream>
#include <cmath>
#include <vector>
#include <stdio.h>
#include "init_params.h"
#include "basic_io.h"
#include <functional>






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

std::ostream& operator<<(std::ostream& os, const primitive& pr){
    os << pr.p_th << ',' << pr.rho << ',' << pr.u[0] << ',' << pr.u[1] << ',' << pr.u[2] << ','<< (1/(1-pow(pr.u[0],2)-pow(pr.u[1],2)-pow(pr.u[2],2)));
    return os;
}



// kanoume klassi gia to grid
class ndvector {
public: 
	// ta int einai oi diastasis
	int i,j,k;
	std::vector<primitive> array;
	std::vector<conservable> conservables_;
	// apaititai mono i 1i diastasi, ta alla mpainoun 1
	ndvector(int ii, int jj=1, int kk=1){
		i = ii;
		j = jj;
		k = kk;
		array.resize(i*j*k);
		conservables_.resize(i*j*k);
	}
	// me tin print tiponoume ta incules tou primitive
	void print(){
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
};




// dim_b 
ndvector main_grid(idim+2*dim_b, (ndim>1) ? jdim : 1, (ndim>2) ? kdim : 1);