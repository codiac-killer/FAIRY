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
		for (int i = 0; i < ndim; ++i)
		{
			u[i]=0;
		}
	}
};

// To return  type std::ostream&: to onoma tis sinartisis einai praktika ta velakia <<, otan kalis ta velakia auta exoun san string oti exei proigithi prin ta velakia 
// san os<<pr , kai an dn exei proigithi tipota praktika einai "keno". Kai epistrefi to idio "string" enw tou exei prosthesi ta variables tou primitive 
std::ostream& operator<<(std::ostream& os, const primitive& pr){
    os << pr.p_th << ',' << pr.rho;
    return os;
}


// kanoume klassi gia to grid
class ndvector {
public: 
	// ta int einai oi diastasis
	int i,j,k;
	std::vector<primitive> array;
	// apaititai mono i 1i diastasi, ta alla mpainoun 1
	ndvector(int ii, int jj=1, int kk=1){
		i = ii;
		j = jj;
		k = kk;
		array.resize(i*j*k);
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


int main(int argc, char const *argv[])
{
	//  (a) ? b : c if a -> b else -> c
	ndvector main_grid(idim, (ndim>1) ? jdim : 1, (ndim>2) ? kdim : 1);


	// na allaksoume to periexomeno tou eswrterikou , alla mono ta misa stixia (wste na tiposi misa alagmena misa default)
	for (int i = 0; i < int(main_grid.array.size()/2); ++i)
	{
		main_grid.array[i].p_th = 0.5;
	}


	main_grid.print();

	return 0;
}
