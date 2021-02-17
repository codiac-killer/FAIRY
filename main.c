#include <iostream>
#include <cmath>
#include <vector>
#include <stdio.h>
#include "init_params.h"
#include "basic_io.h"

class primitive {
public:
	float p_th, rho;
	primitive(){
		p_th = 1;
		rho  = 10;
	}

};

std::ostream& operator<<(std::ostream& os, const primitive& pr){
    os << pr.p_th << ',' << pr.rho;
    return os;
}



class ndvector {
public: 
	int i,j,k;
	std::vector<primitive> array;
	ndvector(int ii, int jj=1, int kk=1){
		i = ii;
		j = jj;
		k = kk;
		array.resize(i*j*k);
	}
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
	//test prin : std::cout<<idim<<jdim<<"\n"<<std::endl;
    //print(idim,jdim,"\n");
	

	//  (a) ? b : c if a -> b else -> c
	ndvector main_grid(idim, (ndim>1) ? jdim : 1, (ndim>2) ? kdim : 1);

	for (auto cell : main_grid.array){
		cell.p_th = 1;
		cell.rho  = 2;
	}
	main_grid.print();

	return 0;
}