#include "../init.c"
#include "../reconstruction_third.cpp"
#include<stdio.h>

int main(){
  startup();
  boundaries();
  //conv_to_prims(main_grid);

  //if(reconstruction_first_order()) printf("success\n");
  reconstruction_third_order();
  main_grid.print_interfaces();
  /*for (int i = 0; i <  main_grid.array.size(); ++i)
  {
  	std::cout<<main_grid.array[i].p_th;
  }*/
 

  return 0;
}
