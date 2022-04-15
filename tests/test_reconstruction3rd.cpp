#include "../init.c"
#include "../reconstruction_third.cpp"
#include<stdio.h>

int main(){
  startup();
  boundaries();
  //conv_to_prims(main_grid);

  //if(reconstruction_first_order()) printf("success\n");
  reconstruction_third_order();
  main_grid.test_interfaces(0.1);
	main_grid.print_grid();

  return 0;
}
