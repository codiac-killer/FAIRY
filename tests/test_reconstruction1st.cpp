#include "../init.c"
#include "../reconstruction_first.cpp"
#include<stdio.h>

int main(){
  startup();
  //conv_to_prims(main_grid);

  if(reconstruction_first_order()) printf("success\n");

  main_grid.print_interfaces();

  return 0;
}
