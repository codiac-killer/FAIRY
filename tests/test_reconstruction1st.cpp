#include "../init.c"
#include "../reconstruction_first.cpp"
#include<stdio.h>

int main(){
  startup();

  if(reconstruction_first_order()) printf("success\n");

  main_grid.print_interfaces();

  return 0;
}
