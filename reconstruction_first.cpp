double tvd_reconstruction(double v, double der_v, double x, double x_i, double delta_x){
  return v + der_v*(x - x_i)/delta_x;
}

template<typename t>
t minmod(t a, t b){
  // sign of number a : 0 when a>=0, 1 when a<0
  bool a_sign = a<0;
  // sign of number b : 0 when b>=0, 1 when b<0
  bool b_sign = b<0;

  // when a and b have same sign and are positive returns their min,
  // when same sign & negative returns their max else 0
  return (1-(a_sign|b_sign))*std::min(a, b) + (a_sign&b_sign)*std::max(a, b);

}

bool reconstruction_first_order(){
  // Iterate over five elements of primitives (Pressure, Density, Velocities)
  // Explain indexes:
  // Iterate over cell centers
  // then we access the interface primitives as 2*(i-dim_b)-1 for the upper limit of the left interface (v^R_{i-1/2})
  // and 2*(i-dim_b)+1 for the lower limit of the right interface (v^L_{i+1/2}) because the array of interfaces doesn't
  // contain ghost cells (-dim_b) and has 2 elements per cell (2*)
  for(int i=0; i<main_grid.i; i++){ 
    // Pressure
    main_grid.interfaces[2*i].p_th = tvd_reconstruction(main_grid.array[i].p_th, minmod(-main_grid.array[i-1].p_th - 1, main_grid.array[i+1].p_th + 1), x[i-1], x[i], d_i);
    main_grid.interfaces[2*i+1].p_th = tvd_reconstruction(main_grid.array[i].p_th, minmod(-main_grid.array[i-1].p_th - 1, main_grid.array[i+1].p_th + 1), x[i+1], x[i], d_i);
    // Density
    main_grid.interfaces[2*i].rho = tvd_reconstruction(main_grid.array[i].rho, minmod(-main_grid.array[i-1].rho - 1, main_grid.array[i+1].rho + 1), x[i-1], x[i], d_i);
    main_grid.interfaces[2*i+1].rho = tvd_reconstruction(main_grid.array[i].rho, minmod(-main_grid.array[i-1].rho - 1, main_grid.array[i+1].rho + 1), x[i+1], x[i], d_i);
    // Velocity (all of its components)
    for (int j = 0; j < n_comp; ++j){
      main_grid.interfaces[2*i].u[j]  = tvd_reconstruction(main_grid.array[i].u[j], minmod(-main_grid.array[i-1].u[j] - 1, main_grid.array[i+1].u[j] + 1), x[i-1], x[i], d_i);
      main_grid.interfaces[2*i+1].u[j]  = tvd_reconstruction(main_grid.array[i].u[j], minmod(-main_grid.array[i-1].u[j] - 1, main_grid.array[i+1].u[j] + 1), x[i+1], x[i], d_i);
    }
  }

  return true;
}
