template <typename T>
T tvd_reconstruction(T v, T der_v, T x, T x_i, T delta_x){
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

// Iterate over five elements of primitives (Pressure, Density, Velocities)
// Explain indexes:
// We stride i by three to iterate over the cell centers
// then we access the interface primitives as i-1 for the upper limit of the left interface (v^R_{i-1/2})
// and i+2 for the lower limit of the right interface (v^L_{i+1/2})
for(int i=2*dim_b+1; i<main_grid.i/3.; i+=3){ 
  // Pressure
  main_grid[i-1].p_th = tvd_reconstruction(v[i].p_th, minmod(-v[i-1].p_th - 1, v[i+1].p_th + 1), x[i-1], x[i], delta_x);
  main_grid[i+1].p_th = tvd_reconstruction(v[i].p_th, minmod(-v[i-1].p_th - 1, v[i+1].p_th + 1), x[i+1], x[i], delta_x);
  // Density
  main_grid[i-1].rho = tvd_reconstruction(v[i].rho, minmod(-v[i-1].rho - 1, v[i+1].rho + 1), x[i-1], x[i], delta_x);
  main_grid[i+1].rho = tvd_reconstruction(v[i].rho, minmod(-v[i-1].rho - 1, v[i+1].rho + 1), x[i+1], x[i], delta_x);
  // Velocity (all of its components)
  for (int j = 0; j < n_comp; ++j){
  main_grid[i-1].u[j]  = tvd_reconstruction(v[i].u[j], minmod(-v[i-1].u[j] - 1, v[i+1].u + 1), x[i-1], x[i], delta_x);
  main_grid[i+1].u[j]  = tvd_reconstruction(v[i].u[j], minmod(-v[i-1].u[j] - 1, v[i+1].u + 1), x[i+1], x[i], delta_x);
  }
}
