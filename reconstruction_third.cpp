template<typename t> t static
minmod(t a, t b)
{
  // sign of number a : 0 when a>=0, 1 when a<0
  bool a_sign = a<0;
  // sign of number b : 0 when b>=0, 1 when b<0
  bool b_sign = b<0;

  // when a and b have same sign and are positive returns their min,
  // when same sign & negative returns their max else 0
  return (1-(a_sign|b_sign))*std::min(a, b) + (a_sign&b_sign)*std::max(a, b);

}

double static
tvd_limit_rec(double v, double der_v, double x, double x_i, double delta_x)
{
	printf("\nL = %lf\n", v + der_v*(x - x_i)/delta_x);
  return v + der_v*(x - x_i)/delta_x;
}

double static
quad_polyonimals(double v, double v_next, double v_prev, double x, double x_j, double delta_x)
{
	printf("Q = %lf\n", v + (v_next - v_prev)*(x - x_j)/(2*delta_x) + (v_next - 2*v + v_prev)*((x - x_j)/delta_x)*((x - x_j)/delta_x)/2);
	return v + (v_next - v_prev)*(x - x_j)/(2*delta_x) + (v_next - 2*v + v_prev)*((x - x_j)/delta_x)*((x - x_j)/delta_x)/2;
}

double static
weighted_diff(double a, double v, double v_next, double v_prev, double x, double x_j, double delta_x, double limiter)
{
	return a*quad_polyonimals(v, v_next, v_prev, x, x_j, delta_x) - limiter;
}

int static
select_k(double d_minus, double d_zero, double d_plus)
{
	// fprintf(stdout, "delta values: %lf %lf %lf\n", d_minus, d_zero, d_plus);
	bool sign_minus, sign_zero, sign_plus;
	/* Check if all d have same sign */
	sign_minus = d_minus > 0;
	sign_zero = d_zero > 0;
	sign_plus = d_plus > 0;
	if (sign_minus == sign_zero && sign_zero == sign_plus) {
		/* find mininmum delta, return its k */
		if (std::abs(d_minus) < std::min(std::abs(d_zero), std::abs(d_plus))) {
			return -1;
		} else if (std::abs(d_zero) < std::min(std::abs(d_minus), std::abs(d_plus))) {
			return 0;
		} else if (std::abs(d_plus) < std::min(std::abs(d_minus), std::abs(d_zero))) {
			return 1;
		} else {
			fprintf(stderr, "select_k returned bad value\n");
			return 2;
		}
	} else {
		return -2;
	}
}

int
reconstruction_third_order()
{
  // Iterate over five elements of primitives (Pressure, Density, Velocities)
  // Explain indexes:
  // Iterate over cell centers
  // then we access the interface primitives as 2*(i-dim_b)-1 for the upper limit of the left interface (v^R_{i-1/2})
  // and 2*(i-dim_b)+1 for the lower limit of the right interface (v^L_{i+1/2}) because the array of interfaces doesn't
  // contain ghost cells (-dim_b) and has 2 elements per cell (2*)

	int i, j, k;
	int n;
	double d_minus, d_zero, d_plus;
	double l;
  double a;
	double d;

  for (i=dim_b; i< idim+dim_b; i++) { // kanoume tin allagi opou to main grid exei pleon knai ta boundaries

    // Pressure
		l = tvd_limit_rec(main_grid.array[i].p_th, minmod(-main_grid.array[i-1].p_th + 1, main_grid.array[i+1].p_th + 1), x[2*i], x[2*i+1], d_i);

		j = i - 1;
		d_minus = weighted_diff(1, main_grid.array[j].p_th, main_grid.array[j+1].p_th, main_grid.array[j-1].p_th, x[2*j], x[2*j+1], d_i, l);
		j = i;
		d_zero = weighted_diff(0.7, main_grid.array[j].p_th, main_grid.array[j+1].p_th, main_grid.array[j-1].p_th, x[2*j], x[2*j+1], d_i, l);
		j = i + 1;
		d_plus = weighted_diff(1, main_grid.array[j].p_th, main_grid.array[j+1].p_th, main_grid.array[j-1].p_th, x[2*j], x[2*j+1], d_i, l);

		k = select_k(d_minus, d_zero, d_plus);
		printf("k = %d\n", k);

		j = i + k;
		if (k == -2) {
			/* use limiter */
			main_grid.interfaces[2*i].p_th = l;
		} else if (k == 2) {
			return -1;
		} else {
			/* Use quadratic */
			main_grid.interfaces[2*i].p_th = quad_polyonimals(main_grid.array[j].p_th, main_grid.array[j+1].p_th, main_grid.array[j-1].p_th, x[2*j], x[2*j+1], d_i);
		}

		l = tvd_limit_rec(main_grid.array[i].p_th, minmod(-main_grid.array[i-1].p_th + 1, main_grid.array[i+1].p_th + 1), x[2*(i+1)], x[i], d_i);

		j = i - 1;
		d_minus = weighted_diff(1, main_grid.array[j].p_th, main_grid.array[j+1].p_th, main_grid.array[j-1].p_th, x[2*(j+1)], x[j], d_i, l);
		j = i;
		d_zero = weighted_diff(0.7, main_grid.array[j].p_th, main_grid.array[j+1].p_th, main_grid.array[j-1].p_th, x[2*(j+1)], x[j], d_i, l);
		j = i + 1;
		d_plus = weighted_diff(1, main_grid.array[j].p_th, main_grid.array[j+1].p_th, main_grid.array[j-1].p_th, x[2*(j+1)], x[j], d_i, l);

		k = select_k(d_minus, d_zero, d_plus);
		printf("k = %d\n", k);
		j = i + k;
		if (k == -2) {
			/* use limiter */
			main_grid.interfaces[2*i+1].p_th = l;
		} else if (k == 2) {
			return -1;
		} else {
			/* Use quadratic */
			main_grid.interfaces[2*i+1].p_th = quad_polyonimals(main_grid.array[j].p_th, main_grid.array[j+1].p_th, main_grid.array[j-1].p_th, x[2*(j+1)], x[j], d_i);
		}
  
//     // Density
// 		j = i - 1;
// 		d_minus = weighted_diff(1, main_grid.array[i].rho, main_grid.array[i+1].rho, main_grid.array[i-1].rho, minmod(-main_grid.array[i-1].rho + 1, main_grid.array[i+1].rho + 1), x[2*i], x[2*i+1], x[2*j+1], d_i);
// 		j = i;
// 		d_zero = weighted_diff(0.7, main_grid.array[i].rho, main_grid.array[i+1].rho, main_grid.array[i-1].rho, minmod(-main_grid.array[i-1].rho + 1, main_grid.array[i+1].rho + 1), x[2*i], x[2*i+1], x[2*j+1], d_i);
// 		j = i + 1;
// 		d_plus = weighted_diff(1, main_grid.array[i].rho, main_grid.array[i+1].rho, main_grid.array[i-1].rho, minmod(-main_grid.array[i-1].rho + 1, main_grid.array[i+1].rho + 1), x[2*i], x[2*i+1], x[2*j+1], d_i);
// 
// 		k = select_k(d_minus, d_zero, d_plus);
// 		j = i + k;
// 		if (k == -2) {
// 			/* use limiter */
// 			main_grid.interfaces[2*i].rho = tvd_limit_rec(main_grid.array[i].rho, minmod(-main_grid.array[i-1].rho + 1, main_grid.array[i+1].rho + 1), x[2*i], x[2*i+1], d_i);
// 			main_grid.interfaces[2*i+1].rho = tvd_limit_rec(main_grid.array[i].rho, minmod(-main_grid.array[i-1].rho + 1, main_grid.array[i+1].rho + 1), x[2*(i+1)], x[i], d_i);
// 		} else if (k == 2) {
// 			return -1;
// 		} else {
// 			/* Use quadratic */
// 			main_grid.interfaces[2*i].rho = quad_polyonimals(main_grid.array[i].rho, main_grid.array[i+1].rho, main_grid.array[i-1].rho, x[2*i], x[2*j+1], d_i);
// 			main_grid.interfaces[2*i+1].rho = quad_polyonimals(main_grid.array[i].rho, main_grid.array[i+1].rho, main_grid.array[i-1].rho, x[2*(i+1)], x[i], d_i);
// 		}
// 
//     // Velocity (all of its components)
//     for (int n = 0; n < n_comp; ++n){
// 			j = i - 1;
// 			d_minus = weighted_diff(1, main_grid.array[i].u[n], main_grid.array[i+1].u[n], main_grid.array[i-1].u[n], minmod(-main_grid.array[i-1].u[n] + 1, main_grid.array[i+1].u[n] + 1), x[2*i], x[2*i+1], x[2*j+1], d_i);
// 			j = i;
// 			d_zero = weighted_diff(0.7, main_grid.array[i].u[n], main_grid.array[i+1].u[n], main_grid.array[i-1].u[n], minmod(-main_grid.array[i-1].u[n] + 1, main_grid.array[i+1].u[n] + 1), x[2*i], x[2*i+1], x[2*j+1], d_i);
// 			j = i + 1;
// 			d_plus = weighted_diff(1, main_grid.array[i].u[n], main_grid.array[i+1].u[n], main_grid.array[i-1].u[n], minmod(-main_grid.array[i-1].u[n] + 1, main_grid.array[i+1].u[n] + 1), x[2*i], x[2*i+1], x[2*j+1], d_i);
// 
// 			k = select_k(d_minus, d_zero, d_plus);
// 			j = i + k;
// 			if (k == -2) {
// 				/* use limiter */
// 				main_grid.interfaces[2*i].u[n] = tvd_limit_rec(main_grid.array[i].u[n], minmod(-main_grid.array[i-1].u[n] + 1, main_grid.array[i+1].u[n] + 1), x[2*i], x[2*i+1], d_i);
// 				main_grid.interfaces[2*i+1].u[n] = tvd_limit_rec(main_grid.array[i].u[n], minmod(-main_grid.array[i-1].u[n] + 1, main_grid.array[i+1].u[n] + 1), x[2*(i+1)], x[i], d_i);
// 			} else if (k == 2) {
// 				return -1;
// 			} else {
// 				/* Use quadratic */
// 				main_grid.interfaces[2*i].u[n] = quad_polyonimals(main_grid.array[i].u[n], main_grid.array[i+1].u[n], main_grid.array[i-1].u[n], x[2*i], x[2*j+1], d_i);
// 				main_grid.interfaces[2*i+1].u[n] = quad_polyonimals(main_grid.array[i].u[n], main_grid.array[i+1].u[n], main_grid.array[i-1].u[n], x[2*(i+1)], x[i], d_i);
// 			}
// 
//     }

			/* TODO name this part of the code */
			main_grid.intf_cons[2*i].build(main_grid.interfaces[2*i]);
			main_grid.intf_cons[2*i+1].build(main_grid.interfaces[2*i+1]);
  }

  return 0;
}
