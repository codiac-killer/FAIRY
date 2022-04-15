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
tvd_limit_rec(double v, double der_v, double x_if, double x_center, double delta_x)
{
	printf("\nL = %lf\n", v + der_v*(x_if - x_center)/delta_x);
  return v + der_v*(x_if - x_center)/delta_x;
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

static int
reconstruct_pressure(const int i)
{
	int j,k;
	double l, d_minus, d_zero, d_plus;

	l = tvd_limit_rec(main_grid.array[2*(i)+1].p_th, minmod(-main_grid.array[2*(i-1)+1].p_th + 1, main_grid.array[2*(i+1)+1].p_th + 1), x[2*i], x[2*i+1], d_i);

	j = i - 1;
	d_minus = weighted_diff(1, main_grid.array[2*(j)+1].p_th, main_grid.array[2*(j+1)+1].p_th, main_grid.array[2*(j-1)+1].p_th, x[2*j], x[2*j+1], d_i, l);
	j = i;
	d_zero = weighted_diff(0.7, main_grid.array[2*(j)+1].p_th, main_grid.array[2*(j+1)+1].p_th, main_grid.array[2*(j-1)+1].p_th, x[2*j], x[2*j+1], d_i, l);
	j = i + 1;
	d_plus = weighted_diff(1, main_grid.array[2*(j)+1].p_th, main_grid.array[2*(j+1)+1].p_th, main_grid.array[2*(j-1)+1].p_th, x[2*j], x[2*j+1], d_i, l);

	k = select_k(d_minus, d_zero, d_plus);
	printf("k = %d, i = %d\n", k, i);

	j = i + k;
	if (k == -2) {
		/* use limiter */
		main_grid.interfaces[2*i].p_th = l;
	} else if (k == 2) {
		return -1;
	} else {
		/* Use quadratic */
		main_grid.interfaces[2*i].p_th = quad_polyonimals(main_grid.array[2*(j)+1].p_th, main_grid.array[2*(j+1)+1].p_th, main_grid.array[2*(j-1)+1].p_th, x[2*j], x[2*j+1], d_i);
	}

	l = tvd_limit_rec(main_grid.array[2*(i)+1].p_th, minmod(-main_grid.array[2*(i-1)+1].p_th + 1, main_grid.array[2*(i+1)+1].p_th + 1), x[2*i+2], x[2*i+1], d_i);

	j = i - 1;
	d_minus = weighted_diff(1, main_grid.array[2*(j)+1].p_th, main_grid.array[2*(j+1)+1].p_th, main_grid.array[2*(j-1)+1].p_th, x[2*j+2], x[2*j+1], d_i, l);
	j = i;
	d_zero = weighted_diff(0.7, main_grid.array[2*(j)+1].p_th, main_grid.array[2*(j+1)+1].p_th, main_grid.array[2*(j-1)+1].p_th, x[2*j+2], x[2*j+1], d_i, l);
	j = i + 1;
	d_plus = weighted_diff(1, main_grid.array[2*(j)+1].p_th, main_grid.array[2*(j+1)+1].p_th, main_grid.array[2*(j-1)+1].p_th, x[2*j+2], x[2*j+1], d_i, l);

	k = select_k(d_minus, d_zero, d_plus);
	printf("k = %d, i = %d\n", k, i);
	j = i + k;
	if (k == -2) {
		/* use limiter */
		main_grid.interfaces[2*i+1].p_th = l;
	} else if (k == 2) {
		return -1;
	} else {
		/* Use quadratic */
		main_grid.interfaces[2*i+1].p_th = quad_polyonimals(main_grid.array[2*(j)+1].p_th, main_grid.array[2*(j+1)+1].p_th, main_grid.array[2*(j-1)+1].p_th, x[2*j+2], x[2*j+1], d_i);
	}

	return 0;

}

static int
reconstruct_density(int i)
{
	int j,k;
	double l, d_minus, d_zero, d_plus;

	l = tvd_limit_rec(main_grid.array[i].rho, minmod(-main_grid.array[i-1].rho + 1, main_grid.array[i+1].rho + 1), x[2*i], x[2*i+1], d_i);

	j = i - 1;
	d_minus = weighted_diff(1, main_grid.array[j].rho, main_grid.array[j+1].rho, main_grid.array[j-1].rho, x[2*j], x[2*j+1], d_i, l);
	j = i;
	d_zero = weighted_diff(0.7, main_grid.array[j].rho, main_grid.array[j+1].rho, main_grid.array[j-1].rho, x[2*j], x[2*j+1], d_i, l);
	j = i + 1;
	d_plus = weighted_diff(1, main_grid.array[j].rho, main_grid.array[j+1].rho, main_grid.array[j-1].rho, x[2*j], x[2*j+1], d_i, l);

	k = select_k(d_minus, d_zero, d_plus);
	printf("k = %d\n", k);

	j = i + k;
	if (k == -2) {
		/* use limiter */
		main_grid.interfaces[2*i].rho = l;
	} else if (k == 2) {
		return -1;
	} else {
		/* Use quadratic */
		main_grid.interfaces[2*i].rho = quad_polyonimals(main_grid.array[j].rho, main_grid.array[j+1].rho, main_grid.array[j-1].rho, x[2*j], x[2*j+1], d_i);
	}

	l = tvd_limit_rec(main_grid.array[i].rho, minmod(-main_grid.array[i-1].rho + 1, main_grid.array[i+1].rho + 1), x[2*i+2], x[2*i+1], d_i);

	j = i - 1;
	d_minus = weighted_diff(1, main_grid.array[j].rho, main_grid.array[j+1].rho, main_grid.array[j-1].rho, x[2*j+2], x[2*j+1], d_i, l);
	j = i;
	d_zero = weighted_diff(0.7, main_grid.array[j].rho, main_grid.array[j+1].rho, main_grid.array[j-1].rho, x[2*j+2], x[2*j+1], d_i, l);
	j = i + 1;
	d_plus = weighted_diff(1, main_grid.array[j].rho, main_grid.array[j+1].rho, main_grid.array[j-1].rho, x[2*j+2], x[2*j+1], d_i, l);

	k = select_k(d_minus, d_zero, d_plus);
	printf("k = %d\n", k);
	j = i + k;
	if (k == -2) {
		/* use limiter */
		main_grid.interfaces[2*i+1].rho = l;
	} else if (k == 2) {
		return -1;
	} else {
		/* Use quadratic */
		main_grid.interfaces[2*i+1].rho = quad_polyonimals(main_grid.array[j].rho, main_grid.array[j+1].rho, main_grid.array[j-1].rho, x[2*j+2], x[2*j+1], d_i);
	}

	return 0;

}

static int
reconstruct_velocity(int i, int n)
{
	int j,k;
	double l, d_minus, d_zero, d_plus;

	l = tvd_limit_rec(main_grid.array[i].u[n], minmod(-main_grid.array[i-1].u[n] + 1, main_grid.array[i+1].u[n] + 1), x[2*i], x[2*i+1], d_i);

	j = i - 1;
	d_minus = weighted_diff(1, main_grid.array[j].u[n], main_grid.array[j+1].u[n], main_grid.array[j-1].u[n], x[2*j], x[2*j+1], d_i, l);
	j = i;
	d_zero = weighted_diff(0.7, main_grid.array[j].u[n], main_grid.array[j+1].u[n], main_grid.array[j-1].u[n], x[2*j], x[2*j+1], d_i, l);
	j = i + 1;
	d_plus = weighted_diff(1, main_grid.array[j].u[n], main_grid.array[j+1].u[n], main_grid.array[j-1].u[n], x[2*j], x[2*j+1], d_i, l);

	k = select_k(d_minus, d_zero, d_plus);
	printf("k = %d\n", k);

	j = i + k;
	if (k == -2) {
		/* use limiter */
		main_grid.interfaces[2*i].u[n] = l;
	} else if (k == 2) {
		return -1;
	} else {
		/* Use quadratic */
		main_grid.interfaces[2*i].u[n] = quad_polyonimals(main_grid.array[j].u[n], main_grid.array[j+1].u[n], main_grid.array[j-1].u[n], x[2*j], x[2*j+1], d_i);
	}

	l = tvd_limit_rec(main_grid.array[i].u[n], minmod(-main_grid.array[i-1].u[n] + 1, main_grid.array[i+1].u[n] + 1), x[2*i+2], x[2*i+1], d_i);

	j = i - 1;
	d_minus = weighted_diff(1, main_grid.array[j].u[n], main_grid.array[j+1].u[n], main_grid.array[j-1].u[n], x[2*j+2], x[2*j+1], d_i, l);
	j = i;
	d_zero = weighted_diff(0.7, main_grid.array[j].u[n], main_grid.array[j+1].u[n], main_grid.array[j-1].u[n], x[2*j+2], x[2*j+1], d_i, l);
	j = i + 1;
	d_plus = weighted_diff(1, main_grid.array[j].u[n], main_grid.array[j+1].u[n], main_grid.array[j-1].u[n], x[2*j+2], x[2*j+1], d_i, l);

	k = select_k(d_minus, d_zero, d_plus);
	printf("k = %d\n", k);
	j = i + k;
	if (k == -2) {
		/* use limiter */
		main_grid.interfaces[2*i+1].u[n] = l;
	} else if (k == 2) {
		return -1;
	} else {
		/* Use quadratic */
		main_grid.interfaces[2*i+1].u[n] = quad_polyonimals(main_grid.array[j].u[n], main_grid.array[j+1].u[n], main_grid.array[j-1].u[n], x[2*j+2], x[2*j+1], d_i);
	}

	return 0;

}

int
reconstruction_third_order(void)
{
  // Iterate over five elements of primitives (Pressure, Density, Velocities)
  // Explain indexes:
  // Iterate over cell centers
  // then we access the interface primitives as 2*(i-dim_b)-1 for the upper limit of the left interface (v^R_{i-1/2})
  // and 2*(i-dim_b)+1 for the lower limit of the right interface (v^L_{i+1/2}) because the array of interfaces doesn't
  // contain ghost cells (-dim_b) and has 2 elements per cell (2*)

	int i;
	int n;

  for (i=dim_b; i< idim+dim_b; i++) { // kanoume tin allagi opou to main grid exei pleon knai ta boundaries

    // Pressure
		if (reconstruct_pressure(i) < 0) 
			return -1;
  
    // Density
// 		if (reconstruct_density(i) < 0)
// 			return -1;

    // Velocity (all of its components)
//     for (int n = 0; n < n_comp; ++n){
// 			if (reconstruct_velocity(i, n) < 0)
// 				return -1;
//     }

			/* Build conservatives ? */
			main_grid.intf_cons[2*i].build(main_grid.interfaces[2*i]);
			main_grid.intf_cons[2*i+1].build(main_grid.interfaces[2*i+1]);
  }

  return 0;
}
