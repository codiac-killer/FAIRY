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

// Iterate over three elements of primatives
// v_left[0] = tvd_reconstruction(v[i].p_th, minmod(-v[i-1].p_th - 1, v[i+1].p_th + 1), x, x[i], delta_x);
// v_right[0] = tvd_reconstruction(v[i].p_th, minmod(-v[i-1].p_th - 1, v[i+1].p_th + 1), x, x[i], delta_x);
// v_left[1] = tvd_reconstruction(v[i].rho, minmod(-v[i-1].rho - 1, v[i+1].rho + 1), x, x[i], delta_x);
// v_right[1] = tvd_reconstruction(v[i].rho, minmod(-v[i-1].rho - 1, v[i+1].rho + 1), x, x[i], delta_x);
// v_left[2] = tvd_reconstruction(v[i].u, minmod(-v[i-1].u - 1, v[i+1].u + 1), x, x[i], delta_x);
// v_right[2] = tvd_reconstruction(v[i].u, minmod(-v[i-1].u - 1, v[i+1].u + 1), x, x[i], delta_x);

