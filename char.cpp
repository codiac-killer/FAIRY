bool rieman_speeds(int i, int j ){
double u_n,h_enth,cs;
double u_p,u_t = 0;
double eig_l[4];
for (int m = 0; m < 2; ++m)
{
	u_p = main_grid.interfaces[2*i+m].u[j];
	for (int l = 0; l < n_comp; ++l)
	{
		u_n = u_p+(l!=j)*pow(main_grid.interfaces[2*i+m].u[l],2);
		u_t = u_t+pow(main_grid.interfaces[2*i+m].u[l],2);
	}
	u_n = pow(u_p,0.5);
	u_t = pow(u_t,0.5);
	h_enth = gamma_par/(gamma_par-1)*main_grid.interfaces[2*i+m].p_th/main_grid.interfaces[2*i+m].rho;
	cs     = gamma_par/h_enth*main_grid.interfaces[2*i+m].p_th/main_grid.interfaces[2*i+m].rho;
	eig_l[2*m] = (u_p*(1-pow(cs,2))+cs*pow((1-u_t*u_t)*(1-u_p*u_p-pow(u_n*cs,2)),0.5))/(1-pow(cs*u_t,2));
	eig_l[2*m+1] = (u_p*(1-pow(cs,2))-cs*pow((1-u_t*u_t)*(1-u_p*u_p-pow(u_n*cs,2)),0.5))/(1-pow(cs*u_t,2));
	

}

//main_grid.eigens[2*i]   = std::max({0.0,eig_l[0],eig_l[2]});
//main_grid.eigens[2*i+1] = std::max({0.0,eig_l[1],eig_l[3]});

main_grid.eigens[2*i]   = std::max({1.,2.,3.});



return true;
}