#include "char.cpp"
bool llf_solver(){
  // we just sum up the crucial components for the required equation of fllf (21)
{  
  double a = 1;
  for(int i=0; i< idim+2*dim_b; i++){ // kanoume tin allagi opou to main grid exei pleon knai ta boundaries

    for (int j = 0; j < n_comp; j++)
    {
    rieman_speeds(i,j);
    a = fmax(main_grid.eigens[2*i],main_grid.eigens[2*i+1]);
    main_grid.solv_flx[i].rhof_[j] = main_grid.intf_flx[2*i+1].rhof_[j] + main_grid.intf_flx[2*i].rhof_[j]-a*(main_grid.intf_cons[2*i+1].rho_ +main_grid.intf_cons[2*i].rho_);
    main_grid.solv_flx[i].ef_[j] = main_grid.intf_flx[2*i+1].ef_[j] + main_grid.intf_flx[2*i].ef_[j]-a*(main_grid.intf_cons[2*i+1].e_ +main_grid.intf_cons[2*i].e_); 

    for (int k = 0; k < n_comp; k++)
    {
      main_grid.solv_flx[i].pf_[j][k] = main_grid.intf_flx[2*i+1].pf_[j][k] + main_grid.intf_flx[2*i].pf_[j][k]-a*(main_grid.intf_cons[2*i+1].p_[k] +main_grid.intf_cons[2*i].p_[k]);

    }
    }
    }
  }

  return true;
}
