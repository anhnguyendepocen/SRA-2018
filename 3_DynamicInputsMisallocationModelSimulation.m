%Created by RM on 2018.05.21 to do dynamic programming for Asker et al.
%model

%parameters
epsilon = 4;
beta_discount = 1 / (1+.065); %pg. 1021
lambda = 1; % pg. 1061

%plug in by sector.  for test, try sector = 1
beta_cap_1996 = .5747824;
beta_labor_1996 =  .3394252;
beta_materials_1996 = .4191257;
ar_tfpr_fe_rho_hat = .6236908;
ar_tfpr_fe_mu_hat =   -1.826996;
ar_tfpr_fe_sigma_hat =   1.188126;
 
%Make grid of capital and omega
cap_grid = 0:10:4.470*10^11;



