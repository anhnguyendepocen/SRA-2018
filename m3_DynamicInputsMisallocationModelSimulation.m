%Created by RM on 2018.05.21 to do dynamic programming for Asker et al.
%model

%parameters
epsilon = -4;
beta_discount = 1 / (1+.065); %pg. 1021
lambda = 1; % pg. 1061
delta = .1; %pg. 1061

%plug in by sector.  for test, try sector = 1
beta_cap_1996 = .5747824;
beta_labor_1996 =  .3394252;
beta_materials_1996 = .4191257;
ar_tfpr_fe_rho_hat = .6236908;
ar_tfpr_fe_mu_hat =   -1.826996;
ar_tfpr_fe_sigma_hat =   1.188126;
 
%%Make grid of capital and omega
lower_cap = 5000;
upper_cap =  4*10^11;
sd_cap =  2.35e+10;
cap_int = 1*10^9;
cap_grid = 2500:cap_int:upper_cap + 1.5*sd_cap; %for testing using 


lower_omega = -10;
upper_omega = .8;
sd_omega = 1.5;
omega_grid = lower_omega -  sd_omega:1:upper_omega +1.5*sd_omega; %for testing

pre_size_cap = size(cap_grid);
size_cap = pre_size_cap(1,2);
pre_size_omega = size(omega_grid);
size_omega = pre_size_omega(1,2);

%find starting values for capital and omega
cap_start = 0;
cap_end = 0;
for i = 1:size_cap;
     if cap_grid(1,i) < lower_cap 
          cap_start = i;
     end
     if cap_grid(1,i) < upper_cap 
         cap_end = i + 1;
     end
end;

omega_start = 0;
omega_end = 0;
for j = 1:size_omega;
     if omega_grid(1,j) < lower_omega
          omega_start = j;
     end
     if omega_grid(1,j) < upper_omega
         omega_end = j + 1;
     end
end;     


%%Inititalize present period profit

profitspace = zeros(size_cap,size_omega);

Omega_exponent = 1 / (beta_cap_1996 + epsilon^(-1));
cap_exponent = beta_cap_1996 / (beta_cap_1996 + epsilon^(-1));

for i = 1:size_cap;
    for j =1:size_omega;
        profit = lambda * exp(omega_grid(1,j))^Omega_exponent*cap_grid(1,i)^cap_exponent;
        profitspace(i,j) = profit;
    end;
end;

%Initialize Cost Function Constant

for costf = 0:.01:.4 
   for costv = 0:.01:.4

cost_fixed = costf;
cost_variable = costv;

print_cost_fixed = cost_fixed * 1000;
print_cost_variable = cost_variable * 1000;

%%Estimation Matri

%%%%%Start Value Function iteration

%%Initalize Value Function Matrix with Guess
v = zeros(size_cap,size_omega) + 10000;

for i = 1:size_cap ;
    for j=1:size_omega ;
        vexp = exp(omega_grid(1,j)+15)*2;
        vadd = sqrt(cap_grid(1,i)^(1.75)*sqrt(vexp));
        v(i,j) = v(i,j) + (vadd - v(i,j))/10000000;
    end;
end;

vinit = v;
v_update_matrix = v;

optimal_invest = zeros(size_cap,size_omega);

for iter = 1:500;
    iter

 v_update_matrix = vinit;
 v = v_update_matrix;

for i = cap_start:cap_end;
%itest= 200;
%jtest = 5;
%for i = itest:itest;
    for j=omega_start:omega_end;
    %for j = jtest:jtest;
       %Find Transition Probabilities
       transition_total = 0;
       omega_current = omega_grid(1,j);
       transition_unscaled = omega_grid * 0;       
       %Now Iterate Over  Omegas
        for k = 1:size_omega;
               transition_unscaled_target = omega_grid(1,k) - ar_tfpr_fe_rho_hat * omega_current - ar_tfpr_fe_mu_hat;
               transition_unscaled(1,k) = (1/sqrt(2*pi*ar_tfpr_fe_sigma_hat*ar_tfpr_fe_sigma_hat))*exp(-(transition_unscaled_target*transition_unscaled_target)/(2*ar_tfpr_fe_sigma_hat*ar_tfpr_fe_sigma_hat));
               transition_total = transition_total + transition_unscaled(1,k);
        end;
        
        transition_scaled = transition_unscaled / transition_total;
       
       %%Now Find Optimal Investment Given v
       cap_current = cap_grid(1,i);     
       best_next_cap = cap_grid * 0;
       best_next_cap_init = best_next_cap; %for testing
       
       best_next_cap_pre_adj_cost = best_next_cap_init; %for testing

       for l = 1:size_cap;
           %start with continuation value at capital l
           %best_next_cap = cap_grid * 0;
           
           for m = 1:size_omega;
              %continuation value for capital choice at index l, over
              %distribution of productivity
             best_next_cap(1,l) =  best_next_cap(1,l) + beta_discount * transition_scaled(1,m) * v(l,m);
              if l < 2
                  if m < 2
                      test =  best_next_cap(1,l);
                  end;
              end;
           end;
           best_next_cont_value = best_next_cap;
               
          %now add in value of current profit less adjustment cost for
          %begin at capital value l
           
           pre_investment = cap_grid(1,l) - (1-delta)*cap_current;
           investment = max(0,pre_investment);
           invest_dummy = 0;
           
           if pre_investment > 0 
               invest_dummy = 1;
           end;
           
           adj_inv_cost = investment + invest_dummy * profitspace(i,j) * cost_fixed + (investment * investment * cost_variable) / cap_current;
           best_next_cap(1,l) = best_next_cap(1,l) + profitspace(i,j);
           best_next_cap_pre_adj_cost(1,l) =  best_next_cap(1,l); %for testing
           best_next_cap(1,l) =  best_next_cap(1,l)  -  adj_inv_cost;     
      
           end;
           
           %%Now find best investment choice
           best_cap_choice = 0;
           best_next_cap_tracker = 0;
          for n= 1:size_cap
              if best_next_cap(1,n) > best_next_cap_tracker
                  best_cap_choice = n;
                  best_next_cap_tracker = best_next_cap(1,n);
              end;
          end;
           
          optimal_invest(i,j) = best_next_cap_tracker;
          
           %%Now Replace the Value Function
           v_update_value = 0;
           %start with continuation value
           for n = 1:size_omega;
             v_update_value =  v_update_value + beta_discount * transition_scaled(1,n) * v(best_cap_choice,n);
           end;
           v_update_cont_value = v_update_value;
           pre_investment_chosen = best_next_cap_tracker - (1-delta)*cap_current;
           investment_chosen = max(0,pre_investment_chosen);
           investment_chosen_dummy = 0;
            if investment_chosen > 0 
               investment_chosen_dummy = 1;
           end;
           v_update_value = v_update_value + profitspace(i,j) - investment_chosen;
           v_update_value = v_update_value - investment_chosen_dummy * profitspace(i,j) * cost_fixed;
           v_update_value = v_update_value - (investment_chosen * investment_chosen * cost_variable) / cap_current;
           
           v_update_matrix(i,j) = max(v_update_value,0);
  
    end;
end;

end;


cd '/Users/russellmorton/Desktop/SRA 2018/STATA/Output/Dynamic Inputs Misallocation';
OutFileName = ['OptimalInvst_Fixed' num2str(print_cost_fixed) '_Variable' num2str(print_cost_variable) '.xlsx'];

xlswrite(OutFileName,optimal_invest);

   end;
end;



%val = 104.0102
%val = v(itest,jtest);
%val_update = v_update_matrix(itest,jtest);



%Remember that picking investment is just like picking next period capital







