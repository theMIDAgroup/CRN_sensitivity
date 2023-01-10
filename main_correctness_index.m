clc
clear
close all

% This code creates the figures representing the correctness of the SSI

set(0, 'defaultAxesTickLabelInterpreter','latex'); 
set(0, 'defaultLegendInterpreter','latex'); 
warning('off', 'all')

addpath('./funcs')

%% Step 1. Load structure containing information on species and reactions
file_crn = fullfile('data', 'CRC_CRN_nodrug.mat');
load(file_crn)

%% Step 2. Import the initial data
vm = new_CMIM.matrix.v;
n_species = numel(new_CMIM.species.names);
n_reactions = size(new_CMIM.matrix.S, 2);
n_cons_laws = size(new_CMIM.matrix.Nl, 1);
ind_one = new_CMIM.matrix.ind_one;

k_values = new_CMIM.rates.std_values;
Sm = new_CMIM.matrix.S;
x_0 = new_CMIM.species.std_initial_values;
idx_basic_species = find(x_0>0);
cons_laws = new_CMIM.matrix.Nl;
rho=cons_laws*x_0;
%% Step 3. Compute the values of the SSI

%% 3.1. Compute equilibrium
max_counter = 300;
ris_phys = f_NLPC_restart(x_0, k_values, Sm, cons_laws, rho, idx_basic_species, ...
           vm, ind_one, max_counter, 0);
x_values = ris_phys.x;

%% 3.2. Analytically compute the sensitivity matrix
% Compute the analytic jacobian of v 
jacobian_v = f_compute_analytic_jacobian_v(vm, n_species, ind_one);

% Calculate the Jacobian of f with respect to x
eval_jac = f_evaluate_jacobian_neworder(k_values, x_values, ...
    Sm, idx_basic_species, jacobian_v, cons_laws);

% Calculate the Jacobian of f with respect to k
eval_jac_k=f_evaluate_jacobian_k(vm, x_values, Sm, idx_basic_species);

% Calculate Jacobian of f with respect to c
eval_jac_c=f_evaluate_jacobian_c(n_species, n_cons_laws,0);

% Calculate sensitivity matrix w.r.t. rate constants k
Jac_k_xeq= -inv(eval_jac)*eval_jac_k;

% Calcolate sensitivity matrix w.r.t. conservation laws' constants c
Jac_c_xeq= -inv(eval_jac)*eval_jac_c;

%% 3.3. Compute SSI w.r.t. rate constants k 
L_k=(1./x_values).*(Jac_k_xeq).*(k_values');
[U_k,D_k,V_k]=svd(L_k);
e_PCA_k =  ((V_k(:,1:n_species).^2)*diag(D_k.^2))/sum(diag(D_k.^2)); %se presi con il quadrato: e_j=(L^T*L)_{jj}/sum_k(lambda_k)

%% 3.4. Compute SSI w.r.t. conservation lawss' constant c
L_c=(1./x_values).*(Jac_c_xeq).*(rho');
[U_c,D_c,V_c]=svd(L_c);
e_PCA_c =  ((V_c(:,1:n_cons_laws).^2)*(diag(D_c).^2))/sum(diag(D_c).^2); %se presi con il quadrato: e_j=(L^T*L)_{jj}/sum_k(lambda_k)

%% Graphs: 
delta_j = [0 1 2 5 10 20];  
n_j = numel(delta_j);
max_counter = 300;

%% Fig 1a:
% Impact of increasing kinetic parameters referred to
% maximum
% minimum
% median
% SSI on global CRC-CRN
[max_e_k, idx_k_high] = max(e_PCA_k); 
[min_e_k, idx_k_low] = min(e_PCA_k);
aux=median(e_PCA_k); idx_k_mean=find(e_PCA_k==aux);

x_eq_high_k = zeros(n_species, n_j);
x_eq_low_k = zeros(n_species, n_j);
x_eq_mean_k = zeros(n_species, n_j);

ind_one=n_species+1;
for ik = 1:n_j % Compute the equilibrium by varying the selected rate constants
   
   fprintf('Value of k num = %d \n', ik) % First rate constant
   rates = new_CMIM.rates.std_values;
   rates(idx_k_high) = rates(idx_k_high)*(delta_j(ik)+1);  
   ris = f_NLPC_restart(x_0, rates, Sm, cons_laws, rho, idx_basic_species, ...
           vm, ind_one, max_counter, 0);
   x_eq_high_k(:, ik) = ris.x;

   clear rates ris

   rates = new_CMIM.rates.std_values;    % Second rate constant
   rates(idx_k_low) = rates(idx_k_low)*(delta_j(ik)+1);  
   ris = f_NLPC_restart(x_0, rates, Sm, cons_laws, rho, idx_basic_species, ...
           vm, ind_one, max_counter, 0);
   x_eq_low_k(:, ik) = ris.x;


   clear rates ris

   rates = new_CMIM.rates.std_values;    %Third rate constant
   rates(idx_k_mean) = rates(idx_k_mean)*(delta_j(ik)+1);   
   ris = f_NLPC_restart(x_0, rates, Sm, cons_laws, rho, idx_basic_species, ...
           vm, ind_one, max_counter, 0);
   x_eq_mean_k(:, ik) = ris.x;

   clear rates ris
   

end

delta_x_eq_high_k = (x_eq_high_k - x_eq_high_k(:, 1)) ./ x_eq_high_k(:, 1);
norm_delta_x_high_k = vecnorm(delta_x_eq_high_k, 2, 1);

delta_x_eq_low_k = (x_eq_low_k - x_eq_low_k(:, 1)) ./ x_eq_low_k(:, 1);
norm_delta_x_low_k = vecnorm(delta_x_eq_low_k, 2, 1);

delta_x_eq_mean_k = (x_eq_mean_k - x_eq_mean_k(:, 1)) ./ x_eq_mean_k(:, 1);
norm_delta_x_mean_k = vecnorm(delta_x_eq_mean_k, 2, 1);

%% Fig 1b:
% Impact of increasing conservation laws' constants referred to
% maximum
% minimum
% median
% SSI on global CRC-CRN
[max_e_c, idx_c_high] = max(e_PCA_c);
[min_e_c, idx_c_low] = min(e_PCA_c);
aux=median(e_PCA_c); idx_c_mean=find(e_PCA_c==aux);

x_eq_high_c = zeros(n_species, n_j);
x_eq_low_c = zeros(n_species, n_j);
x_eq_mean_c = zeros(n_species, n_j);

for ic = 1:n_j % Compute the equilibrium by varying the selected rate constants
   
   fprintf('Value of c num = %d \n', ic) % First rate constant
   rho = cons_laws*x_0;
   rho(idx_c_high) = rho(idx_c_high)*(delta_j(ic)+1);
   ris = f_NLPC_restart(x_0, k_values, Sm, cons_laws, rho, idx_basic_species, ...
           vm, ind_one, max_counter, 0);
   x_eq_high_c(:, ic) = ris.x;

   clear rates ris

   rho = cons_laws*x_0  ;  % Second rate constant
   rho(idx_c_low) = rho(idx_c_low)*(delta_j(ic)+1);
   ris = f_NLPC_restart(x_0, k_values, Sm, cons_laws, rho, idx_basic_species, ...
           vm, ind_one, max_counter, 0);
   x_eq_low_c(:, ic) = ris.x;


clear rates ris

   rho = cons_laws*x_0;    %Third rate constant
   rho(idx_c_mean) = rho(idx_c_mean)*(delta_j(ic)+1);
   ris = f_NLPC_restart(x_0, k_values, Sm, cons_laws, rho, idx_basic_species, ...
           vm, ind_one, max_counter, 0);
   x_eq_mean_c(:, ic) = ris.x;

   clear rates ris
   

end

delta_x_eq_high = (x_eq_high_c - x_eq_high_c(:, 1)) ./ x_eq_high_c(:, 1);
norm_delta_x_high = vecnorm(delta_x_eq_high, 2, 1);

delta_x_eq_low = (x_eq_low_c - x_eq_low_c(:, 1)) ./ x_eq_low_c(:, 1);
norm_delta_x_low = vecnorm(delta_x_eq_low, 2, 1);

delta_x_eq_mean = (x_eq_mean_c - x_eq_mean_c(:, 1)) ./ x_eq_mean_c(:, 1);
norm_delta_x_mean = vecnorm(delta_x_eq_mean, 2, 1);

% Figures
figure('units','normalized','outerposition',[0 0 1 0.7]);
subplot(1, 2 ,1)
plot(delta_j, norm_delta_x_high_k, 'k-*', 'Linewidth', 3)
hold on
plot(delta_j, norm_delta_x_mean_k, 'b-*', 'Linewidth', 3)
hold on
plot(delta_j, norm_delta_x_low_k, 'r-*', 'Linewidth', 3)

xlabel('$\hat{\Delta} k$', 'Interpreter', 'Latex')
ylabel('$\sqrt{Q(\hat{\Delta} k)}$', 'Interpreter', 'Latex')
xticks(delta_j);
my_symlog('y')
set(gca, 'Fontsize', 15)
grid on

subplot(1,2, 2)
plot(delta_j, norm_delta_x_high, 'k-*', 'Linewidth', 3)
hold on
plot(delta_j, norm_delta_x_mean, 'b-*', 'Linewidth', 3)
hold on
plot(delta_j, norm_delta_x_low, 'r-*', 'Linewidth', 3)
my_symlog('y')
% set(gca, 'Yscale', 'log')
xlabel('$\hat{\Delta} c$', 'Interpreter', 'Latex')
ylabel('$\sqrt{Q(\hat{\Delta} c)}$', 'Interpreter', 'Latex')
set(gca, 'Fontsize', 15)
xticks(delta_j)
grid on
legend({'Maximum SSI';'Median SSI';'Minimum SSI'}, 'Location', 'EastOutside')

%% Fig 2:
% Impact of increasing kinetic parameters referred to
% maximum
% minimum
% median
% SSI on TP53 concentration

selected_proteins = {'TP53'};
[aux_, idx_proteins] = ismember(selected_proteins, new_CMIM.species.names);

selpart_L = L_c(idx_proteins, :);
[U, Sigma, V]=svd(selpart_L);
selection=size(selpart_L,1);
eing=diag(Sigma(1:selection, 1:selection));
sum_eig=sum(eing.^2); 
selpart_e_PCA=((V(:,1:selection).^2)*(eing.^2))/sum_eig;

[max_e, idx_c_high] = max(selpart_e_PCA);
[min_e, idx_c_low] = min(selpart_e_PCA);
aux=median(selpart_e_PCA); idx_c_mean=find(selpart_e_PCA==aux);

selpart_x_eq_high = zeros(n_species, n_j);
selpart_x_eq_low = zeros(n_species, n_j);
selpart_x_eq_mean=zeros(n_species, n_j);

for ic = 1:n_j % Compute the equilibrium by varying the selected rate constants
   
   fprintf('Value of c num = %d \n', ic) % First rate constant
   rho = cons_laws*x_0;
   rho(idx_c_high) = rho(idx_c_high)*(delta_j(ic)+1);
   ris = f_NLPC_restart(x_0, k_values, Sm, cons_laws, rho, idx_basic_species, ...
           vm, ind_one, max_counter, 0);
   selpart_x_eq_high(:, ic) = ris.x;

   clear rates ris

   rho = cons_laws*x_0;    % Second rate constant
   rho(idx_c_low) = rho(idx_c_low)*(delta_j(ic)+1);
   ris = f_NLPC_restart(x_0, k_values, Sm, cons_laws, rho, idx_basic_species, ...
           vm, ind_one, max_counter, 0);
   selpart_x_eq_low(:, ic) = ris.x;

clear rates ris

   rho = cons_laws*x_0;   % Third rate constant
   rho(idx_c_mean) = rho(idx_c_mean)*(delta_j(ic)+1);
   ris = f_NLPC_restart(x_0, k_values, Sm, cons_laws, rho, idx_basic_species, ...
           vm, ind_one, max_counter, 0);
   selpart_x_eq_mean(:, ic) = ris.x;

   clear rates ris
   
end

% Visualizzo il risultato.
aux_delta_x_eq_high = (selpart_x_eq_high - selpart_x_eq_high(:, 1)) ...
    ./ selpart_x_eq_high(:, 1);
partsel_delta_x_high = abs(aux_delta_x_eq_high(idx_proteins, :));

aux_delta_x_eq_low = (selpart_x_eq_low - selpart_x_eq_low(:, 1)) ...
    ./ selpart_x_eq_low(:, 1);
partsel_delta_x_low = abs(aux_delta_x_eq_low(idx_proteins, :));

aux_delta_x_eq_mean = (selpart_x_eq_mean - selpart_x_eq_mean(:, 1)) ...
    ./ selpart_x_eq_mean(:, 1);
partsel_delta_x_mean = abs(aux_delta_x_eq_mean(idx_proteins, :));

% Figure
figure('units','normalized','outerposition',[0 0 0.8 0.7]);
plot(delta_j, partsel_delta_x_high, 'k-*', 'Linewidth', 3)
hold on
plot(delta_j, partsel_delta_x_mean, 'b-*', 'Linewidth', 3)
hold on
plot(delta_j, partsel_delta_x_low, 'r-*', 'Linewidth', 3)
my_symlog('y', -4)
grid minor
ylim([0 4.05])
xlabel('$\hat{\Delta} c$', 'Interpreter', 'Latex')
ylabel('$\sqrt{Q(\hat{\Delta} c)}$', 'Interpreter', 'Latex')
set(gca, 'Fontsize', 15)
xticks(delta_j)
legend({'Maximum SSI';'Median SSI';'Minimum SSI'},'Location','NorthEastOutside')
