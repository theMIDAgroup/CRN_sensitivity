clc
clear
close all
% This code creates the figures representing the correctness of the SSI

set(0, 'defaultAxesTickLabelInterpreter','latex'); 
set(0, 'defaultLegendInterpreter','latex'); 

%% Step 1. Load structure containing information on species and reactions
file_crn = fullfile('data', 'CRC_CRN_nodrug.mat');
load(file_crn)
file_eq = fullfile('data', 'results_physiological.mat');
load(file_eq, 'ris_phys')
% file_crc_crn = fullfile('data', 'CRC_CRN.mat');
% load(file_crc_crn)

%% Step 2. Import the initial data

vm = new_CMIM.matrix.v;
n_species = numel(new_CMIM.species.names);
n_reactions = size(new_CMIM.matrix.S, 2);
n_cons_laws = size(new_CMIM.matrix.Nl, 1);
idx_one = new_CMIM.matrix.ind_one;


k_values = new_CMIM.rates.std_values;
x_values=ris_phys.x_eq(1:end-2);
Sm = new_CMIM.matrix.S;
x_0 = new_CMIM.species.std_initial_values;
idx_basic_species = find(x_0>0);
cons_laws = new_CMIM.matrix.Nl;

%% Step 3. Compute the sensitive matrix
% Compute the analytic jacobian of v 
jacobian_v = f_compute_analytic_jacobian_v(vm, n_species, idx_one);

% Calculate the Jacobian of f with respect to x
eval_jac = f_evaluate_jacobian_neworder(k_values, x_values, ...
    Sm, idx_basic_species, jacobian_v, cons_laws);

% Calculate the Jacobian of f with respect to k
eval_jac_k=f_evaluate_jacobian_k(vm, x_values, Sm, idx_basic_species);

% Calculate the k sensitivity matrix
Jac_k_xeq= -inv(eval_jac)*eval_jac_k;

% Calculate Jacobian of f with respect to c
eval_jac_c=f_evaluate_jacobian_c(n_species, n_cons_laws,0);

% Calcolate Jac_k xeq
Jac_c_xeq= -inv(eval_jac)*eval_jac_c;

%% Calculate the c sensitivity matrix and indices 
% Calculate the normalized k sensitivity matrix
L_k=(1./x_values).*(Jac_k_xeq).*(k_values');

% Calculate the SSI for k
[U_k,D_k,V_k]=svd(L_k);
e_PCA_k =  ((V_k(:,1:n_species).^2)*diag(D_k.^2))/sum(diag(D_k.^2)); %se presi con il quadrato: e_j=(L^T*L)_{jj}/sum_k(lambda_k)

% Calculate the normalized c sensitivity matrix
c=cons_laws*x_0;
L_c=(1./x_values).*(Jac_c_xeq).*(c');

% Calculate the SSI for c

[U_c,D_c,V_c]=svd(L_c);
e_PCA_c =  ((V_c(:,1:n_cons_laws).^2)*(diag(D_c).^2))/sum(diag(D_c).^2); %se presi con il quadrato: e_j=(L^T*L)_{jj}/sum_k(lambda_k)

%% Graphs: 
% Fig 1

delta_j = [0 1 2 5 10 20];  
n_j = numel(delta_j);
% Impact of increasing kinetic parameters referred to
% maximum
% minimum
% median
% SSI on global CRC-CRN

[max_e_k, idx_k_high] = max(e_PCA_k); 
[min_e_k, idx_k_low] = min(e_PCA_k);
aux=median(e_PCA_k);idx_k_mean=find(e_PCA_k==aux);

rho = cons_laws*x_0;

x_eq_high_k = zeros(n_species, n_j);
x_eq_low_k = zeros(n_species, n_j);
x_eq_mean_k = zeros(n_species, n_j);

max_counter = 300;

ind_one=n_species+1;
for ik = 1:n_j % Calcolo il punto di equilibrio facendo variare le rate constant selezionate
   
   fprintf('Value of k num = %d \n', ik) % Prima rate constant
   rates = new_CMIM.rates.std_values;
   rates(idx_k_high) = rates(idx_k_high)*(delta_j(ik)+1);
   
   ris = f_PNG_restart(x_0, rates, Sm, cons_laws, rho, idx_basic_species, ...
           vm, ind_one, max_counter);
   x_eq_high_k(:, ik) = ris.x;

   clear rates ris

   rates = new_CMIM.rates.std_values;    % Seconda rate constant
   rates(idx_k_low) = rates(idx_k_low)*(delta_j(ik)+1);
   
   ris = f_PNG_restart(x_0, rates, Sm, cons_laws, rho, idx_basic_species, ...
           vm, ind_one, max_counter);
   x_eq_low_k(:, ik) = ris.x;


   clear rates ris

   rates = new_CMIM.rates.std_values;    % Seconda rate constant
   rates(idx_k_mean) = rates(idx_k_mean)*(delta_j(ik)+1);
   
   ris = f_PNG_restart(x_0, rates, Sm, cons_laws, rho, idx_basic_species, ...
           vm, ind_one, max_counter);
   x_eq_mean_k(:, ik) = ris.x;

   clear rates ris
   

end

delta_x_eq_high_k = (x_eq_high_k - x_eq_high_k(:, 1)) ./ x_eq_high_k(:, 1);
norm_delta_x_high_k = vecnorm(delta_x_eq_high_k, 2, 1);

delta_x_eq_low_k = (x_eq_low_k - x_eq_low_k(:, 1)) ./ x_eq_low_k(:, 1);
norm_delta_x_low_k = vecnorm(delta_x_eq_low_k, 2, 1);

delta_x_eq_mean_k = (x_eq_mean_k - x_eq_mean_k(:, 1)) ./ x_eq_mean_k(:, 1);
norm_delta_x_mean_k = vecnorm(delta_x_eq_mean_k, 2, 1);

% Impact of increasing conservation laws' constants referred to
% maximum
% minimum
% median
% SSI on global CRC-CRN

[max_e_c, idx_c_high] = max(e_PCA_c);
[min_e_c, idx_c_low] = min(e_PCA_c);
aux=median(e_PCA_c);
idx_c_mean=find(e_PCA_c==aux);

x_eq_high_c = zeros(n_species, n_j);
x_eq_low_c = zeros(n_species, n_j);
x_eq_mean_c = zeros(n_species, n_j);

max_counter = 300;

for ic = 1:n_j % Calcolo il punto di equilibrio facendo variare le rate constant selezionate
   
   fprintf('Value of c num = %d \n', ic) % Prima rate constant
   rho = cons_laws*x_0;
   rho(idx_c_high) = rho(idx_c_high)*(delta_j(ic)+1);
   ris = f_PNG_restart(x_0, k_values, Sm, cons_laws, rho, idx_basic_species, ...
           vm, ind_one, max_counter);
   x_eq_high_c(:, ic) = ris.x;

   clear rates ris

   rho = cons_laws*x_0  ;  % Seconda rate constant
   rho(idx_c_low) = rho(idx_c_low)*(delta_j(ic)+1);
   ris = f_PNG_restart(x_0, k_values, Sm, cons_laws, rho, idx_basic_species, ...
           vm, ind_one, max_counter);
   x_eq_low_c(:, ic) = ris.x;


clear rates ris

   rho = cons_laws*x_0;    % Seconda rate constant
   rho(idx_c_mean) = rho(idx_c_mean)*(delta_j(ic)+1);
   ris = f_PNG_restart(x_0, k_values, Sm, cons_laws, rho, idx_basic_species, ...
           vm, ind_one, max_counter);
   x_eq_mean_c(:, ic) = ris.x;

   clear rates ris
   

end

delta_x_eq_high = (x_eq_high_c - x_eq_high_c(:, 1)) ./ x_eq_high_c(:, 1);
norm_delta_x_high = vecnorm(delta_x_eq_high, 2, 1);

delta_x_eq_low = (x_eq_low_c - x_eq_low_c(:, 1)) ./ x_eq_low_c(:, 1);
norm_delta_x_low = vecnorm(delta_x_eq_low, 2, 1);

delta_x_eq_mean = (x_eq_mean_c - x_eq_mean_c(:, 1)) ./ x_eq_mean_c(:, 1);
norm_delta_x_mean = vecnorm(delta_x_eq_mean, 2, 1);

% figures
figure('units','normalized','outerposition',[0 0 1 0.7]);
subplot(1, 3 ,1)
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

subplot(1,3,2)
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
legend({'Maximum SSI';'Median SSI';'Minimum SSI'}, 'Location', 'eastout')

% Fig 2
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


[max_e, idx_c_high] = max(selpart_e_PCA); % Prima rate constant
[min_e, idx_c_low] = min(selpart_e_PCA);     % Seconda rate constant

aux=median(selpart_e_PCA); %terza rate costant intermedia
idx_c_mean=find(selpart_e_PCA==aux);

selpart_x_eq_high = zeros(n_species, n_j);
selpart_x_eq_low = zeros(n_species, n_j);
selpart_x_eq_mean=zeros(n_species, n_j);

for ic = 1:n_j % Calcolo il punto di equilibrio facendo variare le rate constant selezionate
   
   fprintf('Value of c num = %d \n', ic) % Prima rate constant
   rho = cons_laws*x_0;
   rho(idx_c_high) = rho(idx_c_high)*(delta_j(ic)+1);
   ris = f_PNG_restart(x_0, k_values, Sm, cons_laws, rho, idx_basic_species, ...
           vm, ind_one, max_counter);
   selpart_x_eq_high(:, ic) = ris.x;

   clear rates ris

   rho = cons_laws*x_0;    % Seconda rate constant
   rho(idx_c_low) = rho(idx_c_low)*(delta_j(ic)+1);
   ris = f_PNG_restart(x_0, k_values, Sm, cons_laws, rho, idx_basic_species, ...
           vm, ind_one, max_counter);
   selpart_x_eq_low(:, ic) = ris.x;

clear rates ris

   rho = cons_laws*x_0;   % Seconda rate constant
   rho(idx_c_mean) = rho(idx_c_mean)*(delta_j(ic)+1);
   ris = f_PNG_restart(x_0, k_values, Sm, cons_laws, rho, idx_basic_species, ...
           vm, ind_one, max_counter);
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

% figure
figure('units','normalized','outerposition',[0 0 0.8 0.7]);
plot(delta_j, partsel_delta_x_high, 'k-*', 'Linewidth', 3)
hold on
plot(delta_j, partsel_delta_x_mean, 'b-*', 'Linewidth', 3)
hold on
plot(delta_j, partsel_delta_x_low, 'r-*', 'Linewidth', 3)
my_symlog('y')
% set(gca, 'Yscale', 'log')
grid minor
ylim([0 4.05])
xlabel('$\hat{\Delta} c$', 'Interpreter', 'Latex')
ylabel('$\sqrt{Q(\hat{\Delta} c)}$', 'Interpreter', 'Latex')
set(gca, 'Fontsize', 15)
xticks(delta_j)
legend({'Maximum SSI';'Median SSI';'Minimum SSI'},'Location','northeastout')
