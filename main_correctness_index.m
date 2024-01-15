clc
clear
close all

% This code creates the figures representing the correctness of the SSI

set(0, 'defaultAxesTickLabelInterpreter','latex'); 
set(0, 'defaultLegendInterpreter','latex'); 
warning('off', 'all')

addpath('./funcs')

%% 1. Load structure containing information on species and reactions
file_crn = fullfile('data', 'CRC_CRN_nodrug_complete.mat');
folder_figures = './figures';
load(file_crn)

%% 2. Initial data
vm = CMIM.matrix.v;
n_species = numel(CMIM.species.names);
n_reactions = size(CMIM.matrix.S, 2);
n_cons_laws = size(CMIM.matrix.Nl, 1);
ind_one = CMIM.matrix.ind_one;

k_values = CMIM.rates.std_values;
Sm = CMIM.matrix.S;
x_0 = CMIM.species.std_initial_values;
idx_basic_species = find(x_0>0);
cons_laws = CMIM.matrix.Nl;
rho=cons_laws*x_0;

mut_lof=0;

max_counter = 300;

%% 3. Compute the values of the SSI
% 3.1. Compute equilibrium
ris_phys = f_NLPC_restart(x_0, k_values, Sm, cons_laws, rho, idx_basic_species, ...
          vm, ind_one, max_counter, 0);
x_values = ris_phys.x;

% 3.2. Compute the SSI
idx_sp=1:numel(x_values);
[SSI_k, SSI_c]=f_compute_SSI(idx_sp, x_values, k_values, ...
                                 Sm, cons_laws, rho, idx_basic_species, vm);

%% 4.  Computation of equilibrium for different choices of k and c
delta_j = [0 1 2 5 10 20];  
n_j = numel(delta_j);
max_counter = 300;

% Impact of increasing kinetic parameters referred to
% maximum
% minimum
% median
% SSI on global CRC-CRN
[max_e_k, idx_k_high] = max(SSI_k); 
[min_e_k, idx_k_low] = min(SSI_k);
aux=median(SSI_k); idx_k_mean=find(SSI_k==aux);

x_eq_high_k = zeros(n_species, n_j);
x_eq_low_k = zeros(n_species, n_j);
x_eq_mean_k = zeros(n_species, n_j);

ind_one=n_species+1;
for ik = 1:n_j % Compute the equilibrium by varying the selected rate constants
   
   fprintf('Value of k num = %d \n', ik) % First rate constant
   rates = CMIM.rates.std_values;
   rates(idx_k_high) = rates(idx_k_high)*(delta_j(ik)+1);  
   ris = f_NLPC_restart(x_0, rates, Sm, cons_laws, rho, idx_basic_species, ...
           vm, ind_one, max_counter, 0);
   x_eq_high_k(:, ik) = ris.x;

   clear rates ris

   rates = CMIM.rates.std_values;    % Second rate constant
   rates(idx_k_low) = rates(idx_k_low)*(delta_j(ik)+1);  
   ris = f_NLPC_restart(x_0, rates, Sm, cons_laws, rho, idx_basic_species, ...
           vm, ind_one, max_counter, 0);
   x_eq_low_k(:, ik) = ris.x;


   clear rates ris

   rates = CMIM.rates.std_values;    %Third rate constant
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

%%%%%

[max_e_c, idx_c_high] = max(SSI_c);
[min_e_c, idx_c_low] = min(SSI_c);
aux=median(SSI_c); idx_c_mean=find(SSI_c==aux);

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

delta_x_eq_high_c = (x_eq_high_c - x_eq_high_c(:, 1)) ./ x_eq_high_c(:, 1);
norm_delta_x_high_c = vecnorm(delta_x_eq_high_c, 2, 1);

delta_x_eq_low_c = (x_eq_low_c - x_eq_low_c(:, 1)) ./ x_eq_low_c(:, 1);
norm_delta_x_low_c = vecnorm(delta_x_eq_low_c, 2, 1);

delta_x_eq_mean_c = (x_eq_mean_c - x_eq_mean_c(:, 1)) ./ x_eq_mean_c(:, 1);
norm_delta_x_mean_c = vecnorm(delta_x_eq_mean_c, 2, 1);


%% Figure on sensitivity analysis of the CRC-CRN
fprintf('# Rate constant with SSI > 0.001  ---> %d \n', sum(SSI_k > 0.001))
fprintf('# Rate constant with SSI in [0.001, 10^-4]  ---> %d \n', ...
    sum((10^-4 < SSI_k & SSI_k< 0.001)))
fprintf('# Rate constant with SSI < 10^-4  ---> %d \n', sum(10^-4 > SSI_k))


fprintf('# c with SSI > 10^-2  ---> %d \n', sum(SSI_c > 10^(-2)))
fprintf('# c with SSI in [10^-3, 10^-2]  ---> %d \n', ...
    sum((10^-3 < SSI_c & SSI_c< 10^-2)))
fprintf('# c with SSI < 10^-3  ---> %d \n', sum(SSI_c < 10^(-3)))

figure_SSI_corrctness=figure('units','normalized','outerposition',[0 0  0.75 1]);

% 1. Histograms
% 1.1 rates
subplot(2,2,1)
binedges=[-Inf, -6:1:-1];
histogram(log10(SSI_k), 'BinEdges', binedges, 'Normalization', 'count')
hold on 
plot(log10(max_e_k), 1, 'k*', 'Markersize', 10, 'Linewidth', 10)
plot(log10(SSI_k(idx_k_mean)), 1, 'b*', 'Markersize', 10, 'Linewidth', 10)
plot(-7, 1, 'r*', 'Markersize', 10, 'Linewidth', 10)
xlim([-7,-1])
xticks(-6:1:-1)
xticklabels({'$10^{-6}$', '$10^{-5}$', '$10^{-4}$', '$10^{-3}$', '$10^{-2}$', '$10^{-1}$'})
set(gca, 'Fontsize', 18)
ylabel('Counts', 'Fontsize', 22)
xlabel('$e_j^\mathbf{k}$', 'Interpreter', 'Latex', 'Fontsize', 22)


% 1.2 conservation laws
subplot(2,2,2)
binedges=[-Inf,-3.5:0.5:-0.5];
histogram(log10(SSI_c), 'BinEdges', binedges, 'Normalization', 'count')
hold on 
plot(log10(max_e_c), 0, 'k*', 'Markersize', 10, 'Linewidth', 10)
plot(log10(SSI_c(idx_c_mean)), 0, 'b*', 'Markersize', 10, 'Linewidth', 10)
plot(log10(min_e_c), 0, 'r*', 'Markersize', 10, 'Linewidth', 10)
xlim([-3.5,-0.5])
xticks(-3:0.5:-0.5)
xticklabels({'$10^{-3}$', '$10^{-2.5}$', '$10^{-2}$', ...
    '$10^{-1.5}$', '$10^{-1}$', '$10^{-0.5}$'})
set(gca, 'Fontsize', 18)
ylabel('Counts', 'Fontsize', 22)
xlabel('$e_j^\mathbf{c}$', 'Interpreter', 'Latex', 'Fontsize', 22)


% Plots
% 2.1 rates
subplot(2, 2 ,3)
plot(delta_j, norm_delta_x_high_k, 'k-*', 'Linewidth', 3)
hold on
plot(delta_j, norm_delta_x_mean_k, 'b-*', 'Linewidth', 3)
hold on
plot(delta_j, norm_delta_x_low_k, 'r-*', 'Linewidth', 3)
set(gca, 'Fontsize', 18)
xlabel('$\hat{\Delta} k$', 'Interpreter', 'Latex', 'Fontsize', 22)
ylabel('$\sqrt{Q(\hat{\Delta} k)}$', 'Interpreter', 'Latex', 'Fontsize', 22)
xticks(delta_j);
my_symlog('y')
grid on

% 2.2 conservation laws
subplot(2,2, 4)
plot(delta_j, norm_delta_x_high_c, 'k-*', 'Linewidth', 3)
hold on
plot(delta_j, norm_delta_x_mean_c, 'b-*', 'Linewidth', 3)
hold on
plot(delta_j, norm_delta_x_low_c, 'r-*', 'Linewidth', 3)
my_symlog('y')
set(gca, 'Fontsize', 18)
xlabel('$\hat{\Delta} c$', 'Interpreter', 'Latex', 'Fontsize', 22)
ylabel('$\sqrt{Q(\hat{\Delta} c)}$', 'Interpreter', 'Latex', 'Fontsize', 22)
xticks(delta_j)
grid on
L=legend({'Max SSI';'Median SSI';'Min SSI'},'Location','NorthOutside', 'Orientation','horizontal');
L.Position=[0.5, 0.479, 0.001, 0.001];
L.FontSize = 20;

saveas(figure_SSI_corrctness, fullfile(folder_figures, 'figure2.png'))

%% Analysis restricted on TP53

%% Step 1. Compute SSI by considering only TP53's concentration
selected_proteins = {'TP53'};
[aux_, idx_proteins] = ismember(selected_proteins, CMIM.species.names);
[~, selpart_SSI]=f_compute_SSI(idx_proteins, x_values, k_values, ...
                                 Sm, cons_laws, rho, idx_basic_species, vm);

%% Step 2. Table on the conservation laws' constant                             
[temp_col, temp_raw] = find(cons_laws(:, idx_basic_species)');
bsp_names = CMIM.species.names(idx_basic_species(temp_col));
    % First we reorder the elemental species based on the conservation
    % laws.
[selpart_SSI_sorted, selpart_SSI_idx] = sort(selpart_SSI, 'descend');
bsp_names{selpart_SSI_idx};

aux = fullfile(folder_figures, 'SSI_TP53_c.txt');
fileID = fopen(aux, 'w');
n_selpart_SSI = numel(find(selpart_SSI >= 0.015));
for ii = 1:n_selpart_SSI
    fprintf(fileID, '%1.3f & \\verb| %s |  \\\\ \n', ...
        selpart_SSI_sorted(ii), ...
        bsp_names{selpart_SSI_idx(ii)});
end
fclose(fileID);   

%% Additional analysis
[max_e, idx_c_high] = max(selpart_SSI);
[min_e, idx_c_low] = min(selpart_SSI);
aux=median(selpart_SSI); idx_c_mean=find(selpart_SSI==aux);

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
legend({'Maximum SSI';'Median SSI';'Minimum SSI'},'Location','NorthOutside')
