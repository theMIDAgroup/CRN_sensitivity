clc
clear
close all

%% FAI INDICI TOTALI!!!

% This code computes the SSI for a CRC-CRN when drugs DBF and TMT are inserted

set(0, 'defaultAxesTickLabelInterpreter','latex');
set(0, 'defaultLegendInterpreter','latex');
warning('off', 'all')

addpath('./funcs')

%% Step 1. Define general parameters
% 1.1. Data
target_folder = 'data';
file_mim_clean = fullfile(target_folder, 'CRC_CRN_nodrug_complete.mat');

% 1.2. Folders and files
folder_results = 'results/drugs';

if exist(folder_results, 'dir')==0
    mkdir(folder_results)
end

% 1.3. Starting mutation
mut_prot = 'Ras';
drug2 = 'TMT';
init_drug2 = 240;
perc = 0;

% 1.4 Fix a toll
toll=10^-14; % change it

% 1.5 Set NLPC parameter
max_counter = 500;
%% Step 2. Load and store data
load(file_mim_clean)

% physiological case
x_0_phys = CMIM.species.std_initial_values;
idx_basic_species = find(x_0_phys);
rate_constants_phys = CMIM.rates.std_values;
vm = CMIM.matrix.v;
cons_laws = CMIM.matrix.Nl;
Sm = CMIM.matrix.S;
n_species = numel(x_0_phys);
n_cons_laws = size(cons_laws, 1);
ind_one = CMIM.matrix.ind_one;


%% Step 2. Compute physiological equilibrium
rho=cons_laws*x_0_phys;

ris_phys = f_NLPC_restart(x_0_phys, rate_constants_phys, Sm, cons_laws, rho, idx_basic_species, ...
           vm, ind_one, max_counter, 0);

x_eq_phys = ris_phys.x;
%% Step 3. Simulate mutation + drugs
% 3.1. Add drugs to CRN

[CMIM_drug, n_new_species] = f_add_drug_Raf_from_file(CMIM, drug2);

[~, idx_k1] = ismember('cd_1', CMIM_drug.rates.names); [~, idx_k2] = ismember('cd_2', CMIM_drug.rates.names);
[~, idx_k3] = ismember('cd_3', CMIM_drug.rates.names); [~, idx_k4] = ismember('cd_4', CMIM_drug.rates.names);
[~, idx_k5] = ismember('cd_5', CMIM_drug.rates.names); [~, idx_k6] = ismember('cd_6', CMIM_drug.rates.names);
idx_k = [idx_k1 idx_k2 idx_k3 idx_k4 idx_k5 idx_k6];

k1_drug = 0.106 * 1e-3; k2_drug = 1.2296 * 1e-3;
k3_drug = k1_drug; k4_drug = 1.2296 * 1e-3;
k5_drug = 0.1 * 1e-1; k6_drug = 0.33 * 1e-2;
k = [k1_drug k2_drug k3_drug k4_drug k5_drug k6_drug];

[~, idx_d2] = ismember(drug2, CMIM_drug.species.names);
idx_basic_species_drug = [idx_basic_species; idx_d2];

x_0_drug = [x_eq_phys; zeros(n_new_species,1)];
x_0_drug(idx_d2) = init_drug2;
rate_constants_drug = rate_constants_phys;
rate_constants_drug(idx_k) = k;

ind_one_drug = CMIM_drug.matrix.ind_one;

cons_laws_drug=CMIM_drug.matrix.Nl;
n_basic_drug=size(cons_laws_drug,1);
n_species_drug=size(cons_laws_drug,2);

S_drug=CMIM_drug.matrix.S;
k_drug=rate_constants_drug;
n_reactions=numel(k_drug);
v_drug=CMIM_drug.matrix.v;


% 3.2. Add RAS mutation

MIM_mut=f_compute_eq_mutated_CRN("Ras", CMIM_drug, idx_basic_species_drug, x_0_drug, rate_constants_drug);
react_rem=MIM_mut.info.react_rem;

%% Step 4. SSI computation
x_e_drug_2=MIM_mut.species.x_eq;
idx_sp=find(x_e_drug_2>10^-14);

k_mut_drug=MIM_mut.rates.std_values;
S_mut_drug=MIM_mut.matrix.S;
cons_laws_mut_drug=MIM_mut.matrix.Nl;
idx_basic_species_mut_drug=MIM_mut.species.idx_basic_species;
v_mut_drug=MIM_mut.matrix.v;
rho_mut_drug=cons_laws_mut_drug*x_0_drug;

SSI=f_compute_SSI_tot(idx_sp, x_e_drug_2, k_mut_drug,...
                                    S_mut_drug, cons_laws_mut_drug, rho_mut_drug,...
                                    idx_basic_species_mut_drug, v_mut_drug);

[value_tot, ind_tot]=sort(SSI, 'descend');

% ERKPP
 selected_proteins = {'ERKPP'};
[aux_, idx_proteins] = ismember(selected_proteins, CMIM_drug.species.names);

selpart_SSI=f_compute_SSI_tot(idx_proteins,  x_e_drug_2, k_mut_drug,...
                                    S_mut_drug, cons_laws_mut_drug, rho_mut_drug, idx_basic_species_drug, v_mut_drug);

%% Step 5. Tables

kk=append('R', string(1:numel(k_mut_drug)+numel(react_rem)));
cc=append('CL ', string(1:numel(idx_basic_species_drug)));

param=[(kk(end-5:end))'; cc(end)'];
SSI_param=[SSI(numel(k_mut_drug)-5:numel(k_mut_drug)); SSI(end)];


[value_SSI, ind_ord_SSI]=sort(SSI_param, 'descend');
[a, new_ind]=ismember([numel(k_mut_drug)-5:numel(k_drug), numel(SSI)], ind_tot);


SSI_param=SSI_param(ind_ord_SSI, :);
param_ord=param(ind_ord_SSI);

table(param_ord, SSI_param, 'VariableNames',{'Kinetic Parameter h', 'e_j^h'})

table_file=fullfile(folder_results, 'TMT_SSI.txt');
fileID = fopen(table_file, 'w');
disp(['Writing on ', table_file, '...'])

row_table=numel(param_ord);
for ii = 1:row_table
    fprintf(fileID, '%s &  %1.2e  \\\\ \\hline \n', ...
        string(param_ord(ii)), string(SSI_param(ii)));
end
fclose(fileID);

% ERKPP
selpart_SSI_param=[selpart_SSI(numel(k_mut_drug)-5:numel(k_mut_drug)); selpart_SSI(end)];

[value_SSI, selpart_ind_ord_SSI]=sort(selpart_SSI_param, 'descend');
selpart_SSI_param=selpart_SSI_param(selpart_ind_ord_SSI, :);
selpart_param_ord=param(selpart_ind_ord_SSI);

table(selpart_param_ord, selpart_SSI_param, 'VariableNames',{'Kinetic Parameter h', 'e_j^h'})

table_file=fullfile(folder_results, 'TMT_selpart_SSI.txt');
fileID = fopen(table_file, 'w');
disp(['Writing on ', table_file, '...'])

row_table=numel(selpart_param_ord);
for ii = 1:row_table
    fprintf(fileID, '%s &  %1.2e  \\\\ \\hline \n', ...
        string(selpart_param_ord(ii)), string(selpart_SSI_param(ii)));
end
fclose(fileID);

%% GRAFICO
% indice di sensitività in ordine 
figure
semilogy(new_ind, value_tot(new_ind), '*')
hold on
semilogy(value_tot, 'k', 'linewidth', 1)

