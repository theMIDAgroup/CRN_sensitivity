clc
clear
close all

%% This code computes the SSIs for a CRC-CRN when drugs DBF and TMT are inserted

set(0, 'defaultAxesTickLabelInterpreter','latex');
set(0, 'defaultLegendInterpreter','latex');
warning('off', 'all')

addpath('./funcs')

%% Step 1. Define general parameters
% 1.1. Data
target_folder = 'data';
file_mim_clean = fullfile(target_folder, 'CRC_CRN_nodrug_complete.mat');

% 1.2. Folders and files
folder_results = 'results';
file_ris_phys = fullfile(folder_results, 'x_phys.mat');

folder_results_drugs = 'results/drugs';
file_ris_drugs = fullfile(folder_results_drugs, 'nlpc_DBF_TMT_on_mut_Ras_50.00_240.00.mat');

% 1.3. Starting mutation and drugs
mut_prot = 'Ras';
drug1 = 'DBF'; drug2 = 'TMT';
drug = strcat(drug1, '_', drug2);
mut_lof=0;
perc = 0;

%% Step 2. Load and store data
load(file_mim_clean)

% physiological case
rate_constants_phys = CMIM.rates.std_values;
x_0_phys = CMIM.species.std_initial_values;
idx_basic_species = find(x_0_phys);
cons_laws = CMIM.matrix.Nl;
vm = CMIM.matrix.v;
Sm = CMIM.matrix.S;

n_species = numel(x_0_phys);
n_cons_laws = size(cons_laws, 1);
ind_one = CMIM.matrix.ind_one;

% max counter for NLPC
max_counter = 500;

% drug initial concentration
init_drug1 = 50;
init_drug2 = 240;

%% Step 2. Compute physiological equilibrium
rho=cons_laws*x_0_phys;

ris_phys = f_NLPC_restart(x_0_phys, rate_constants_phys, Sm, cons_laws, rho, idx_basic_species, ...
           vm, ind_one, max_counter, 0);
x_eq_phys=ris_phys.x;

%% Step 3. Simulate mutation + drugs
% 3.1. Add drugs to CRN

[CMIM_drug, n_new_species1] = f_add_drug_Raf_from_file(CMIM, drug1);
[CMIM_drug, n_new_species2] = f_add_drug_Raf_from_file(CMIM_drug, drug2);
n_new_species = n_new_species1 + n_new_species2;

[~, idx_k1] = ismember('cd_1', CMIM_drug.rates.names); [~, idx_k2] = ismember('cd_2', CMIM_drug.rates.names);
[~, idx_k3] = ismember('cd_3', CMIM_drug.rates.names); [~, idx_k4] = ismember('cd_4', CMIM_drug.rates.names);
[~, idx_k5] = ismember('cd_5', CMIM_drug.rates.names); [~, idx_k6] = ismember('cd_6', CMIM_drug.rates.names);
[~, idx_k7] = ismember('cd_7', CMIM_drug.rates.names); [~, idx_k8] = ismember('cd_8', CMIM_drug.rates.names);
idx_k = [idx_k1 idx_k2 idx_k3 idx_k4 idx_k5 idx_k6 idx_k7 idx_k8];

k1_drug = 0.106 * 1e-3; k2_drug = 0.593 * 1e-4; k3_drug = k1_drug; k4_drug = 1.2296 * 1e-3;
k5_drug = k1_drug; k6_drug = k4_drug; k7_drug = 0.1 * 1e-1; k8_drug = 0.33 * 1e-2;
k = [k1_drug k2_drug k3_drug k4_drug k5_drug k6_drug k7_drug k8_drug];

[~, idx_d1] = ismember(drug1, CMIM_drug.species.names);
[~, idx_d2] = ismember(drug2, CMIM_drug.species.names);
idx_basic_species_drug = [idx_basic_species; idx_d1; idx_d2];

x_0_combo = padarray(x_eq_phys,[n_new_species 0],0,'post');
x_0_combo(idx_d1) = init_drug1; x_0_combo(idx_d2) = init_drug2;
rate_constants_combo = rate_constants_phys; rate_constants_combo(idx_k) = k;

ind_one_drug = CMIM_drug.matrix.ind_one;

cons_laws_drug=CMIM_drug.matrix.Nl;
n_basic_drug=size(cons_laws_drug,1);
n_species_drug=size(cons_laws_drug,2);

S_drug=CMIM_drug.matrix.S;
k_drug=rate_constants_combo;
n_reactions=numel(k_drug);
v_drug=CMIM_drug.matrix.v;

% 3.2. Add RAS mutation
% deleting columns from S and entries from v and k refereed to all the reactions ivolved into RAS mutation

MIM_mut=f_compute_eq_mutated_CRN("Ras", CMIM_drug, idx_basic_species_drug, x_0_combo, rate_constants_combo);
x_e_drug_2=MIM_mut.species.x_eq;
null_species=MIM_mut.species.null_species;

%% Step 4. SSI computation

x_e_drug_2(null_species)=0;
idx_sp=find(x_e_drug_2);

k_mut_drug=MIM_mut.rates.std_values;
S_mut_drug=MIM_mut.matrix.S;
cons_laws_mut_drug=MIM_mut.matrix.Nl;
idx_basic_species_mut_drug=MIM_mut.species.idx_basic_species;
v_mut_drug=MIM_mut.matrix.v;
rho_mut_drug=cons_laws_mut_drug*x_0_combo;

[SSI_k, SSI_c]=f_compute_SSI(idx_sp, x_e_drug_2, k_mut_drug,...
                                    S_mut_drug, cons_laws_mut_drug, rho_mut_drug, ...
                                    idx_basic_species_mut_drug, v_mut_drug);

%% Step 5. Figure
f_ssi_phys_vs_drug = figure('units','normalized','outerposition',[0 0  0.75 1]);

% compute physiological SSI
react_rem=MIM_mut.info.react_rem;

x_values=x_eq_phys;
idx_sp_phys=1:numel(x_values);
rho_phys=cons_laws*x_0_phys;

[SSI_k_phys, SSI_c_phys]=f_compute_SSI(idx_sp_phys, x_values, rate_constants_phys, ...
                                           Sm, cons_laws, rho_phys, idx_basic_species, vm);

SSI_k_phys(react_rem)=[];

% figures

[aux_PCA_k_phys_sort, order_phys]=sort(SSI_k_phys, 'descend');

subplot(1,2,1)
semilogy(SSI_k(order_phys), 'k', 'linewidth', 3)
hold on
semilogy(aux_PCA_k_phys_sort, 'r--', 'linewidth', 3)
ylim([10^-16, 10^0])

[aux_PCA_c_phys_sort, order_phys]=sort(SSI_c_phys, 'descend');
subplot(1,2,2)
semilogy(SSI_c(order_phys), 'k', 'linewidth', 3)
hold on
semilogy([aux_PCA_c_phys_sort;0;0], 'r--', 'linewidth', 3)

%% tables

arrow=CMIM_drug.reactions.arrow;
react2flux=CMIM_drug.reactions.reactions2flux_rates;
arrow(react2flux(react_rem))=[];
react2flux(react_rem)=[];
details=CMIM_drug.reactions.details;
details(react_rem,:)=[];
for i=1:numel(react2flux)
    sumflux=sum(react2flux(i)>react_rem);
    if(sumflux>0)
        react2flux(i)=react2flux(i)-sumflux;
    end
end

[max_values_SSI_k, ord_max_SSI_k]=sort(SSI_k,'descend');
[max_values_SSI_c, ord_max_SSI_c]=sort(SSI_c,'descend');
str=strings(10,2);
for i=1:10
    str(i,1)=convertCharsToStrings(arrow{react2flux(ord_max_SSI_k(i))});
    str(i,2)=convertCharsToStrings(details{ord_max_SSI_k(i),2});
end
T_k=table(str, max_values_SSI_k(1:10), 'VariableNames', {'Reaction List', 'SSI'});

str=strings(10,1);
for i=1:10
    str(i,1)=convertCharsToStrings(CMIM_drug.species.names{idx_basic_species_drug(ord_max_SSI_c(i))});
end
T_c=table(str, max_values_SSI_c(1:10), 'VariableNames', {'Basic Species List', 'SSI'});

%% ERKPP
selected_proteins = {'ERKPP'};
[aux_, idx_proteins] = ismember(selected_proteins, CMIM_drug.species.names);


[selpart_SSI_k, selpart_SSI_c]=f_compute_SSI(idx_proteins,  x_e_drug_2, k_mut_drug,...
                                    S_mut_drug, cons_laws_mut_drug, rho_mut_drug, ...
                                    idx_basic_species_mut_drug, v_mut_drug);
