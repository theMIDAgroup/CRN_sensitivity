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
protein_selected="Ras";
perc=0;
[~, ind_react_rem] = ...
    f_define_mutated_condition(protein_selected, ...
    x_0_combo, rate_constants_combo, CMIM_drug, perc);
react_rem=find(ind_react_rem==0);

null_species_names=["RP_GAP_Ras_GTP", "RP_G_GABP_GAP_Ras_GTP", ...
    "RP_ShP_G_GABP_GAP_Ras_GTP", "ERBP_GAP_Ras_GTP", ...
    "ERBP_G_GABP_GAP_Ras_GTP", "ERBP_ShP_G_GABP_GAP_Ras_GTP", ...
    "ERB3P_GAP_Ras_GTP", ...
    "ERB3P_G_GABP_GAP_Ras_GTP", "ERB3P_ShP_G_GABP_GAP_Ras_GTP"];

[~, null_species]=ismember(null_species_names, CMIM_drug.species.names);

S_drug(:,react_rem)=[];
k_drug(react_rem)=[];
n_react_drug=n_reactions-numel(react_rem);

% check react_rem is in order in v!
%[~,react_to_rem]=ismember(react_rem,v_drug(:,1));

v_drug(react_rem,:)=[];
first_column=v_drug(:,1);
first_reagent=v_drug(:,2);
second_reagent=v_drug(:,3);
idx_reactions=setdiff(1:n_reactions, react_rem);

for ir=1:n_react_drug
    disp('idx old')
    disp(CMIM_drug.matrix.v(ir,1))
    disp(v_drug(ir,1))
    disp(find(v_drug(ir,1)==idx_reactions))
    first_column(ir)=find(v_drug(ir,1)==idx_reactions);
end

v_drug(:,1)=first_column;

% 3.3. Compute equilibrium point of the new MIM
rho_drug=cons_laws_drug*x_0_combo;
ris=f_NLPC_restart(x_0_combo, k_drug, S_drug, cons_laws_drug, rho_drug,...
    idx_basic_species_drug, v_drug, ind_one_drug, max_counter, 1);
x_e_drug_2=ris.x;

%% Step 4. SSI computation

x_e_drug_2(null_species)=0;
n_basic_species_drug=numel(idx_basic_species_drug);

% 1. Compute the analic jacobian of v
jacobian_v = f_compute_analytic_jacobian_v(v_drug, n_species_drug, ind_one_drug);

% 2. Evaluate the Jacobian in x of f
eval_jac = f_evaluate_jacobian_neworder(k_drug, x_e_drug_2, ...
    S_drug, idx_basic_species_drug, jacobian_v, cons_laws_drug);

eval_jac_k=f_evaluate_jacobian_k(v_drug, x_e_drug_2, S_drug, idx_basic_species_drug);
eval_jac_c=f_evaluate_jacobian_c(n_species+n_new_species, n_basic_species_drug, 0);

% 3. Calcolate Jac_k_xeq e Jac_c_eq
Jac_k_xeq= -inv(eval_jac)*eval_jac_k;
Jac_c_xeq= -inv(eval_jac)*eval_jac_c;

x_values=x_e_drug_2;
x_values(null_species)=1;
L_k=(1./x_values).*(Jac_k_xeq).*(k_drug');
L_c=(1./x_values).*(Jac_c_xeq).*(rho_drug');

[U_k,D_k,V_k]=svd(L_k);
e_PCA_k =  ((V_k(:,1:n_species_drug).^2)*(diag(D_k).^2))/sum(diag(D_k).^2); %se presi con il quadrato: e_j=(L^T*L)_{jj}/sum_k(lambda_k)
[U_c,D_c,V_c]=svd(L_c);
e_PCA_c =((V_c(:,1:n_basic_species_drug).^2)*(diag(D_c).^2))/sum(diag(D_c).^2); %se presi con il quadrato: e_j=(L^T*L)_{jj}/sum_k(lambda_k)

%% Step 5. Figure
f_ssi_phys_vs_drug = figure('units','normalized','outerposition',[0 0  0.75 1]);

% compute physiological SSI
x_values=x_eq_phys;

jacobian_v_phys = f_compute_analytic_jacobian_v(vm, n_species, ind_one);

eval_jac_phys = f_evaluate_jacobian_neworder(rate_constants_phys, x_values, ...
    Sm, idx_basic_species, jacobian_v_phys, cons_laws);

eval_jac_k_phys=f_evaluate_jacobian_k(vm, x_values, Sm, idx_basic_species);
eval_jac_c_phys=f_evaluate_jacobian_c(n_species, n_cons_laws,0);

Jac_k_xeq_phys= -inv(eval_jac_phys)*eval_jac_k_phys;
Jac_c_xeq_phys= -inv(eval_jac_phys)*eval_jac_c_phys;

L_c_phys=(1./x_values).*(Jac_c_xeq_phys).*(rho');
[U_c_phys,D_c_phys,V_c_phys]=svd(L_c_phys);
e_PCA_c_phys =  ((V_c_phys(:,1:n_cons_laws).^2)*(diag(D_c_phys).^2))/sum(diag(D_c_phys).^2); %se presi con il quadrato: e_j=(L^T*L)_{jj}/sum_k(lambda_k)

L_k_phys=(1./x_values).*(Jac_k_xeq_phys).*(rate_constants_phys');

[U_k_phys,D_k_phys,V_k_phys]=svd(L_k_phys);
e_PCA_k_phys =  ((V_k_phys(:,1:n_species).^2)*(diag(D_k_phys).^2))/sum(diag(D_k_phys).^2); %se presi con il quadrato: e_j=(L^T*L)_{jj}/sum_k(lambda_k)

e_PCA_k_phys(react_rem)=[];

% figures

[aux_PCA_k_phys_sort, order_phys]=sort(e_PCA_k_phys, 'descend');

subplot(1,2,1)
semilogy(e_PCA_k(order_phys), 'k', 'linewidth', 3)
hold on
semilogy(aux_PCA_k_phys_sort, 'r--', 'linewidth', 3)
ylim([10^-16, 10^0])

[aux_PCA_c_phys_sort, order_phys]=sort(e_PCA_c_phys, 'descend');
subplot(1,2,2)
semilogy(e_PCA_c(order_phys), 'k', 'linewidth', 3)
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

[max_values_SSI_k, ord_max_SSI_k]=sort(e_PCA_k,'descend');
[max_values_SSI_c, ord_max_SSI_c]=sort(e_PCA_c,'descend');
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

selpart_L = L_c(idx_proteins, :);
[U, Sigma, V]=svd(selpart_L);
selection=size(selpart_L,1);
eing=diag(Sigma(1:selection, 1:selection));
sum_eig=sum(eing.^2);
selpart_e_PCA_c=((V(:,1:selection).^2)*(eing.^2))/sum_eig;

[max_selpart, order_selpart]=sort(selpart_e_PCA_c, 'descend');

selpart_L = L_k(idx_proteins, :);
[U, Sigma, V]=svd(selpart_L);
selection=size(selpart_L,1);
eing=diag(Sigma(1:selection, 1:selection));
sum_eig=sum(eing.^2);
selpart_e_PCA_k=((V(:,1:selection).^2)*(eing.^2))/sum_eig;
