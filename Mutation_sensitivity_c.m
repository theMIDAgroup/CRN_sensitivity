%% Clear all
clc
clear
close all

addpath(fullfile('.', 'funcs'))
folder_results = fullfile('.', 'results');
if ~exist(folder_results, 'dir')
   mkdir(folder_results)
end
warning('off', 'all')

%% Step 1. Load data
% Load structure containing information on species and reactions
file_crn = fullfile('data', 'CRC_CRN_nodrug.mat'); % Da cambiare in base a dove metti i file
load(file_crn, 'new_CMIM')
file_crc_crn = fullfile('data', 'CRC_CRN.mat');
load(file_crc_crn, 'CMIM')
list_reactions = CMIM.reactions; clear CRC
% file_eq_phys = fullfile('data', 'results_physiological.mat');
% load(file_eq_phys, 'ris_phys')
% x_eq_phys = ris_phys.x_eq(1:end-2); clear ris_phys
% Load equilibrium for the physiological network computed as in Sommariva
% et al, Scientific Reports, 2021.

%% Step 2. Import the initial data
% number of species, reactions, conservation laws
n_species = numel(new_CMIM.species.names);
n_reactions = size(new_CMIM.matrix.S, 2);
n_cons_laws = size(new_CMIM.matrix.Nl, 1);


% define the new matrices and vectors of the system
Sm = new_CMIM.matrix.S;
vm = new_CMIM.matrix.v;
k_values=new_CMIM.rates.std_values;
cons_laws=new_CMIM.matrix.Nl;
x_0_phys=new_CMIM.species.std_initial_values;
idx_basic_species = find(x_0_phys>0);
max_counter=300;
idx_one=new_CMIM.matrix.ind_one;
% ris_phys=f_NLPC_restart(x_0_phys, k_values, Sm, cons_laws, cons_laws*x_0_phys,...
%         idx_basic_species, vm,idx_one, max_counter, 0);
% x_eq_phys = ris_phys.x;

%% Step 3. Create the mutatated CRC
lof_mutation = {'APC', 'SMAD4'};
gof_mutation = {'Ras'};
lof_mutation_type2 = {'TP53'};
all_proteins = [gof_mutation, lof_mutation, lof_mutation_type2];

protein=["Ras", "APC", "SMAD4", "TP53"];
f_eff_mut_log_c = figure('units','normalized','outerposition',[0 0 1 1]);

for im=1:numel(protein)
    protein_selected=protein(im);
    

    switch protein_selected
        case 'APC'
            mut_lof=1;
            null_species=[];
            const_species=[];
        case 'SMAD4'
            mut_lof=1;
            null_species=[];
            const_species='TFBIS';
        case 'TP53'
            mut_lof=1;
            null_species=[];
            const_species=[];
        case 'Ras'
            mut_lof=0;
            null_species_names=["RP_GAP_Ras_GTP", "RP_G_GABP_GAP_Ras_GTP", ...
                "RP_ShP_G_GABP_GAP_Ras_GTP", "ERBP_GAP_Ras_GTP", ...
                "ERBP_G_GABP_GAP_Ras_GTP", "ERBP_ShP_G_GABP_GAP_Ras_GTP", ...
                "ERB3P_GAP_Ras_GTP", ...
                "ERB3P_G_GABP_GAP_Ras_GTP", "ERB3P_ShP_G_GABP_GAP_Ras_GTP"];
            [~,null_species]=ismember(null_species_names, new_CMIM.species.names);
            const_species=[];
    end
    
    % 2. Define the mutated initial condition
    % lof: the removable conservation law and the species involved
    % lof_type2: the removable conservation law and species
    % gof: redefine
    switch protein_selected
        case lof_mutation
            x_0_mut=x_0_phys;
            basic_rem=find(new_CMIM.species.names==string(protein_selected));
            idx_law_rem=find(cons_laws(:, basic_rem));
            sp_rem=find(cons_laws(idx_law_rem,:));
            sp_rem=sort(sp_rem);
            x_0_mut(sp_rem)=[];
            idx_basic_mut=find(x_0_mut>0);
 
        case lof_mutation_type2
            x_0_mut=x_0_phys;
            basic_rem=find(new_CMIM.species.names=="TP53_generator");
            idx_law_rem=find(cons_laws(:, basic_rem));
            sp_rem_names=["TP53_generator", "TP53", "TP53U", "MDM2P_TP53", "TP53_TFBSIV"];
            [a, sp_rem]=ismember(sp_rem_names, new_CMIM.species.names);
            sp_rem=sort(sp_rem);
            x_0_mut(sp_rem)=[];
            idx_basic_mut=find(x_0_mut>0);
            
        case gof_mutation
            idx_law_rem=[];
            sp_rem=[];
            idx_basic_mut=idx_basic_species;
            x_0_mut=x_0_phys;
            
    end
    
    % 3. Define the new conservation laws matrix
    
    switch protein_selected
        case [lof_mutation, lof_mutation_type2]
            cons_laws_mut=cons_laws;
            cons_laws_mut(idx_law_rem,:)=[];
            cons_laws_mut(:, sp_rem)=[];
            n_species_mut=n_species-numel(sp_rem);
            ind_one_mut=n_species_mut+1;
            n_basic_mut=size(cons_laws_mut,1);
            
        case gof_mutation
            cons_laws_mut=cons_laws;
            ind_one_mut=n_species+1;
            n_basic_mut=n_cons_laws;
            n_species_mut=n_species;
    end
    
    % 4. Define the  new reactions set, S, v, k
    
    S_mut=Sm;
    v_mut=vm;
    k_mut=k_values;
    S_mut(sp_rem,:)=[];
    S_rem=Sm(sp_rem,:);
    
    switch protein_selected
        
        case [lof_mutation, lof_mutation_type2]
            
            react_rem=find(sum(abs(new_CMIM.matrix.S(sp_rem,:))));
           
        case gof_mutation
            [~, aux] = ...
                f_define_mutated_condition(protein_selected, ...
                                    x_0_phys, k_values, new_CMIM, 0);
            react_rem = find(aux==0);                   
            % Define GoF as in Sommariva et al. Scientific Reports 2021
    end
     
    S_mut(:,react_rem)=[];
    v_mut(react_rem,:)=[];
    k_mut(react_rem)=[];
    n_react_mut=n_reactions-numel(react_rem);
    
    % 4.1: redefine the vector v
    first_column=v_mut(:,1);
    first_reagent=v_mut(:,2);
    second_reagent=v_mut(:,3);
    
    switch protein_selected
        case [lof_mutation, lof_mutation_type2]
            idx_species_old=setdiff(1:n_species+1, sp_rem);
            idx_reactions=setdiff(1:n_reactions, react_rem);
            for ir=1:n_react_mut
                first_column(ir)=find(first_column(ir)==idx_reactions);
                first_reagent(ir)=find(first_reagent(ir)==idx_species_old);
                second_reagent(ir)=find(second_reagent(ir)==idx_species_old);
            end
            
            v_mut(:,1)=first_column;
            v_mut(:,2)=first_reagent;
            v_mut(:,3)=second_reagent;
            
        case gof_mutation
            
            idx_reactions=setdiff(1:n_reactions, react_rem);
            
            first_column=v_mut(:,1);
            
            for ir=1:n_react_mut
                first_column(ir)=find(first_column(ir)==idx_reactions);
            end
            
            v_mut(:,1)=first_column;
    end
    
    
    rho_mut=cons_laws_mut*x_0_mut;
    ris=f_NLPC_restart(x_0_mut, k_mut, S_mut, cons_laws_mut, rho_mut,...
        idx_basic_mut, v_mut, ind_one_mut, max_counter, 0);
    x_e_mut=ris.x;
   
    
    %% Sensitivity matrices for the mutated CRC-CRN
    x_e_mut(null_species)=0;
    
    % 1. Compute the analic jacobian of v
    jacobian_v = f_compute_analytic_jacobian_v(v_mut, n_species_mut, ind_one_mut);
    
    % 2. Evaluate the Jacobian in x of f
    
    eval_jac = f_evaluate_jacobian_neworder(k_mut, x_e_mut, ...
        S_mut, idx_basic_mut, jacobian_v, cons_laws_mut);
    
    eval_jac_c=f_evaluate_jacobian_c(n_species_mut, n_cons_laws, mut_lof);
    eval_jac_k=f_evaluate_jacobian_k(v_mut, x_e_mut, S_mut, idx_basic_mut);
    
    % 3. Calcolate Jac_k_xeq e Jac_c_eq
    Jac_k_xeq= -inv(eval_jac)*eval_jac_k;
    Jac_c_xeq= -inv(eval_jac)*eval_jac_c;
    
    % 3. Calcolate the senitivity matrices
    
    x_values=x_e_mut;
    x_values(null_species)=1;
    L_c=(1./x_values).*(Jac_c_xeq).*(rho_mut');
    L_k=(1./x_values).*(Jac_k_xeq).*(k_mut');
    
    %% c-SSI for the mutated CRC-CRN
    
    [U_c,D_c,V_c]=svd(L_c);
    e_PCA_c =((V_c(:,1:n_basic_mut).^2)*(diag(D_c).^2))/sum(diag(D_c).^2); %se presi con il quadrato: e_j=(L^T*L)_{jj}/sum_k(lambda_k)
    
   
    %% Comparison with the physiological CRC-CRN
    jacobian_v_phys = f_compute_analytic_jacobian_v(vm, n_species, idx_one);
    x_values=x_eq_phys;
    
    eval_jac_phys = f_evaluate_jacobian_neworder(k_values, x_values, ...
        Sm, idx_basic_species, jacobian_v_phys, cons_laws);
    
    mut_lof=0;
    eval_jac_c_phys=f_evaluate_jacobian_c(n_species, n_cons_laws,mut_lof);
    Jac_c_xeq_phys= -inv(eval_jac_phys)*eval_jac_c_phys;
    c=cons_laws*x_0_phys;
    L_c_phys=(1./x_values).*(Jac_c_xeq_phys).*(c');
    [U_c_phys,D_c_phys,V_c_phys]=svd(L_c_phys);
    e_PCA_c_phys =  ((V_c_phys(:,1:n_cons_laws).^2)*(diag(D_c_phys).^2))/sum(diag(D_c_phys).^2); %se presi con il quadrato: e_j=(L^T*L)_{jj}/sum_k(lambda_k)

    %% Graphs
    

    aux_PCA_c_phys=e_PCA_c_phys;
    aux_PCA_c=e_PCA_c;
    aux_PCA_c_phys(idx_law_rem)=[]; 
    [aux_PCA_c_phys_sort, order_phys]=sort(aux_PCA_c_phys, 'descend'); 
 
    
    subplot(2,2,im)
    
    semilogy(aux_PCA_c(order_phys), 'k', 'linewidth', 2.5)
    hold on
    semilogy(aux_PCA_c_phys_sort, 'r--', 'linewidth', 2)
    
    
    switch protein_selected
        case "Ras"
            aux_title = 'GoF Ras (a)';
        case "APC"
            aux_title = 'LoF APC (b)';
        case "SMAD4"
            aux_title = 'LoF SMAD4 (c)';
        case "TP53"
            aux_title = 'LoF TP53 (d)';
    end
    
    switch im
        case 1
            ylabel('e^c_j')
        case 3
            ylabel('e^c_j')
            xlabel('j')
        case 4
            xlabel('j')
    end
    
    title(aux_title)
    xlim([1 n_basic_mut])
    xticks(0:20:n_basic_mut)
    
    
    legend('Mutated', 'Phys', 'Location', 'northeast')
      
    
end
saveas(f_eff_mut_log_c, fullfile(folder_results, 'single_gene_mutations_c.png'))
