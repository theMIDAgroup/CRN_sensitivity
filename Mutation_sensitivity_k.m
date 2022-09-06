%% Clear all
clc
clear
close all

%% Load data
% Load structure containing information on species and reactions
file_crn = fullfile('data', 'CRC_CRN_nodrug.mat'); % Da cambiare in base a dove metti i file
load(file_crn)
file_eq = fullfile('data', 'results_physiological.mat');
load(file_eq, 'ris_phys')
file_crc_crn = fullfile('data', 'CRC_CRN.mat');
load(file_crc_crn)

%% Initial data

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

%% Mutations
lof_mutation = {'APC', 'SMAD4'};
gof_mutation = {'Ras'};
lof_mutation_type2 = {'TP53'};
all_proteins = [gof_mutation, lof_mutation, lof_mutation_type2];

protein=["Ras", "APC", "SMAD4", "TP53"];


f_eff_mut_log_k = figure('units','normalized','outerposition',[0 0  0.75 1]);


for im=1:numel(protein)
    protein_selected=protein(im);
    
    % 1. Load results, specify the kind of mutation and the null species at the
    % equilibrium
    switch protein_selected
        
        case 'APC'
            ris_mutated= fullfile('results', 'results_mutation_APC_perc_0.0.mat');
            load(ris_mutated)
            mut_lof=1;
            null_species=[];
            const_species=[];
        case 'SMAD4'
            ris_mutated= fullfile('results', 'results_mutation_SMAD4_perc_0.0.mat');
            load(ris_mutated)
            mut_lof=1;
            null_species=[];
            const_species='TFBIS';
        case 'TP53'
            ris_mutated= fullfile('results', 'results_mutation_TP53_perc_0.0.mat');
            load(ris_mutated)
            mut_lof=1;
            null_species=[];
            const_species=[];
        case 'Ras'
            ris_mutated= fullfile('results', 'results_mutation_Ras_perc_0.0.mat');
            load(ris_mutated)
            mut_lof=0;
            null_species_names=["RP_GAP_Ras_GTP", "RP_G_GABP_GAP_Ras_GTP", ...
                "RP_ShP_G_GABP_GAP_Ras_GTP", "ERBP_GAP_Ras_GTP", ...
                "ERBP_G_GABP_GAP_Ras_GTP", "ERBP_ShP_G_GABP_GAP_Ras_GTP", ...
                "ERB3P_GAP_Ras_GTP", ...
                "ERB3P_G_GABP_GAP_Ras_GTP", "ERB3P_ShP_G_GABP_GAP_Ras_GTP"];
            [~,null_species]=ismember(null_species_names, CMIM.species.names);
            const_species=[];
    end
    
    % 2. Define the mutated initial condition
    % lof: the removable conservation law and the species involved
    % lof_type2: the removable conservation law and species
    % gof: redefine
    
    switch protein_selected
        case lof_mutation
            x_0_mut=x_0_phys;
            basic_rem=find(CMIM.species.names==string(protein_selected));
            idx_law_rem=find(cons_laws(:, basic_rem));
            sp_rem=find(cons_laws(idx_law_rem,:));
            sp_rem=sort(sp_rem);
            x_0_mut(sp_rem)=[];
            idx_basic_mut=find(x_0_mut>0);
            
        case lof_mutation_type2
            x_0_mut=x_0_phys;
            basic_rem=find(CMIM.species.names=="TP53_generator");
            idx_law_rem=find(cons_laws(:, basic_rem));
            sp_rem_names=["TP53_generator", "TP53", "TP53U", "MDM2P_TP53", "TP53_TFBSIV"];
            [a, sp_rem]=ismember(sp_rem_names, CMIM.species.names);
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
            
            
            react_rem=find(sum(abs(CMIM.matrix.S(sp_rem,:))));
            %         react_rem=zeros(n_reactions,1);
            %         for ir=1:n_reactions
            %             if sum(abs(S_rem(:, ir))) ~= 0
            %                 react_rem(ir)=ir;
            %             end
            %         end
            %         react_rem=react_rem(react_rem~=0);
            %         react_rem=sort(react_rem);
            
        case gof_mutation
            react_rem=sort(ris_mutated.index_rc);
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
    max_counter=300;
    
    x_e_mut=ris_mutated.x_e_mut;
    %     switch protein_selected
    %         case 'TP53'
    %             ris=f_PNG_restart(x_0_mut, k_mut, S_mut, cons_laws_mut, rho_mut,...
    %                 idx_basic_mut, v_mut, ind_one_mut, max_counter);
    %             x_e_mut=ris.x;
    %             ris_mutated.x_e_mut=x_e_mut;
    %         otherwise
    %
    %
    %     end
    
    
    %% Sensitivity matrices
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
    
    %% Sensitivity indices
    
    [U_c,D_c,V_c]=svd(L_c);
    e_PCA_c =((V_c(:,1:n_basic_mut).^2)*(diag(D_c).^2))/sum(diag(D_c).^2); %se presi con il quadrato: e_j=(L^T*L)_{jj}/sum_k(lambda_k)
    
    % figure
    % plot(e_PCA_c)
    % xlabel('c')
    % ylabel('e_j')
    % xlim([1 n_cons_laws-1])
    % title('Indici di sensitività rete con LoF di APC')
    
    [U_k,D_k,V_k]=svd(L_k);
    e_PCA_k =  ((V_k(:,1:n_species_mut).^2)*(diag(D_k).^2))/sum(diag(D_k).^2); %se presi con il quadrato: e_j=(L^T*L)_{jj}/sum_k(lambda_k)
    % figure
    % plot(e_PCA_k)
    % xlabel('k')
    % ylabel('e_j')
    % xlim([1 numel(k_mut)])
    % title('Indici di sensitività rete con LoF di APC')
    %
    
    %% Confronto con rete sana
    idx_one=new_CMIM.matrix.ind_one;
    jacobian_v_phys = f_compute_analytic_jacobian_v(vm, n_species, idx_one);
    x_values=ris_phys.x_eq(1:end-2);
    
    eval_jac_phys = f_evaluate_jacobian_neworder(k_values, x_values, ...
        Sm, idx_basic_species, jacobian_v_phys, cons_laws);
    
    mut_lof=0;
    eval_jac_c_phys=f_evaluate_jacobian_c(n_species, n_cons_laws,mut_lof);
    eval_jac_k_phys=f_evaluate_jacobian_k(vm, x_values, Sm, idx_basic_species);
    
    Jac_c_xeq_phys= -inv(eval_jac_phys)*eval_jac_c_phys;
    Jac_k_xeq_phys= -inv(eval_jac_phys)*eval_jac_k_phys;
    
    c=cons_laws*x_0_phys;
    L_c_phys=(1./x_values).*(Jac_c_xeq_phys).*(c');
    L_k_phys=(1./x_values).*(Jac_k_xeq_phys).*(k_values');
    
    [U_c_phys,D_c_phys,V_c_phys]=svd(L_c_phys);
    e_PCA_c_phys =  ((V_c_phys(:,1:n_cons_laws).^2)*(diag(D_c_phys).^2))/sum(diag(D_c_phys).^2); %se presi con il quadrato: e_j=(L^T*L)_{jj}/sum_k(lambda_k)
    
    [U_k_phys,D_k_phys,V_k_phys]=svd(L_k_phys);
    e_PCA_k_phys =  ((V_k_phys(:,1:n_species).^2)*(diag(D_k_phys).^2))/sum(diag(D_k_phys).^2); %se presi con il quadrato: e_j=(L^T*L)_{jj}/sum_k(lambda_k)
    
    
    
    %% Maggiori differenze
    aux_PCA_k=zeros(n_reactions,1);
    
    % creo vettore di e_PCA_k dove ho i valori di e_PCA e nelle entrate di
    % react rem metto 0.
    % in questo modo posso usare i dati (recations .arrow e .reactions2fluxes)
    
    for ik=1:n_reactions
        if ismember(ik, react_rem)==0
            b=0;
            if ik>react_rem(1)
                b=numel(find(ik>react_rem));
            end
            aux_PCA_k(ik)=e_PCA_k(ik-b);
        end
    end
    
    aux_PCA_k_p=e_PCA_k_phys;
    aux_PCA_k_p(react_rem)=0;
    
    diff_k=abs(aux_PCA_k_p-aux_PCA_k);
    index_max_k=zeros(n_reactions, 1);
    
    % confronto del massimo
    
%     ik=1;
    
    [max, index_max_k] = sort(diff_k, 'descend');
    
% se voglio solo quelle ali che hanno differenza >10^-2: 
% index_max_k = index_max_k(max > 10^-2);


% se non uso sort:
%     while ik<n_reactions+1
%         [m,i]=max(diff_k);
%         if numel(i)==1
%             index_max_k(ik)=i;
%             ik=ik+1;
%         else
%             index_max_k(ik:ik+numel(aux)-1)=i;
%             ik=ik+numel(aux);
%         end
%         diff_k(i)=-1;
%     end
    
    e_max_phys=aux_PCA_k_p(index_max_k);
    e_max=aux_PCA_k(index_max_k);
    
    
    
    table(string(CMIM.reactions.details(index_max_k,1)), e_max_phys-e_max, 'VariableNames', {'Reactions', 'Diff e_j'})
    
    table_file=fullfile('.\results', 'reactions_Ras.txt');
    fileID = fopen(table_file, 'w');
    disp(['Writing on ', table_file, '...'])
%     aux_idx = find(species.is_constant);
    

for ii = 1:numel(index_max_k)
        fprintf(fileID, '\\verb| %s | & %1.4f \\\\ \\hline \n', ...
        string(CMIM.reactions.details(index_max_k(ii),1)), ...
        e_max_phys(ii)-e_max(ii));  
end
    fclose(fileID);
    %% Grafici
    
    % elimino la legge di APC
    e_PCA_c_phys(idx_law_rem)=[];
    e_PCA_k_phys(react_rem)=[];
    
    % Sortiamo le entrate di e_PCA_c_phys e e_PCA_k_phys
    e_PCA_c_phys_sort=sort(e_PCA_c_phys, 'descend');
    e_PCA_k_phys_sort=sort(e_PCA_k_phys, 'descend');
    
    e_PCA_c_sort=zeros(n_basic_mut,1);
    e_PCA_k_sort=zeros(n_react_mut,1);
    
    ic=1;
    while ic<=n_basic_mut
        aux=e_PCA_c(find(e_PCA_c_phys==e_PCA_c_phys_sort(ic)));
        if numel(aux)==1
            e_PCA_c_sort(ic)=aux;
            ic=ic+1;
        else
            e_PCA_c_sort(ic:ic+numel(aux)-1)=aux;
            ic=ic+numel(aux);
        end
    end
    ik=1;
    while ik<=n_react_mut
        aux=e_PCA_k(find(e_PCA_k_phys==e_PCA_k_phys_sort(ik)));
        if numel(aux)==1
            e_PCA_k_sort(ik)=aux;
            ik=ik+1;
        else
            e_PCA_k_sort(ik:ik+numel(aux)-1)=aux;
            ik=ik+numel(aux);
        end
    end
    %
    % figure
    % semilogy(e_PCA_c_phys_sort)
    % hold on
    % semilogy(e_PCA_c_sort)
    % xlabel('c')
    % ylabel('e_j')
    % xlim([1 n_basic_mut])
    % legend('Physiological Network','Mutated Network')
    % title('Sensitivity map')
    
    
    
    
    subplot(2,2,im)
    
    semilogy(e_PCA_k_sort, 'k', 'linewidth', 2.5)
    hold on
    semilogy(e_PCA_k_phys_sort, 'r--', 'linewidth', 2)
    
    
    switch protein_selected
        case "Ras"
            aux_title = 'GoF KRAS (a)';
        case "APC"
            aux_title = 'LoF APC (b)';
        case "SMAD4"
            aux_title = 'LoF SMAD4 (c)';
        case "TP53"
            aux_title = 'LoF TP53 (d)';
    end
    
    switch im
        case 1
            ylabel('indici di sensitività')
            
        case 3
            ylabel('indici di sensitività')
            xlabel('k')
        case 4
            xlabel('k')
    end
    
    if im < 3
        set(gca,'xticklabel',{[]}, 'FontSize', 20)
    else
        set(gca, 'FontSize', 20)
    end
    title(aux_title)
    xlim([1 n_react_mut])
    
    legend('Mutated', 'Phys', 'Location', 'southwest')
    
    
end

saveas(f_eff_mut_log_k, fullfile('.\results', 'single_gene_mutations_k.png'))