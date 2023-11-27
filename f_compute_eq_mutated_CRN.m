function MIM_mut=f_compute_eq_mutated_CRN(protein_selected, MIM, idx_basic_species, x_0_phys, k_values)

% This function allows to compute the equilibrium point ofa mutated CRC-CRN affected by
% - Gain of function of KRAS
% - Loss of function of TP53
% - Loss of function of APC
% - Loss of function of SMAD4

% The inputs are
% - the selected protein (among KRAS, TP53, APC or SMAD4)
% - original MIM

% The output is the result computed by NLPC algorithm described in
% (add ref)

%% 1. reading the MIM
species_names=MIM.species.names;
cons_laws=MIM.matrix.Nl;
Sm=MIM.matrix.S;
vm=MIM.matrix.v;
n_species=size(Sm,1);
n_cons_laws=size(cons_laws,2);
n_reactions=numel(k_values);


max_counter=500;

%% Step 3. Create the mutatated CRC
lof_mutation = {'APC', 'SMAD4'};
gof_mutation = {'Ras'};
lof_mutation_type2 = {'TP53'};
all_proteins = [gof_mutation, lof_mutation, lof_mutation_type2];


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
        [~,null_species]=ismember(null_species_names, MIM.species.names);
        const_species=[];
end

% 2. Define the mutated initial condition
% lof: the removable conservation law and the species involved
% lof_type2: the removable conservation law and species
% gof: redefine
switch protein_selected
    case lof_mutation
        x_0_mut=x_0_phys;
        basic_rem=find(species_names==string(protein_selected));
        idx_law_rem=find(cons_laws(:, basic_rem));
        sp_rem=find(cons_laws(idx_law_rem,:));
        sp_rem=sort(sp_rem);
        x_0_mut(sp_rem)=[];
        idx_basic_mut=find(x_0_mut>0);
        
    case lof_mutation_type2
        x_0_mut=x_0_phys;
        basic_rem=find(species_names=="TP53_generator");
        idx_law_rem=find(cons_laws(:, basic_rem));
        sp_rem_names=["TP53_generator", "TP53", "TP53U", "MDM2P_TP53", "TP53_TFBSIV"];
        [~, sp_rem]=ismember(sp_rem_names, species_names);
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
        
        react_rem=find(sum(abs(MIM.matrix.S(sp_rem,:))));
        
    case gof_mutation
        [~, ind_react_rem] = ...
            f_define_mutated_condition(protein_selected, ...
            x_0_phys, k_values, MIM, 0);
        react_rem = find(ind_react_rem==0);
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

disp(size(cons_laws_mut))
disp(size(x_0_mut))
rho_mut=cons_laws_mut*x_0_mut;
disp(size(rho_mut))
ris=f_NLPC_restart(x_0_mut, k_mut, S_mut, cons_laws_mut, rho_mut,...
    idx_basic_mut, v_mut, ind_one_mut, max_counter, 0);
x_e_mut=ris.x;


%% Step 4. Output 

MIM_mut.species.x_eq=x_e_mut;
MIM_mut.species.x_0=x_0_mut;
MIM_mut.species.idx_basic_species=idx_basic_mut;
MIM_mut.species.null_species=null_species;
MIM_mut.species.const_species=const_species;
MIM_mut.species.idx_law_rem=idx_law_rem;
MIM_mut.species.n_basic_mut=n_basic_mut;
MIM_mut.rates.std_values=k_mut;
MIM_mut.matrix.S=S_mut;
MIM_mut.matrix.v=v_mut;
MIM_mut.matrix.Nl=cons_laws_mut;
MIM_mut.matrix.ind_one=ind_one_mut;
MIM_mut.matrix.rho=rho_mut;
MIM_mut.info.name=protein_selected;
MIM_mut.info.mut_lof=mut_lof;
MIM_mut.info.react_rem=react_rem;

if(protein_selected=="Ras")
    MIM_mut.species.null_species_names=null_species_names;
end





