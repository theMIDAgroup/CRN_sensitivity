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

%% 1. Reading the MIM

Sm=MIM.matrix.S;
vm=MIM.matrix.v;
species_names=MIM.species.names;
cons_laws=MIM.matrix.Nl;
ind_one_mut=MIM.matrix.ind_one;
n_species=size(Sm,1);
n_cons_laws=size(cons_laws,1);
n_reactions=numel(k_values);

%% 2. Initialize initial values

react_rem_2=[];
react_rem=[];
const_species=[];

max_counter=500;

%% 3. Create the mutatated CRC
lof_mutation = {'APC', 'SMAD4'};
gof_mutation = {'Ras'};
lof_mutation_type2 = {'TP53'};

switch protein_selected
    case lof_mutation
        
        mut_lof=1;
        basic_rem=find(species_names==string(protein_selected));
        idx_law_rem=find(cons_laws(:, basic_rem));
        null_species=find(cons_laws(idx_law_rem,:));
        null_species=sort(null_species);
        
    case lof_mutation_type2
        
        mut_lof=1;
        basic_rem=find(species_names=="TP53_generator");
        idx_law_rem=find(cons_laws(:, basic_rem));
        null_species_names=["TP53_generator", "TP53", "TP53U", "MDM2P_TP53", "TP53_TFBSIV"];
        [~, null_species]=ismember(null_species_names, species_names);
        null_species=sort(null_species);
        
    case gof_mutation
        
        mut_lof=0;
        idx_law_rem=[];
        null_species_names=["RP_GAP_Ras_GTP", "RP_G_GABP_GAP_Ras_GTP", ...
            "RP_ShP_G_GABP_GAP_Ras_GTP", "ERBP_GAP_Ras_GTP", ...
            "ERBP_G_GABP_GAP_Ras_GTP", "ERBP_ShP_G_GABP_GAP_Ras_GTP", ...
            "ERB3P_GAP_Ras_GTP", ...
            "ERB3P_G_GABP_GAP_Ras_GTP", "ERB3P_ShP_G_GABP_GAP_Ras_GTP"];
        [~,null_species]=ismember(null_species_names, MIM.species.names);
        
end

S_mut=Sm;
v_mut=vm;
k_mut=k_values;

switch protein_selected
    
    case [lof_mutation, lof_mutation_type2]
        
        x_0_mut=x_0_phys;
        x_0_mut(basic_rem)=0;
        react_rem_2=find(sum(abs(Sm(null_species,:))));
        idx_basic_mut=idx_basic_species;
        n_basic_mut=numel(find(x_0_mut));

    case gof_mutation
        [~, ind_react_rem] = ...
            f_define_mutated_condition(protein_selected, ...
            x_0_phys, k_values, MIM, 0);
        react_rem = find(ind_react_rem==0);
        idx_basic_mut=idx_basic_species;
        x_0_mut=x_0_phys;
        n_basic_mut=numel(find(x_0_mut));

end

S_mut(:,react_rem)=[];
v_mut(react_rem,:)=[];
k_mut(react_rem)=[];
n_react_mut=n_reactions-numel(react_rem);

% 4.1: redefine the vector v
first_column=v_mut(:,1);
switch protein_selected
    
    case gof_mutation
        idx_reactions=setdiff(1:n_reactions, react_rem);

        for ir=1:n_react_mut
            first_column(ir)=find(first_column(ir)==idx_reactions);
        end
        
        v_mut(:,1)=first_column;
end

cons_laws_mut=f_compute_semipositive_conservations(S_mut);

rho_mut=cons_laws_mut*x_0_mut;

ris=f_NLPC_restart(x_0_mut, k_mut, S_mut, cons_laws_mut, rho_mut,...
    idx_basic_mut, v_mut, ind_one_mut, max_counter, 0);
x_e_mut=ris.x;


%% Step 4. Output 

MIM_mut.species.x_eq=x_e_mut;
MIM_mut.species.x_0=x_0_mut;
MIM_mut.species.idx_basic_species=idx_basic_mut;
MIM_mut.species.const_species=const_species;
MIM_mut.species.idx_law_rem=idx_law_rem;
MIM_mut.rates.std_values=k_mut;
MIM_mut.matrix.S=S_mut;
MIM_mut.matrix.v=v_mut;
MIM_mut.matrix.Nl=cons_laws_mut;
MIM_mut.matrix.ind_one=ind_one_mut;
MIM_mut.matrix.rho=rho_mut;
MIM_mut.info.name=protein_selected;
MIM_mut.info.mut_lof=mut_lof;
MIM_mut.info.react_rem=react_rem;
MIM_mut.info.react_rem_lof=react_rem_2;
MIM_mut.info.null_species=null_species;
MIM_mut.info.n_basic_mut=n_basic_mut;

if(protein_selected=="Ras")
    MIM_mut.species.null_species_names=null_species_names;
end





